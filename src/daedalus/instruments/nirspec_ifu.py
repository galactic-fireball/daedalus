from astroquery.mast import Observations
from datetime import datetime, UTC
import json
import multiprocessing as mp

import daedalus.mast.mast as mast
from daedalus.instruments.instrument_common import Instrument

import jwst
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline
from jwst.pipeline.calwebb_spec3 import Spec3Pipeline


class NIRSpec_IFU(Instrument):

    def __init__(self, config):
        super().__init__(config)
        self.obs = self.config.target.instrument.obs


    def run_download(self, opts):
        if not opts.mast_token is None:
            mast.login(opts.mast_token)

        print('Downloading uncalibrated data to: %s'%str(self.config.uncal_dir))
        program_id = self.config.target.program

        df = mast.get_program_data(str(program_id), 'NIRSPEC/IFU')

        # uncal
        uncal_df = Observations.filter_products(df, calib_level=[1], productSubGroupDescription='UNCAL').to_pandas()
        obs_fmt = 'jw%05d%03d' % (program_id, self.obs)
        obs_df = uncal_df[uncal_df.obs_id.str.contains(obs_fmt)]

        for file_name in obs_df.productFilename.values:
            dest = self.config.uncal_dir.joinpath(file_name)
            if dest.exists():
                continue
            print('downloading %s' % file_name)
            mast.download_file(file_name, dest=dest)

        # asn2
        asn2_df = Observations.filter_products(df, calib_level=[2], productSubGroupDescription='ASN').to_pandas()
        obs_fmt = 'jw%05d%03d' % (program_id, self.obs)
        obs_df = asn2_df[(asn2_df.obs_id.str.contains(obs_fmt)) & (asn2_df.productFilename.str.contains('spec'))]

        for file_name in obs_df.productFilename.values:
            dest = self.config.stage2_dir.joinpath(file_name)
            if dest.exists():
                continue
            print('downloading %s' % file_name)
            mast.download_file(file_name, dest=dest)

        # asn3
        asn3_df = Observations.filter_products(df, calib_level=[3], productSubGroupDescription='ASN').to_pandas()
        obs_fmt = 'jw%05d-o%03d' % (program_id, self.obs)
        obs_df = asn3_df[asn3_df.obs_id.str.contains(obs_fmt)]

        for file_name in obs_df.productFilename.values:
            dest = self.config.stage3_dir.joinpath(file_name)
            if dest.exists():
                continue
            print('downloading %s' % file_name)
            mast.download_file(file_name, dest=dest)


    def update_asn(self, asn_temp):
        with open(asn_temp) as f:
            data = json.load(f)

        name = data['products'][0]['name']
        data['code_version'] = jwst.__version__
        data['version_id'] = datetime.now(UTC).strftime('%Y%m%dt%H%M%S')
        # data['asn_pool'] = ''

        asn_file = asn_temp.parent.joinpath('%s_updated.json' % asn_temp.stem)
        with open(asn_file, 'w') as f:
            json.dump(data, f)

        return asn_file, name

    def run_stage2_single(self, infile, opts):
        print('Processing: {}'.format(str(infile)))

        build_args = {
            'output_dir': str(self.config.stage3_dir),
            'save_results': True,
            'steps': opts.stage2.steps,
        }

        asn_in = infile.stem.split('_')[-1] == 'asn'
        if asn_in:
            infile, name = self.update_asn(infile)
            out_file = self.config.stage3_dir.joinpath('%s_cal.fits'%name)
            if out_file.exists() and not opts.overwrite:
                return None
        else:
            if self.config.stage3_dir.joinpath(infile.name.replace('rate', 'cal')).exists() and not opts.overwrite:
                return None


        config, _ = Spec2Pipeline.build_config(infile, **build_args)
        print('Running Spec 2 Pipeline with config: \n'+json.dumps(config.dict(), indent=4))
        Spec2Pipeline.call(infile, **build_args)


    def run_stage2(self, opts):

        # TODO: move all option parsing up to Instrument.__init__
        # As of JWST 1.20.2 / CRDS 13.1.10, this step is skipped by default, turn on here
        if opts.stage1.clean_noise: # clean_noise option is set in stage 1
            # This step is still called 'nsclean' as of JWST 1.20.2, but a change seems to be coming...
            if not 'nsclean' in opts.stage2.steps: opts.stage2.steps['nsclean'] = {}
            opts.stage2.steps['nsclean']['skip'] = False

            # if not 'clean_flicker_noise' in opts.stage2.steps: opts.stage2.steps['clean_flicker_noise'] = {}
            # opts.stage2.steps['clean_flicker_noise']['skip'] = False

        input_dir = self.config.stage2_dir
        in_files = list(input_dir.glob('*_spec2_*_asn.json'))
        if len(in_files) == 0:
            print('No stage 2 association files, running directly on rate files')
            in_files = list(input_dir.glob('*_rate.fits'))

        if not opts.multiprocess:
            for in_file in in_files:
                self.run_stage2_single(in_file, opts)
            return

        proc_args = [(in_file, opts) for in_file in in_files]
        pool = mp.Pool(processes=opts.nprocesses, maxtasksperchild=1)
        pool.starmap(self.run_stage2_single, proc_args, chunksize=1)
        pool.close()
        pool.join()


    def run_stage3_single(self, asn_temp, opts):
        print('Processing: {}'.format(str(asn_temp)))

        asn_file, name = self.update_asn(asn_temp)

        build_args = {
            'output_dir': str(self.config.output_dir),
            'save_results': True,
            'steps': opts.stage3.steps,
        }

        config, _ = Spec3Pipeline.build_config(asn_file, **build_args)
        print('Running Spec 3 Pipeline with config: \n'+json.dumps(config.dict(), indent=4))
        Spec3Pipeline.call(asn_file, **build_args)


    def run_stage3(self, opts):
        asn_files = list(self.config.stage3_dir.glob('*_spec3_*_asn.json'))
        if not opts.multiprocess:
            for asn_file in asn_files:
                self.run_stage3_single(asn_file, opts)
            return

        proc_args = [(asn_file, opts) for asn_file in asn_files]
        pool = mp.Pool(processes=opts.nprocesses, maxtasksperchild=1)
        pool.starmap(self.run_stage3_single, proc_args, chunksize=1)
        pool.close()
        pool.join()

