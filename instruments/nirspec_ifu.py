import astropy.io.fits as fits
import astropy.units as u
from astroquery.mast import Observations
from datetime import datetime
import json
import multiprocessing as mp
import numpy as np

import mast.mast as mast
from instruments.instrument_common import Instrument

import jwst
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline
from jwst.pipeline.calwebb_spec3 import Spec3Pipeline


class NIRSpec_IFU(Instrument):

    def download(self, context, args):
        if 'mast_token' in args:
            mast.login(args['mast_token'])

        program_id = context.target.program_id
        output_dir = args.get('output_dir', context.pipeline_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        df = mast.get_program_data(str(program_id), 'NIRSPEC/IFU')

        # uncal
        uncal_df = Observations.filter_products(df, calib_level=[1], productSubGroupDescription='UNCAL').to_pandas()
        obs_fmt = 'jw%05d%03d' % (program_id, self.obs)
        obs_df = uncal_df[uncal_df.obs_id.str.contains(obs_fmt)]

        for file_name in obs_df.productFilename.values:
            dest = output_dir.joinpath(file_name)
            if dest.exists():
                continue
            print('downloading %s' % file_name)
            mast.download_file(file_name, dest=dest)

        # asn2
        asn2_df = Observations.filter_products(df, calib_level=[2], productSubGroupDescription='ASN').to_pandas()
        obs_fmt = 'jw%05d%03d' % (program_id, self.obs)
        obs_df = asn2_df[(asn2_df.obs_id.str.contains(obs_fmt)) & (asn2_df.productFilename.str.contains('spec'))]

        for file_name in obs_df.productFilename.values:
            dest = output_dir.joinpath(file_name)
            if dest.exists():
                continue
            print('downloading %s' % file_name)
            mast.download_file(file_name, dest=dest)

        # asn3
        asn3_df = Observations.filter_products(df, calib_level=[3], productSubGroupDescription='ASN').to_pandas()
        obs_fmt = 'jw%05d-o%03d' % (program_id, self.obs)
        obs_df = asn3_df[asn3_df.obs_id.str.contains(obs_fmt)]

        for file_name in obs_df.productFilename.values:
            dest = output_dir.joinpath(file_name)
            if dest.exists():
                continue
            print('downloading %s' % file_name)
            mast.download_file(file_name, dest=dest)


    def update_asn(self, asn_temp):
        data = json.load(open(asn_temp))

        name = data['products'][0]['name']
        data['code_version'] = jwst.__version__
        data['version_id'] = datetime.utcnow().strftime('%Y%m%dt%H%M%S')
        # data['asn_pool'] = ''

        asn_file = asn_temp.parent.joinpath('%s_updated.json' % asn_temp.stem)
        json.dump(data, open(asn_file, 'w'))

        return asn_file, name

    def run_stage2_single(self, infile, output_dir, context, args):
        print('Processing: {}'.format(str(infile)))

        stage_args = args.get('stage2', {})
        overwrite = stage_args.get('overwrite', False)
        step_opts = stage_args.get('steps', {})

        asn_in = infile.stem.split('_')[-1] == 'asn'

        if asn_in:
            asn_file, name = self.update_asn(infile)
            out_file = output_dir.joinpath('%s_cal.fits'%name)
            if out_file.exists() and not overwrite:
                return None
        else:
            if output_dir.joinpath(infile.name.replace('rate', 'cal')).exists() and not overwrite:
                return None

        spec2 = Spec2Pipeline()
        spec2.output_dir = str(output_dir)
        for step, options in step_opts.items():
            for opt, val in options.items():
                setattr(getattr(spec2, step), opt, val)

        if asn_in:
            spec2.save_results = True
            result = spec2(str(asn_file))
        else: # Run stage 2 directly on rate files
            spec2.output_file = str(infile.with_suffix(''))
            spec2.run(infile)
            spec2.suffix = 'spec2'
        return None


    def run_stage2_all(self, context, args):
        stage2_opts = args.get('stage2', {})
        input_dir = stage2_opts.get('input_dir', context.pipeline_dir)
        output_dir = stage2_opts.get('output_dir', context.pipeline_dir)
        multiprocess = args.get('multiprocess', False)
        nprocesses = args.get('nprocesses', 1)

        in_files = input_dir.glob('*_spec2_*_asn.json')
        if len(list(in_files)) == 0:
            print('No stage 2 association files, running directly on rate files')
            in_files = input_dir.glob('*_rate.fits')

        if not multiprocess:
            for in_file in in_files:
                self.run_stage2_single(in_file, output_dir, context, args)
            return

        proc_args = [(in_file, output_dir, context, args) for in_file in in_files]
        pool = mp.Pool(processes=nprocesses, maxtasksperchild=1)
        pool.starmap(self.run_stage2_single, proc_args, chunksize=1)
        pool.close()
        pool.join()


    def run_stage3_single(self, asn_temp, out_dir, context, args):
        print('Processing: {}'.format(str(asn_temp)))
        stage_args = args.get('stage3', {})
        step_opts = stage_args.get('steps', {})

        asn_file, name = self.update_asn(asn_temp)
        crds_config = Spec3Pipeline.get_config_from_reference(str(asn_file))

        spec3 = Spec3Pipeline.from_config_section(crds_config)
        spec3.output_dir = str(out_dir)
        spec3.save_results = True

        for step, options in step_opts.items():
            for opt, val in options.items():
                setattr(getattr(spec3, step), opt, val)

        spec3.run(asn_file)
        return None


    def run_stage3_all(self, context, args):
        stage3_opts = args.get('stage3', {})
        input_dir = stage3_opts.get('input_dir', context.pipeline_dir)
        output_dir = stage3_opts.get('output_dir', context.pipeline_dir)
        multiprocess = args.get('multiprocess', False)
        nprocesses = args.get('nprocesses', 1)

        asn_files = input_dir.glob('*_spec3_*_asn.json')

        if not multiprocess:
            for asn_file in asn_files:
                self.run_stage3_single(asn_file, output_dir, context, args)
            return

        proc_args = [(asn_file, output_dir, context, args) for asn_file in asn_files]
        pool = mp.Pool(processes=nprocesses, maxtasksperchild=1)
        pool.starmap(self.run_stage3_single, proc_args, chunksize=1)
        pool.close()
        pool.join()

