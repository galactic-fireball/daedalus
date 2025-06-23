from astropy.io import fits
import copy
import multiprocessing as mp
import pandas as pd

import mast.mast as mast
from instruments.instrument_common import Instrument

from jwst.associations import asn_from_list as afl
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase
from jwst.associations import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline
from jwst.pipeline.calwebb_spec3 import Spec3Pipeline


class MIRI_IFU(Instrument):

    def post_setup(self, context):
        self.pipeline_sci_dir = context.pipeline_dir.joinpath('sci')
        self.pipeline_sci_dir.mkdir(parents=True, exist_ok=True)
        self.pipeline_bkgd_dir = context.pipeline_dir.joinpath('bkgd')
        self.pipeline_bkgd_dir.mkdir(parents=True, exist_ok=True)


    def download(self, context, args):
        # cached_file = self.data_dir.joinpath('miri_{pid}_uncal.csv'.format(pid=self.program_id))
        # if cached_file.exists():
        #     df = pd.read_csv(cached_file)
        # else:
            
        #     df.to_csv(cached_file, index=False)

        if 'mast_token' in args:
            mast.login(args['mast_token'])

        df = mast.get_data_products(str(context.target.program_id), 'MIRI/IFU', 1, 'UNCAL').to_pandas()

        # Get only IFU data
        df = df[(df.obs_id.str.contains('mirifulong')) | (df.obs_id.str.contains('mirifushort'))]

        # Download science files
        sci_fmt = 'jw%05d%03d' % (context.target.program_id, self.sci_obs)
        sci_df = df[df.obs_id.str.contains(sci_fmt)]

        for file_name in sci_df.productFilename.values:
            dest = self.pipeline_sci_dir.joinpath(file_name)
            if dest.exists():
                continue
            print('downloading %s' % file_name)
            mast.download_file(file_name, dest=dest)

        # Download background files
        bkgd_fmt = 'jw%05d%03d' % (context.target.program_id, self.bkgd_obs)
        bkgd_df = df[df.obs_id.str.contains(bkgd_fmt)]

        for file_name in bkgd_df.productFilename.values:
            dest = self.pipeline_bkgd_dir.joinpath(file_name)
            if dest.exists():
                continue
            print('downloading %s' % file_name)
            mast.download_file(file_name, dest=dest)


    # TODO: fix and add as download option
    def download_mast_final():
        cached_file = self.data_dir.joinpath('miri_{pid}_s3d.csv'.format(pid=self.program_id))
        if cached_file.exists():
            df = pd.read_csv(cached_file)
        else:
            df = mast.get_data_products(str(self.program_id), 'MIRI', 3, 'S3D').to_pandas()
            df.to_csv(cached_file, index=False)

        sci_fmt = 'jw%05d-o%03d' % (self.program_id, self.sci_obs)
        sci_df = df[df.obs_id.str.contains(sci_fmt)]

        self.mast_dir = self.data_dir.joinpath('mast')
        self.mast_dir.mkdir(parents=True, exist_ok=True)

        for file_name in sci_df.productFilename.values:
            if not self.mast_cube_prefix:
                self.mast_cube_prefix = file_name.name.split('-shortmediumlong_s3d.fits')[0].split('_ch')[0]
            mast.download_file(file_name, dest=self.mast_dir.joinpath(file_name))


    def run_pipeline(self, context, args):
        # crds_utils.cache_crds('miri')
        self.run_stage1_all(context, args, indir=self.pipeline_sci_dir, outdir=self.pipeline_sci_dir)
        self.run_stage1_all(context, args, indir=self.pipeline_bkgd_dir, outdir=self.pipeline_bkgd_dir)
        self.run_stage2_all(context, args, self.pipeline_sci_dir, self.pipeline_bkgd_dir, self.pipeline_sci_dir, self.pipeline_bkgd_dir, skip_cubes=True)
        self.run_stage3(context, args, self.pipeline_sci_dir, self.pipeline_bkgd_dir, context.pipeline_dir)
        print('Pipeline for {pname} complete!'.format(pname=context.config.product_name))


    def run_stage2_single_orig(self, context, args, rfile, output_dir, skip_cubes):
        print('Processing: {}'.format(str(rfile)))

        stage_args = args.get('stage2', {})
        overwrite = stage_args.get('overwrite', False)
        step_opts = stage_args.get('steps', {})

        if output_dir.joinpath(rfile.name.replace('rate', 'cal')).exists() and not overwrite:
            return None

        spec2 = Spec2Pipeline()
        spec2.output_dir = str(output_dir)
        spec2.output_file = str(rfile.with_suffix(''))

        if (skip_cubes):
            spec2.cube_build.skip = True
            spec2.extract_1d.skip = True

        if stage_args.get('fringe_corr', False):
            spec2.residual_fringe.skip = False

        for step, options in step_opts.items():
            for opt, val in options.items():
                setattr(getattr(spec2, step), opt, val)

        spec2.run(rfile)
        spec2.suffix = 'spec2'
        return None


    def run_stage2_all_orig(self, context, args, input_dir, output_dir, skip_cubes=False):
        rate_files = input_dir.glob('*_rate.fits')
        multiprocess = args.get('multiprocess', False)
        nprocesses = args.get('nprocesses', 1)

        if not multiprocess:
            for rfile in rate_files:
                self.run_stage2_single(context, args, rfile, output_dir, skip_cubes)
            return

        args = [(context, args, rfile, output_dir, skip_cubes) for rfile in rate_files]
        pool = mp.Pool(processes=nprocesses, maxtasksperchild=1)
        pool.starmap(self.run_stage2_single, args, chunksize=1)
        pool.close()
        pool.join()


    # The following functions for stage 2 and 3 were modified from David Law and Kirsten Larson's jupyter notebook:
    #   https://github.com/spacetelescope/jwst-pipeline-notebooks/blob/main/notebooks/MIRI/MRS/JWPipeNB-MIRI-MRS.ipynb
    def create_stage2_asn(self, context, rfile, bg_files, cal_files):
        asn = afl.asn_from_list([str(rfile),], rule=DMSLevel2bBase, product_name=context.config.product_name+'_level2')

        hdu = fits.open(rfile)
        hdr = hdu[0].header
        ch, band = hdr['CHANNEL'], hdr['BAND']
        hdu.close()

        for file in bg_files:
            hdu = fits.open(file)
            if (hdu[0].header['CHANNEL'] == ch) and (hdu[0].header['BAND'] == band):
                asn['products'][0]['members'].append({'expname':str(file), 'exptype':'background'})
            hdu.close()

        for file in cal_files:
            hdu = fits.open(file)
            if hdu[0].header['CHANNEL'] == ch:
                asn['products'][0]['members'].append({'expname':str(file), 'exptype':'selfcal'})
            hdu.close()

        asn_out = rfile.parent.joinpath('l2asn.json')

        _,asn_serialized = asn.dump()
        asn_file = open(asn_out, 'w')
        asn_file.write(asn_serialized)
        asn_file.close()
        return asn_out


    def run_stage2_single(self, context, args, step_dict, rfile, bg_rate_files, cal_files, output_dir):
        asnfile = str(self.create_stage2_asn(context, rfile, bg_rate_files, cal_files))
        Spec2Pipeline.call(asnfile, steps=step_dict, save_results=True, output_dir=str(output_dir))


    def run_stage2_all(self, context, args, sci_input_dir, bg_input_dir, sci_output_dir, bg_output_dir, skip_cubes=False):
        sci_rate_files = list(sci_input_dir.glob('*_rate.fits'))
        bg_rate_files = list(bg_input_dir.glob('*_rate.fits'))

        cal_files = sci_rate_files.copy() + bg_rate_files.copy()

        spec2_steps = [
            'assign_wcs', 'badpix_selfcal', 'bkg_subtract', 'flat_field', 'srctype', 'straylight',
            'fringe', 'photom', 'residual_fringe', 'pixel_replace', 'cube_build', 'extract_1d'
        ]
        spec2_dict = {s: {} for s in spec2_steps}

        stage_args = args.get('stage2', {})
        if stage_args.get('pixel_bg', False):
            spec2_dict['bkg_subtract']['skip'] = False
        else:
            spec2_dict['bkg_subtract']['skip'] = True

        spec2_sci_dict = copy.deepcopy(spec2_dict)
        spec2_sci_dict['cube_build']['skip'] = True
        spec2_sci_dict['extract_1d']['skip'] = True

        for rfile in sci_rate_files:
            self.run_stage2_single(context, args, spec2_sci_dict, rfile, bg_rate_files, cal_files, sci_output_dir)

        for rfile in bg_rate_files:
            self.run_stage2_single(context, args, spec2_dict, rfile, [], cal_files, bg_output_dir)


    def create_stage3_asn(self, context, sci_files, bg_files, output_dir):
        asn = afl.asn_from_list(sci_files, rule=DMS_Level3_Base, product_name=context.target.name)#context.config.product_name+'_level3')

        for file in bg_files:
            asn['products'][0]['members'].append({'expname':file, 'exptype':'background'})

        asn_out = output_dir.joinpath('l3asn.json')

        _,asn_serialized = asn.dump()
        asn_file = open(asn_out, 'w')
        asn_file.write(asn_serialized)
        asn_file.close()
        return asn_out


    def run_stage3(self, context, args, sci_input, bk_input, output_dir):
        spec3_steps = [
            'assign_mtwcs', 'master_background', 'outlier_detection', 'mrs_imatch',
            'cube_build', 'pixel_replace', 'extract_1d', 'spectral_leak'
        ]

        spec3_dict = {s:{} for s in spec3_steps}

        stage3_opts = args.get('stage3', {})

        if stage3_opts.get('master_bg', True):
            spec3_dict['master_background']['skip'] = False
        else:
            spec3_dict['master_background']['skip'] = True

        spec3_dict['pixel_replace']['skip'] = False
        spec3_dict['pixel_replace']['algorithm'] = 'mingrad'
        spec3_dict['cube_build']['output_type'] = 'channel' # 'channel', 'band', 'multi'
        # spec3_dict['cube_build']['list_par1'] = ['1', '1', '1', '2', '2', '2', '3', '3', '3', '4', '4', '4',]
        # spec3_dict['cube_build']['list_par1'] = ['short', 'medium', 'long', 'short', 'medium', 'long', 'short', 'medium', 'long', 'short', 'medium', 'long',]
        spec3_dict['extract_1d']['ifu_autocen'] = True

        sci_cal_files = [str(s) for s in sci_input.glob('*_cal.fits')]
        bg_x1d_files = [str(s) for s in bk_input.glob('*_x1d.fits')]

        asn_file = self.create_stage3_asn(context, sci_cal_files, bg_x1d_files, output_dir)
        Spec3Pipeline.call(asn_file, steps=spec3_dict, save_results=True, output_dir=str(output_dir))

