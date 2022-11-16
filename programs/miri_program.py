import pandas as pd
import pathlib
import sys

JWST_UTILS_DIR = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(JWST_UTILS_DIR))

import mast.mast as mast
import pipeline.miri as miri
import badass.badass_miri as badass_miri
import badass.miri_consts as miri_consts

class MiriProgram:

    download_funcs = {
        'uncal': 'download_uncal',
        'cube': 'download_mast_final',
    }

    def __init__(self, attr_dict):
        self.__dict__.update(attr_dict)

        for attr in ['program_id', 'redshift', 'sci_obs', 'bkgd_obs', 'data_dir', 'product_name']:
            if not hasattr(self, attr) or getattr(self, attr) is None:
                raise Exception('Program input missing expected value: {attr}'.format(attr=attr))

        self.data_dir.mkdir(exist_ok=True, parents=True)
        if not hasattr(self, 'mast_cube_prefix'):
            self.mast_cube_prefix = None


    def download_uncal(self):
        cached_file = self.data_dir.joinpath('miri_{pid}_uncal.csv'.format(pid=self.program_id))
        if cached_file.exists():
            df = pd.read_csv(cached_file)
        else:
            df = mast.get_data_products(str(self.program_id), 'MIRI', 1, 'UNCAL').to_pandas()
            df.to_csv(cached_file, index=False)

        # Get only IFU data
        df = df[(df.obs_id.str.contains('mirifulong')) | (df.obs_id.str.contains('mirifushort'))]

        # Download science files
        sci_fmt = 'jw%05d%03d' % (self.program_id, self.sci_obs)
        sci_df = df[df.obs_id.str.contains(sci_fmt)]

        dl_data_dir = self.data_dir.joinpath('input', 'stage1', 'sci')
        dl_data_dir.mkdir(exist_ok=True, parents=True)

        for file_name in sci_df.productFilename.values:
            mast.download_file(file_name, dest=dl_data_dir.joinpath(file_name))

        # Download background files
        bkgd_fmt = 'jw%05d%03d' % (self.program_id, self.bkgd_obs)
        bkgd_df = df[df.obs_id.str.contains(bkgd_fmt)]

        dl_data_dir = self.data_dir.joinpath('input', 'stage1', 'bkgd')
        dl_data_dir.mkdir(exist_ok=True, parents=True)

        for file_name in bkgd_df.productFilename.values:
            mast.download_file(file_name, dest=dl_data_dir.joinpath(file_name))


    def download_mast_final():
        cached_file = self.data_dir.joinpath('miri_{pid}_s3d.csv'.format(pid=self.program_id))
        if cached_file.exists():
            df = pd.read_csv(cached_file)
        else:
            df = mast.get_data_products(str(self.program_id), 'MIRI', 3, 'S3D').to_pandas()
            df.to_csv(cached_file, index=False)

        sci_fmt = 'jw%05d-o%03d' % (self.program_id, self.sci_obs)
        sci_df = df[df.obs_id.str.contains(sci_fmt)]

        dl_data_dir = self.data_dir.joinpath('output', 'stage3', 'sci')
        dl_data_dir.mkdir(exist_ok=True, parents=True)

        for file_name in sci_df.productFilename.values:
            if not self.mast_cube_prefix:
                self.mast_cube_prefix = file_name.name.split('-shortmediumlong_s3d.fits')[0].split('_ch')[0]
            mast.download_file(file_name, dest=dl_data_dir.joinpath(file_name))


    def run_pipeline(self):
        # STAGE 1 SCI
        sci_input_dir = self.data_dir.joinpath('input', 'stage1', 'sci')
        sci_output_dir = self.data_dir.joinpath('output', 'stage1', 'sci')
        sci_output_dir.mkdir(exist_ok=True, parents=True)
        miri.run_stage1_all(sci_input_dir, sci_output_dir)

        # STAGE 1 BKGD
        bkgd_input_dir = self.data_dir.joinpath('input', 'stage1', 'bkgd')
        bkgd_output_dir = self.data_dir.joinpath('output', 'stage1', 'bkgd')
        bkgd_output_dir.mkdir(exist_ok=True, parents=True)
        miri.run_stage1_all(bkgd_input_dir, bkgd_output_dir)

        # STAGE 2 SCI
        sci_input_dir = self.data_dir.joinpath('input', 'stage2', 'sci')
        # Move the stage 1 output to stage 2 input
        sci_input_dir.mkdir(exist_ok=True, parents=True)
        for f in sci_output_dir.glob('*'):
          f.rename(sci_input_dir.joinpath(f.name))

        sci_output_dir = self.data_dir.joinpath('output', 'stage2', 'sci')
        sci_output_dir.mkdir(exist_ok=True, parents=True)
        miri.run_stage2_all(sci_input_dir, sci_output_dir, skip_cubes=True)

        # STAGE 2 BKGD
        bkgd_input_dir = self.data_dir.joinpath('input', 'stage2', 'bkgd')
        # Move the stage 1 output to stage 2 input
        bkgd_input_dir.mkdir(exist_ok=True, parents=True)
        for f in bkgd_output_dir.glob('*'):
          f.rename(bkgd_input_dir.joinpath(f.name))

        bkgd_output_dir = self.data_dir.joinpath('output', 'stage2', 'bkgd')
        bkgd_output_dir.mkdir(exist_ok=True, parents=True)
        miri.run_stage2_all(input_dir, bkgd_output_dir, skip_cubes=False)

        sci_input_dir = self.data_dir.joinpath('input', 'stage3', 'sci')
        # Move the stage 1 output to stage 2 input
        sci_input_dir.mkdir(exist_ok=True, parents=True)
        for f in sci_output_dir.glob('*'):
          f.rename(sci_input_dir.joinpath(f.name))

        bkgd_input_dir = self.data_dir.joinpath('input', 'stage3', 'bkgd')
        # Move the stage 1 output to stage 2 input
        bkgd_input_dir.mkdir(exist_ok=True, parents=True)
        for f in bkgd_output_dir.glob('*'):
          f.rename(bkgd_input_dir.joinpath(f.name))

        output_dir = self.data_dir.joinpath('output', 'stage3', 'sci')
        output_dir.mkdir(exist_ok=True, parents=True)
        miri.run_stage3(sci_input_dir, bkgd_input_dir, self.product_name, output_dir)

        print('Pipeline for {pname} complete!'.format(pname=self.product_name))
        self.pipeline_out_dir = output_dir


    def init_for_badass(self):
        if not hasattr(self, 'pipeline_out_dir'):
            raise Exception('Set pipeline output directory!')

        self.badass_dir = self.data_dir.joinpath('badass')
        self.badass_dir.mkdir(exist_ok=True, parents=True)

        cubes = list(self.pipeline_out_dir.glob('*-shortmediumlong_s3d.fits'))
        if len(cubes) == 0:
            print('WARNING: Found 0 prepared IFU cubes!')
            return

        print('Found {n} cubes'.format(n=len(cubes)))
        for cube in cubes:
            cube.rename(self.badass_dir.joinpath(cube.name))


    def run_badass_line_region(self, line_name, options_file):
        if not hasattr(self, 'badass_dir') or not self.badass_dir.exists():
            raise Exception('Badass directory not set or initialized!')

        line = badass_miri.get_line(line_name)

        cube_prefix = self.mast_cube_prefix if (self.product_name == 'mast_direct') else self.product_name
        cube_fits = self.badass_dir.joinpath('%s_ch%d-shortmediumlong_s3d.fits' % (cube_prefix, line.channel))
        if not cube_fits.exists():
            raise Exception('Could not find cube fits file: %s' % str(cube_fits))

        spaxel_dir = cube_fits.with_suffix('')
        if not spaxel_dir.exists():
            specres = badass_miri.get_channel_resolution(line.channel, line.subarray)
            badass_miri.prepare_cube(cube_fits, specres, z=self.redshift)

        badass_miri.run_badass(spaxel_dir, options_file, '%s_region'%line_name, fit_reg=(line.wave-2000,line.wave+2000))


    def run_cube_reconstruct(self, line_name):
        line = badass_miri.get_line(line_name)

        cube_prefix = self.mast_cube_prefix if (self.product_name == 'mast_direct') else self.product_name
        cube_fits = self.badass_dir.joinpath('%s_ch%d-shortmediumlong_s3d.fits' % (cube_prefix, line.channel))
        if not cube_fits.exists():
            raise Exception('Could not find cube fits file: %s' % str(cube_fits))

        badass_miri.reconstruct_cube(cube_fits, '%s_region'%line_name)
        badass_miri.plot_cube(cube_fits, '%s_region'%line_name)



