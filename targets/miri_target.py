import pandas as pd
import pathlib
import re
import sys

import mast.mast as mast
import pipeline.miri as miri
import badass.badass_miri as badass_miri
import badass.miri_consts as miri_consts
import pipeline.crds_utils as crds_utils
from targets.target_common import Target

class MiriTarget(Target):

    def init(self):
        super().init()

        for attr in ['redshift', 'sci_obs', 'bkgd_obs']:
            if not hasattr(self, attr) or getattr(self, attr) is None:
                raise Exception('Program input missing expected value: {attr}'.format(attr=attr))

        self.pipeline_sci_dir = self.pipeline_dir.joinpath('sci')
        self.pipeline_sci_dir.mkdir(parents=True, exist_ok=True)
        self.pipeline_bkgd_dir = self.pipeline_dir.joinpath('bkgd')
        self.pipeline_bkgd_dir.mkdir(parents=True, exist_ok=True)


    def download(self):
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

        for file_name in sci_df.productFilename.values:
            print('downloading %s' % file_name)
            mast.download_file(file_name, dest=self.pipeline_sci_dir.joinpath(file_name))

        # Download background files
        bkgd_fmt = 'jw%05d%03d' % (self.program_id, self.bkgd_obs)
        bkgd_df = df[df.obs_id.str.contains(bkgd_fmt)]

        for file_name in bkgd_df.productFilename.values:
            print('downloading %s' % file_name)
            mast.download_file(file_name, dest=self.pipeline_bkgd_dir.joinpath(file_name))


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


    def pipeline(self):
        crds_utils.cache_crds('miri')

        miri.run_stage1_all(self.pipeline_sci_dir, self.pipeline_sci_dir)
        miri.run_stage1_all(self.pipeline_bkgd_dir, self.pipeline_bkgd_dir)
        miri.run_stage2_all(self.pipeline_sci_dir, self.pipeline_sci_dir, skip_cubes=True)
        miri.run_stage2_all(self.pipeline_bkgd_dir, self.pipeline_bkgd_dir, skip_cubes=False)

        miri.run_stage3(self.pipeline_sci_dir, self.pipeline_bkgd_dir, self.product_name, self.pipeline_dir)

        print('Pipeline for {pname} complete!'.format(pname=self.product_name))


    def get_cubes(self):
        cubes = list(self.pipeline_dir.glob('*-shortmediumlong_s3d.fits'))
        print('Found {n} cubes'.format(n=len(cubes)))
        for cube in cubes:
            cube.rename(self.badass_dir.joinpath(cube.name))

        # grab the list this way since we might not have actually moved anything
        return list(self.badass_dir.glob('*_s3d.fits'))


    def extract_1D(self):
        cubes = self.get_cubes()
        for cube in cubes:
            aperture_dir = cube.with_suffix('').joinpath('aperture')
            aperture_dir.mkdir(parents=True, exist_ok=True)


    def extract_spaxels(self):
        cubes = self.get_cubes()

        for cube in cubes:
            spaxel_dir = cube.with_suffix('')
            if not spaxel_dir.exists():
                channel = int(re.findall(r'_ch\d+-shortmediumlong_s3d', cube.name)[0].split('-')[0][3:])
                specres = badass_miri.get_channel_resolution(channel, 'short')
                badass_miri.prepare_cube(cube, specres, z=self.redshift)


    # def run_badass_line_region_extracted(self, line_name, options_file):
    #     self.badass_dir = self.data_dir.joinpath('badass')
    #     extracted_dir = self.badass_dir.joinpath('extracted')
    #     if not extracted_dir.exists():
    #         raise Exception('Extracted directory not found')

    #     line = badass_miri.get_line(line_name)

    #     cube_prefix = self.mast_cube_prefix if (self.product_name == 'mast_direct') else self.product_name
    #     spec_fits = extracted_dir.joinpath('%s_ch%d_extract.fits' % (cube_prefix, line.channel))
    #     if not spec_fits.exists():
    #         raise Exception('Could not find extracted spectrum file: %s' % str(spec_fits))

    #     specres = badass_miri.get_channel_resolution(line.channel, line.subarray)

    #     # TODO: make wavelength range an option
    #     badass_miri.run_badass_extracted(spec_fits, options_file, '%s_region'%line_name, specres, self.redshift, fit_reg=(line.wave-2000,line.wave+2000))


    # def run_badass_line_region(self, line_name, options_file):
    #     if not hasattr(self, 'badass_dir') or not self.badass_dir.exists():
    #         raise Exception('Badass directory not set or initialized!')

    #     line = badass_miri.get_line(line_name)

    #     cube_prefix = self.mast_cube_prefix if (self.product_name == 'mast_direct') else self.product_name
    #     cube_fits = self.badass_dir.joinpath('%s_ch%d-shortmediumlong_s3d.fits' % (cube_prefix, line.channel))
    #     if not cube_fits.exists():
    #         raise Exception('Could not find cube fits file: %s' % str(cube_fits))

    #     spaxel_dir = cube_fits.with_suffix('')
    #     if not spaxel_dir.exists():
    #         raise Exception('Spaxel directories not initialized (%s)' % (str(spaxel_dir)))

    #     badass_miri.run_badass(spaxel_dir, options_file, '%s_region'%line_name, fit_reg=(line.wave-2000,line.wave+2000))


    # def run_cube_reconstruct(self, line_name):
    #     line = badass_miri.get_line(line_name)

    #     cube_prefix = self.mast_cube_prefix if (self.product_name == 'mast_direct') else self.product_name
    #     cube_fits = self.badass_dir.joinpath('%s_ch%d-shortmediumlong_s3d.fits' % (cube_prefix, line.channel))
    #     if not cube_fits.exists():
    #         raise Exception('Could not find cube fits file: %s' % str(cube_fits))

    #     badass_miri.reconstruct_cube(cube_fits, '%s_region'%line_name)
    #     badass_miri.plot_cube(cube_fits, '%s_region'%line_name)


    # def create_ratio_map(self, line1, line2):
    #     cube_prefix = self.mast_cube_prefix if (self.product_name == 'mast_direct') else self.product_name
    #     cube_fits1 = self.badass_dir.joinpath('%s_ch%d-shortmediumlong_s3d.fits' % (cube_prefix, line1.channel))
    #     cube1 = self.badass_dir.joinpath(cube_fits1.stem, 'full_cube', '%s_region'%line1.name)
    #     if not cube1.exists():
    #         raise Exception('Could not find spaxel directory: %s' % str(cube1))

    #     cube_fits2 = self.badass_dir.joinpath('%s_ch%d-shortmediumlong_s3d.fits' % (cube_prefix, line2.channel))
    #     cube2 = self.badass_dir.joinpath(cube_fits2.stem, 'full_cube', '%s_region'%line2.name)
    #     if not cube2.exists():
    #         raise Exception('Could not find spaxel directory: %s' % str(cube2))

    #     plot_dir = self.badass_dir.joinpath('ratio_plots')
    #     plot_dir.mkdir(exist_ok=True, parents=True)
    #     out_plot = plot_dir.joinpath('%s_vs_%s.pdf' % (line1.name, line2.name))

    #     badass_miri.plot_ratio(cube1, 'NA_%s_FLUX'%line1.name.upper(), cube2, 'NA_%s_FLUX'%line2.name.upper(), out_plot)


    # def create_all_ratio_maps(self):
    #     for i, line1 in enumerate(miri_consts.lines):
    #         for line2 in miri_consts.lines[i+1:]:
    #             self.create_ratio_map(line1, line2)


