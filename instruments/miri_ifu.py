from astropy import units as u
from astropy.io import fits
import multiprocessing as mp
import numpy as np
import pandas as pd
import pathlib
import re
import sys

import mast.mast as mast
import badass.badass_miri as badass_miri
import badass.miri_consts as miri_consts
import pipeline.crds_utils as crds_utils
from instruments.instrument_common import Instrument, Pipeline, Extractor

from jwst.associations import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from jwst.associations.lib.rules_level3 import Asn_Lv3SpectralSource
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline
from jwst.pipeline.calwebb_spec3 import Spec3Pipeline

class MIRI_IFU(Instrument):

    def init(self):
        super().init()

        for attr in ['redshift', 'sci_obs', 'bkgd_obs']:
            if not hasattr(self, attr) or getattr(self, attr) is None:
                raise Exception('Program input missing expected value: {attr}'.format(attr=attr))

        self.pipeline_sci_dir = self.pipeline_dir.joinpath('sci')
        self.pipeline_sci_dir.mkdir(parents=True, exist_ok=True)
        self.pipeline_bkgd_dir = self.pipeline_dir.joinpath('bkgd')
        self.pipeline_bkgd_dir.mkdir(parents=True, exist_ok=True)


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

        for file_name in sci_df.productFilename.values:
            dest = self.pipeline_sci_dir.joinpath(file_name)
            if dest.exists():
                continue
            print('downloading %s' % file_name)
            mast.download_file(file_name, dest=dest)

        # Download background files
        bkgd_fmt = 'jw%05d%03d' % (self.program_id, self.bkgd_obs)
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


    def run_pipeline(self):
        crds_utils.cache_crds('miri')
        pipeline_config = getattr(self, 'pipeline', {})
        pipeline = MIRI_IFU_Pipeline(pipeline_config)

        pipeline.run_stage1_all(self.pipeline_sci_dir, self.pipeline_sci_dir)
        pipeline.run_stage1_all(self.pipeline_bkgd_dir, self.pipeline_bkgd_dir)
        pipeline.run_stage2_all(self.pipeline_sci_dir, self.pipeline_sci_dir, skip_cubes=True)
        pipeline.run_stage2_all(self.pipeline_bkgd_dir, self.pipeline_bkgd_dir, skip_cubes=False)

        pipeline.run_stage3(self.pipeline_sci_dir, self.pipeline_bkgd_dir, self.product_name, self.pipeline_dir)

        print('Pipeline for {pname} complete!'.format(pname=self.product_name))


    def get_cubes(self):
        cubes = list(self.pipeline_dir.glob('*-shortmediumlong_s3d.fits'))
        print('Found {n} cubes'.format(n=len(cubes)))
        for cube in cubes:
            cube.rename(self.badass_dir.joinpath(cube.name))

        # grab the list this way since we might not have actually moved anything
        return list(self.badass_dir.glob('*_s3d.fits'))


    def extract_1D(self):
        # TODO: use default config
        # self.extract = getattr(self, 'extract', {})
        # self.extract['aperture_name'] = getattr(self.extract, 'aperture_name', 'aperture')
        # self.extract['aperture_radius'] = getattr(self.extract, 'aperture_radius', 'psf')

        cubes = self.get_cubes()
        for cube in cubes:
            aperture_dir = cube.parent.joinpath(cube.stem + '_extracted', self.extract['aperture_name'])
            aperture_dir.mkdir(parents=True, exist_ok=True)

            extractor = MIRI_IFU_Extractor(cube, aperture_dir.joinpath(cube.stem + '.fits'), ap_r=self.extract['aperture_radius'])
            extractor.extract_aperture()


    def extract_spaxels(self):
        raise Exception('UNIMPLEMENTED')
        # cubes = self.get_cubes()

        # for cube in cubes:
        #     spaxel_dir = cube.with_suffix('')
        #     if not spaxel_dir.exists():
        #         channel = int(re.findall(r'_ch\d+-shortmediumlong_s3d', cube.name)[0].split('-')[0][3:])
        #         specres = badass_miri.get_channel_resolution(channel, 'short')
        #         badass_miri.prepare_cube(cube, specres, z=self.redshift)


    def run_badass(self):
        if (not hasattr(self, 'badass')) or (not 'target' in self.badass):
            raise Exception('badass config expected')

        # TODO: fix
        if self.badass['target'] == '1D':
            self.run_badass_1D()
        elif self.badass['target'] == 'cube':
            self.run_badass_cube()
        else:
            raise Exception('unknown badass target: %s' % self.badass['target'])


    def run_badass_1D(self):
        if not 'line' in self.badass:
            raise Exception('Expected line')

        line = badass_miri.get_line(self.badass['line'])

        extracted_dirs = list(self.badass_dir.glob('%s_ch%d*_extracted' % (self.product_name, line.channel)))
        if len(extracted_dirs) != 1:
            raise Exception('Could not find extracted_dir')
        spec_fits = list(extracted_dirs[0].joinpath(self.extract['aperture_name']).glob('*.fits'))
        if len(spec_fits) != 1:
            raise Exception('Could not find spec fits file for channel %d' % line.channel)
        spec_fits = spec_fits[0]

        specres = badass_miri.get_channel_resolution(line.channel, line.subarray)
        diff = self.badass['wave_range']/2 if 'wave_range' in self.badass else 2000
        fit_reg = (line.wave-diff, line.wave+diff)

        badass_miri.run_badass_extracted(spec_fits, self.badass['options_file'], '%s_region'%self.badass['line'], specres, self.redshift, fit_reg=fit_reg)


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


class MIRI_IFU_Pipeline(Pipeline):
    def init(self):
        super().init()
        self.residual_fringe_correct = getattr(self, 'residual_fringe_correct', False)
        self.include_bkgd = getattr(self, 'include_bkgd', True)


    def run_stage2_single(self, rfile, output_dir, skip_cubes):
        print('Processing: {}'.format(str(rfile)))
        if output_dir.joinpath(rfile.name.replace('rate', 'cal')).exists():
            return None

        spec2 = Spec2Pipeline()
        spec2.output_dir = str(output_dir)
        spec2.output_file = str(rfile.with_suffix(''))

        if (skip_cubes):
            spec2.cube_build.skip = True
            spec2.extract_1d.skip = True

        if self.residual_fringe_correct:
            spec2.residual_fringe.skip = False

        spec2.run(rfile)
        spec2.suffix = 'spec2'
        return None


    def run_stage2_all(self, input_dir, output_dir, skip_cubes=False):
        rate_files = input_dir.glob('*_rate.fits')

        if not self.multiprocess:
            for rfile in rate_files:
                run_stage2_single(rfile, output_dir, skip_cubes)
            return

        args = [(rfile, output_dir, skip_cubes) for rfile in rate_files]
        pool = mp.Pool(processes=self.nprocesses, maxtasksperchild=1)
        pool.starmap(self.run_stage2_single, args, chunksize=1)
        pool.close()
        pool.join()


    def run_stage3(self, sci_input, bk_input, product_name, output_dir):
        sci_files = [str(f) for f in sci_input.glob('*_cal.fits')]
        asn = asn_from_list.asn_from_list(sci_files, rule=DMS_Level3_Base, product_name=product_name)

        if self.include_bkgd:
            bk_files = bk_input.glob('*_x1d.fits')
            for bkfile in bk_files:
                asn['products'][0]['members'].append({'expname': str(bkfile), 'exptype': 'background'})

        asnfile = output_dir.joinpath('%s_asn.json' % product_name)
        with open(asnfile, 'w') as afile:
            afile.write(asn.dump()[1])

        crds_config = Spec3Pipeline.get_config_from_reference(str(asnfile))
        spec3 = Spec3Pipeline.from_config_section(crds_config)

        spec3.output_dir = str(output_dir)
        spec3.save_results = True
        spec3.run(asnfile)


class MIRI_IFU_Extractor(Extractor):
    def __init__(self, file, outfile, ap_r='psf'):
        super().__init__(file, outfile, ap_r)

        self.scope_diam = 6.5 * u.m

        hdu = fits.open(self.file)
        hdr = hdu['SCI'].header

        # Linear wavelength solution in per-channel cubes
        self.wave = (np.arange(hdr['NAXIS3'])*hdr['CDELT3']+hdr['CRVAL3']) * u.Unit(hdr['CUNIT3'])
        pxar = hdr['PIXAR_SR'] * u.sr # pixel area in steradians
        self.cube_spec = (hdu['SCI'].data * u.Unit(hdr['BUNIT'])) * pxar
        self.var = ((hdu['ERR'].data * u.Unit(hdr['BUNIT'])) * pxar) ** 2
        self.pix_size = (hdr['CDELT1'] * u.Unit(hdr['CUNIT1'])).to(u.arcsec)
        hdu.close()

