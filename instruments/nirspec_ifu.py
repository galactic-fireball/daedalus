import astropy.io.fits as fits
import astropy.units as u
from astroquery.mast import Observations
from datetime import datetime
import json
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np

import mast.mast as mast
import pipeline.crds_utils as crds_utils
import badass.badass_nirspec as badass_nirspec
from instruments.instrument_common import Instrument, Pipeline, Extractor

import jwst
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline
from jwst.pipeline.calwebb_spec3 import Spec3Pipeline


STAGE2_ASN = True


class NIRSpec_IFU(Instrument):

    def download_all(self):
        mast.login(self.mast_token)
        df = mast.get_program_data(str(self.program_id), 'NIRSPEC/IFU')

        # uncal
        uncal_df = Observations.filter_products(df, calib_level=[1], productSubGroupDescription='UNCAL').to_pandas()
        obs_fmt = 'jw%05d%03d' % (self.program_id, self.obs)
        obs_df = uncal_df[uncal_df.obs_id.str.contains(obs_fmt)]

        for file_name in obs_df.productFilename.values:
            dest = self.pipeline_dir.joinpath(file_name)
            if dest.exists():
                continue
            print('downloading %s' % file_name)
            mast.download_file(file_name, dest=dest)

        # asn2
        asn2_df = Observations.filter_products(df, calib_level=[2], productSubGroupDescription='ASN').to_pandas()
        obs_fmt = 'jw%05d%03d' % (self.program_id, self.obs)
        obs_df = asn2_df[(asn2_df.obs_id.str.contains(obs_fmt)) & (asn2_df.productFilename.str.contains('spec'))]

        for file_name in obs_df.productFilename.values:
            dest = self.pipeline_dir.joinpath(file_name)
            if dest.exists():
                continue
            print('downloading %s' % file_name)
            mast.download_file(file_name, dest=dest)

        # asn3
        asn3_df = Observations.filter_products(df, calib_level=[3], productSubGroupDescription='ASN').to_pandas()
        obs_fmt = 'jw%05d-o%03d' % (self.program_id, self.obs)
        obs_df = asn3_df[asn3_df.obs_id.str.contains(obs_fmt)]

        for file_name in obs_df.productFilename.values:
            dest = self.pipeline_dir.joinpath(file_name)
            if dest.exists():
                continue
            print('downloading %s' % file_name)
            mast.download_file(file_name, dest=dest)


    def run_pipeline(self):
        # crds_utils.cache_crds('nirspec')

        pipeline_config = getattr(self, 'pipeline', {})
        pipeline = NIRSpec_IFU_Pipeline(pipeline_config)

        # pipeline.run_stage1_all(self.pipeline_dir, self.pipeline_dir)
        # pipeline.run_stage2_all(self.pipeline_dir, self.pipeline_dir)
        pipeline.run_stage3_all(self.pipeline_dir, self.pipeline_dir)
        print('Pipeline for {pname} complete!'.format(pname=self.product_name))


    def get_cubes(self):
        # TODO: use product_name once the asn file is updated correctly
        cubes = list(self.pipeline_dir.glob('*_nirspec_*s3d.fits'))
        print('Found {n} cubes'.format(n=len(cubes)))
        for cube in cubes:
            cube.rename(self.badass_dir.joinpath(cube.name))

        # grab the list this way since we might not have actually moved anything
        return list(self.badass_dir.glob('*_s3d.fits'))
        # return list(self.badass_dir.glob('*g140h-f100lp_s3d.fits'))


    def extract_1D(self):
        cubes = self.get_cubes()
        for cube in cubes:
            print('Extracting %s' % str(cube))
            aperture_dir = cube.parent.joinpath(cube.stem + '_extracted', self.extract['aperture_name'])
            aperture_dir.mkdir(parents=True, exist_ok=True)

            extractor = NIRSpec_IFU_Extractor(cube, aperture_dir.joinpath(cube.stem + '.fits'), ap_r=self.extract['aperture_radius'], redshift=self.redshift)
            extractor.extract_aperture()


    def get_spec_res(self, fits_file):
        DISP_FILE_FMT = 'jwst_nirspec_%s_disp.fits'
        grating = fits_file.name.split('_')[-2].split('-')[0]
        disp_file = self.inst_data_dir.joinpath(DISP_FILE_FMT % grating)
        if not disp_file.exists():
            raise Exception('Could not find dispersion file: %s' % str(disp_file))

        hdu = fits.open(disp_file)
        return dict(wavelength=hdu[1].data['WAVELENGTH'], R=hdu[1].data['R'])


    def run_badass(self):
        if (not hasattr(self, 'badass')) or (not 'target' in self.badass):
            raise Exception('badass config expected')

        # TODO: fix
        if self.badass['target'] == '1D':
            raise Exception('unimplemented')
        elif self.badass['target'] == 'cube':
            self.run_badass_cube()
        else:
            raise Exception('unknown badass target: %s' % self.badass['target'])


    def run_badass_cube(self):
        cubes = self.get_cubes()
        output_name = 'cube_run1'
        for cube in cubes:
            print('Preparing cube: %s' % str(cube))

            specres = self.get_spec_res(cube)
            badass_nirspec.prepare_cube(cube, specres, z=self.redshift, plot=True)
            spaxel_dir = cube.with_suffix('')
            badass_nirspec.run_badass(spaxel_dir, self.badass['options_file'], output_name)
            badass_nirspec.run_cube_build(cube, output_name)
            # badass_nirspec.run_badass_spaxel(spaxel_dir.joinpath('spaxel_21_23'), self.badass['options_file'])
            # 'spaxel_21_23/spaxel_21_23.fits'
            # badass_nirspec.run_badass(spaxel_dir, self.badass['options_file'], output_name, test_line=None, fit_reg=None)



class NIRSpec_IFU_Pipeline(Pipeline):
    def update_asn(self, asn_temp):
        data = json.load(open(asn_temp))

        name = data['products'][0]['name']
        data['code_version'] = jwst.__version__
        data['version_id'] = datetime.utcnow().strftime('%Y%m%dt%H%M%S')
        # data['asn_pool'] = ''

        asn_file = asn_temp.parent.joinpath('%s_updated.json' % asn_temp.stem)
        json.dump(data, open(asn_file, 'w'))

        return asn_file, name

    def run_stage2_single(self, infile, output_dir):
        print('Processing: {}'.format(str(infile)))

        # stage2_params = getattr(self, 'stage2', None)
        # step_params = getattr(stage2_params, 'steps', {})
        # breakpoint()

        if STAGE2_ASN:
            asn_file, name = self.update_asn(infile)

            out_file = output_dir.joinpath('%s_cal.fits'%name)
            if out_file.exists():
                return None

            # crds_config = Spec2Pipeline.get_config_from_reference(str(asn_file))
            # spec2 = Spec2Pipeline.from_config_section(crds_config)

            # spec2.output_dir = str(output_dir)
            # spec2.save_results = True
            # spec2.run(asn_file)

            spec2 = Spec2Pipeline()
            spec2.save_results = True
            spec2.output_dir = str(output_dir)

            # for step, params in step_params.items():
            #     for key, value in params.items():
            #         setattr(spec2, step)

            spec2.cube_build.weighting = 'emsm'

            result = spec2(str(asn_file))
            print('spec2 result: %s' % result)

            plt.figure()
            hdu = fits.open(out_file)
            plt.imshow(hdu['SCI'].data)
            plt.colorbar()
            plt.savefig(output_dir.joinpath('test_plots', out_file.stem + '.png'))

            return None

        # Run stage 2 directly on rate files
        if output_dir.joinpath(infile.name.replace('rate', 'cal')).exists():
            return None

        spec2 = Spec2Pipeline()
        spec2.output_dir = str(output_dir)
        spec2.output_file = str(infile.with_suffix(''))

        spec2.run(infile)
        spec2.suffix = 'spec2'
        return None


    def run_stage2_all(self, indir, output_dir):
        if STAGE2_ASN:
            in_files = indir.glob('*_spec2_*_asn.json')
        else:
            in_files = indir.glob('*_rate.fits')

        if not self.multiprocess:
            for in_file in in_files:
                self.run_stage2_single(in_file, output_dir)
            return

        args = [(in_file, output_dir) for in_file in in_files]
        pool = mp.Pool(processes=self.nprocesses, maxtasksperchild=1)
        pool.starmap(self.run_stage2_single, args, chunksize=1)
        pool.close()
        pool.join()


    def run_stage3_single(self, asn_temp, out_dir):
        print('Processing: {}'.format(str(asn_temp)))
        asn_file, name = self.update_asn(asn_temp)

        # if len(list(out_dir.glob('%s*_s3d.fits'%name))) == 1:
        #     return None

        crds_config = Spec3Pipeline.get_config_from_reference(str(asn_file))
        spec3 = Spec3Pipeline.from_config_section(crds_config)

        spec3.outlier_detection.skip = True
        # spec3.extract_1d.center_xy = (29, 26)

        spec3.output_dir = str(out_dir)
        spec3.save_results = True
        spec3.run(asn_file)
        return None


    def run_stage3_all(self, indir, output_dir):
        asn_files = indir.glob('*_spec3_*_asn.json')

        if not self.multiprocess:
            for asn_file in asn_files:
                self.run_stage3_single(asn_file, output_dir)
            return

        args = [(asn_file, output_dir) for asn_file in asn_files]
        pool = mp.Pool(processes=self.nprocesses, maxtasksperchild=1)
        pool.starmap(self.run_stage3_single, args, chunksize=1)
        pool.close()
        pool.join()


class NIRSpec_IFU_Extractor(Extractor):
    def __init__(self, file, outfile, ap_r='psf', redshift=0.0):
        super().__init__(file, outfile, ap_r=ap_r, redshift=redshift)

        self.scope_diam = 6.5 * u.m

        hdu = fits.open(self.file)
        hdr = hdu['SCI'].header

        # Linear wavelength solution in per-channel cubes
        self.wave = (np.arange(hdr['NAXIS3'])*hdr['CDELT3']+hdr['CRVAL3']) * u.Unit(hdr['CUNIT3'])
        # pxar = hdr['PIXAR_SR'] * u.sr # pixel area in steradians
        # self.cube_spec = (hdu['SCI'].data * u.Unit(hdr['BUNIT'])) * pxar
        # self.cube_spec = (hdu['SCI'].data * u.Unit(hdr['BUNIT'])) * u.sr
        # self.var = ((hdu['ERR'].data * u.Unit(hdr['BUNIT'])) * pxar) ** 2

        # TARGET_FLUX_UNIT = u.erg / u.s / (u.cm**2) / u.um
        self.cube_spec = (hdu['SCI'].data * u.Unit(hdr['BUNIT']))
        # # self.cube_spec = self.cube_spec.to(TARGET_FLUX_UNIT, equivalencies=u.spectral_density(self.wave))
        self.cube_spec[np.isnan(self.cube_spec)] = 0.0
        self.cube_spec[self.cube_spec < 0.0] = 0.0

        self.var = ((hdu['ERR'].data * u.Unit(hdr['BUNIT']))) ** 2
        self.pix_size = (hdr['CDELT1'] * u.Unit(hdr['CUNIT1'])).to(u.arcsec)
        hdu.close()




