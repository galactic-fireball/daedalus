import matplotlib.pyplot as plt
import multiprocessing as mp
import os
import pathlib
from prodict import Prodict

from astropy.io import fits
from astropy.stats import sigma_clip
from astropy import units as u

from photutils.aperture import CircularAperture, CircularAnnulus, ApertureStats
from photutils.detection import DAOStarFinder

import numpy as np

import ssl
# Needed to work around ssl certificate verification during crds downloads
ssl._create_default_https_context = ssl._create_unverified_context

from utils import plotly_plot

from jwst.pipeline.calwebb_detector1 import Detector1Pipeline

INSTRUMENTS_DIR = pathlib.Path(__file__).parent


def configure_crds(cache=None, use_ops=True):
    if use_ops:
        os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu'
        if cache:
            os.environ['CRDS_PATH'] = str(pathlib.Path(cache).joinpath('ops'))
    else:
        os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds-pub.stsci.edu'
        if cache:
            os.environ['CRDS_PATH'] = str(pathlib.Path(cache).joinpath('pub'))


class Instrument(Prodict):

    def run_pipeline(self, context, args):
        self.run_stage1_all(context, args)
        self.run_stage2_all(context, args)
        self.run_stage3_all(context, args)
        print('Pipeline for {pname} complete!'.format(pname=context.config.product_name))


    def run_stage1_single(self, ufile, output_dir, context, args):
        print('Processing: {}'.format(str(ufile)))
        out_file = output_dir.joinpath(ufile.name.replace('uncal', 'rate'))
        if out_file.exists():
            # TODO
            # if CLEAN_RATES:
            #     self.clean_rate(out_file)
            return None

        detector1 = Detector1Pipeline()
        detector1.output_dir = str(output_dir)
        detector1.output_file = str(output_dir.joinpath(ufile.stem))

        step_opts = args.get('stage1', {}).get('steps', {})
        for step, options in step_opts.items():
            for opt, val in options.items():
                setattr(getattr(detector1, step), opt, val)

        detector1.run(ufile)

        # TODO
        # if CLEAN_RATES:
        #     self.clean_rate(out_file)

        return None


    def run_stage1_all(self, context, args):
        stage1_opts = args.get('stage1', {})
        input_dir = stage1_opts.get('input_dir', context.pipeline_dir)
        output_dir = stage1_opts.get('output_dir', context.pipeline_dir)
        multiprocess = args.get('multiprocess', False)
        nprocesses = args.get('nprocesses', 1)

        uncal_files = input_dir.glob('*_uncal.fits')

        if not multiprocess:
            for ufile in uncal_files:
                self.run_stage1_single(ufile, output_dir, context, args)
            return

        proc_args = [(ufile, output_dir, context, args) for ufile in uncal_files]
        pool = mp.Pool(processes=self.nprocesses, maxtasksperchild=1)
        pool.starmap(self.run_stage1_single, proc_args, chunksize=1)
        pool.close()
        pool.join()



class Pipeline(Prodict):
    def init(self):
        for attr in []:
            if not hasattr(self, attr):
                raise Exception('\'%s\' not provided' % attr)

        self.multiprocess = getattr(self, 'multiprocess', False)
        self.nprocesses = getattr(self, 'nprocesses', 4)


    



# TARGET_FLUX_UNIT = 1e-17 * u.erg / u.s / (u.cm**2) / u.AA
TARGET_FLUX_UNIT = u.erg / u.s / (u.cm**2) / u.AA
# TARGET_FLUX_UNIT = u.Unit('Jy')

class Extractor:
    def __init__(self, file, outfile, ap_r='psf', plot=True, redshift=0.0):
        self.file = file
        self.outfile = outfile
        self.aperture_radius = ap_r
        self.redshift = redshift

        if isinstance(self.aperture_radius, u.Quantity):
            self.aperture_radius.to(u.arcsec)

        if isinstance(self.aperture_radius, (int, float)):
            self.aperture_radius *= u.arcsec

        self.ap_figure = plt.figure() if plot else None
        self.spec_figure = plt.figure() if plot else None


    def find_source(self):
        medcube = np.nanmedian(self.cube_spec.value,axis=0)
        medcube[np.isnan(medcube)] = 0.0

        if self.ap_figure:
            plt.figure(self.ap_figure)
            plt.imshow(medcube)

        daofind = DAOStarFinder(fwhm=3.0, threshold=1e-17)
        sources = daofind(medcube)
        if not sources:
            raise Exception('Failed to find centroid source for: %s' % str(self.file))
        source = sources[np.argmax(sources['flux'])]

        # TODO: configuration to set a different source?
        # source = sources.to_pandas().sort_values(by=['flux']).iloc[-2]

        print(source['xcentroid','ycentroid','flux'])
        position = np.transpose((source['xcentroid'], source['ycentroid']))
        return position


    def extract_aperture(self):
        source = self.find_source()

        spec = np.zeros(len(self.wave))*self.cube_spec.unit
        bkgd = np.zeros(len(self.wave))*self.cube_spec.unit
        var = np.zeros(len(self.wave))*self.var.unit

        # for i in range(self.cube_spec.shape[1]):
        #     for j in range(self.cube_spec.shape[2]):
        #         self.cube_spec[:,i,j] = sigma_clip(self.cube_spec[:,i,j].value)*self.cube_spec.unit

        # self.cube_spec = sigma_clip(self.cube_spec.value, axis=0, sigma=5)*self.cube_spec.unit

        if self.aperture_radius != 'psf':
            radius = (self.aperture_radius / self.pix_size).value
            aperture = CircularAperture(source, r=radius)
            annulus_aperture = CircularAnnulus(source, r_in=radius*2, r_out=radius*2.5)
            if self.ap_figure:
                plt.figure(self.ap_figure)
                aperture.plot()
                annulus_aperture.plot()

        mid_idx = int(len(self.wave)/2)
        for i in range(0, len(self.wave)):
            if self.aperture_radius == 'psf':
                res = (1.22 * (self.wave[i].to(self.scope_diam.unit) / self.scope_diam)) * u.rad
                radius = res.to(u.arcsec).value
                aperture = CircularAperture(source, r=radius)
                annulus_aperture = CircularAnnulus(source, r_in=radius*2, r_out=radius*2.5)
                if (i == mid_idx) and (self.ap_figure):
                    plt.figure(self.ap_figure)
                    aperture.plot()
                    annulus_aperture.plot()

            ap_tot = ApertureStats(self.cube_spec[i,:,:], aperture).sum
            bkgd[i] = ApertureStats(self.cube_spec[i,:,:], annulus_aperture).mean * aperture.area
            # bkgd = ApertureStats(self.cube_spec[i,:,:], annulus_aperture).mean
            bkgd[i] = bkgd[i] if not np.isnan(bkgd[i]) else 0.0
            bkgd[i] = bkgd[i] if bkgd[i].value < 1e-6 else bkgd[i-1]
            spec[i] = ap_tot# - bkgd[i]

            ap_var_tot = ApertureStats(self.var[i,:,:],aperture).sum
            bkgd_var_tot = ApertureStats(self.var[i,:,:], annulus_aperture).sum
            var[i] = ap_var_tot+(bkgd_var_tot*aperture.area/annulus_aperture.area)

        spec = spec.T.to(TARGET_FLUX_UNIT, equivalencies=u.spectral_density(self.wave)).T.value
        bkgd = bkgd.T.to(TARGET_FLUX_UNIT, equivalencies=u.spectral_density(self.wave)).T.value
        err = np.sqrt(var).T.to(TARGET_FLUX_UNIT, equivalencies=u.spectral_density(self.wave)).T.value

        spec = sigma_clip(spec, sigma=4)

        fits_data = fits.BinTableHDU.from_columns(fits.ColDefs([
            fits.Column(name='WAVE', array=self.wave, format='D'),
            fits.Column(name='SPEC', array=spec, format='D'),
            fits.Column(name='ERR', array=err, format='D'),
        ]))

        fits_data.header.append(('WAVEUNIT', self.wave.unit.to_string(), 'wavelength unit'))
        fits_data.header.append(('SPECUNIT', TARGET_FLUX_UNIT.to_string(), 'spectral flux unit'))

        # TODO: create an actual primary HDU with RA, DEC, etc. from original fits
        # TODO: include aperture radius
        hdu = fits.HDUList([fits.PrimaryHDU(), fits_data])
        hdu.writeto(self.outfile, overwrite=True)

        if self.ap_figure:
            plt.figure(self.ap_figure)
            plt.savefig(self.outfile.parent.joinpath(self.outfile.stem+'_aperture.png'), dpi=300)

        if self.spec_figure:
            plt.figure(self.spec_figure)
            # plt.yscale('log')
            nz_mask = np.where(spec > 1e-18)
            plt.plot(self.wave.value[nz_mask], spec[nz_mask], linewidth=0.5, label='Data')
            # plt.plot(self.wave.value, bkgd, linewidth=0.5, label='Background')
            # plt.legend()
            plt.xlabel('Wavelength (%s)' % self.wave.unit.to_string())
            plt.ylabel('Flux (%s)' % TARGET_FLUX_UNIT.to_string())
            plt.savefig(self.outfile.parent.joinpath(self.outfile.stem+'_spec.png'), dpi=300)

            plotly_plot.plot_fits(self.outfile, z=self.redshift)

