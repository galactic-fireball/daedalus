import matplotlib.pyplot as plt
import multiprocessing as mp
import os
import pathlib
from prodict import Prodict

from astropy.io import fits
from astropy import units as u

from photutils.aperture import CircularAperture, CircularAnnulus, ApertureStats
from photutils.detection import DAOStarFinder

import numpy as np

import ssl
# Needed to work around ssl certificate verification during crds downloads
ssl._create_default_https_context = ssl._create_unverified_context

USE_CRDS_OPS = True
# Needs to be set before crds/jwst imports
if USE_CRDS_OPS:
    os.environ['CRDS_PATH'] = str(pathlib.Path(__file__).resolve().parent.parent.joinpath('pipeline', 'crds_cache', 'ops'))
    os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu'
else:
    os.environ['CRDS_PATH'] = str(pathlib.Path(__file__).resolve().parent.parent.joinpath('pipeline', 'crds_cache', 'pub'))
    os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds-pub.stsci.edu'

import crds
import jwst
from jwst.pipeline.calwebb_detector1 import Detector1Pipeline

JWST_VERSION_STR = 'pipeline_%s' % jwst.__version__

PROGRAMS_DIR = pathlib.Path(__file__).parent.parent.joinpath('programs')

class Instrument(Prodict):
    def from_config(config, instruments):

        for attr in ['target_name', 'program_id', 'instrument', 'actions']:
            if attr not in config:
                raise Exception('No \'%s\' in configuration' % attr)

        instrument = config['instrument']
        if not instrument in instruments:
            raise Exception('Unknown instrument: %s' % instrument)

        return instruments[instrument](config)


    def init(self):
        self.product_name = '%s_%s' % (self.target_name, JWST_VERSION_STR)
        self.data_dir = PROGRAMS_DIR.joinpath(str(self.program_id), 'data_sets', self.target_name, self.instrument, JWST_VERSION_STR)
        self.pipeline_dir = self.data_dir.joinpath('pipeline')
        self.pipeline_dir.mkdir(parents=True, exist_ok=True)
        self.badass_dir = self.data_dir.joinpath('badass')
        self.badass_dir.mkdir(parents=True, exist_ok=True)


    def run(self):
        print(self.__dict__)
        if not all(hasattr(self, action) for action in self.actions):
            raise Exception('Unknown action in: %s' % self.actions)

        for action in self.actions:
            getattr(self, action)()


class Pipeline(Prodict):
    def init(self):
        for attr in []:
            if not hasattr(self, attr):
                raise Exception('\'%s\' not provided' % attr)

        self.multiprocess = getattr(self, 'multiprocess', False)
        self.nprocesses = getattr(self, 'nprocesses', 4)


    def run_stage1_single(self, ufile, output_dir):
        print('Processing: {}'.format(str(ufile)))
        if output_dir.joinpath(ufile.name.replace('uncal', 'rate')).exists():
            return None

        detector1 = Detector1Pipeline()
        detector1.output_dir = str(output_dir)
        detector1.output_file = str(output_dir.joinpath(ufile.stem))
        detector1.run(ufile)
        return None


    def run_stage1_all(self, uncal_dir, output_dir):
        uncal_files = uncal_dir.glob('*_uncal.fits')

        if not self.multiprocess:
            for ufile in uncal_files:
                run_stage1_single(ufile, output_dir)
            return

        args = [(ufile, output_dir) for ufile in uncal_files]
        pool = mp.Pool(processes=self.nprocesses, maxtasksperchild=1)
        pool.starmap(self.run_stage1_single, args, chunksize=1)
        pool.close()
        pool.join()



TARGET_FLUX_UNIT = 1e-17 * u.erg / u.s / (u.cm**2) / u.AA

class Extractor:
    def __init__(self, file, outfile, ap_r='psf', plot=True):
        self.file = file
        self.outfile = outfile
        self.aperture_radius = ap_r

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
        print(source['xcentroid','ycentroid','flux'])
        position = np.transpose((source['xcentroid'], source['ycentroid']))
        return position


    def extract_aperture(self):
        source = self.find_source()

        spec = np.zeros(len(self.wave))*self.cube_spec.unit
        var = np.zeros(len(self.wave))*self.var.unit

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
            bkgd = ApertureStats(self.cube_spec[i,:,:], annulus_aperture).mean * aperture.area
            bkgd = bkgd if not np.isnan(bkgd) else 0.0
            spec[i] = ap_tot - bkgd

            ap_var_tot = ApertureStats(self.var[i,:,:],aperture).sum
            bkgd_var_tot = ApertureStats(self.var[i,:,:], annulus_aperture).sum
            var[i] = ap_var_tot+(bkgd_var_tot*aperture.area/annulus_aperture.area)

        spec = spec.T.to(TARGET_FLUX_UNIT, equivalencies=u.spectral_density(self.wave)).T.value
        err = np.sqrt(var).T.to(TARGET_FLUX_UNIT, equivalencies=u.spectral_density(self.wave)).T.value

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
            plt.plot(self.wave.value, spec, linewidth=0.5)
            plt.xlabel('Wavelength (%s)' % self.wave.unit.to_string())
            plt.ylabel('Flux (%s)' % TARGET_FLUX_UNIT.to_string())
            plt.savefig(self.outfile.parent.joinpath(self.outfile.stem+'_spec.png'), dpi=300)

