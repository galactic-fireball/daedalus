import astropy.units as u
from astropy.io import fits
import numpy as np
from scipy import interpolate
from specutils import Spectrum1D

import matplotlib.pyplot as plt
from matplotlib import colors
import plotly.express as px
import plotly.subplots

USE_SPECUTILS_READ = False

DEFAULT_WAVE_UNIT = u.um
DEFAULT_FLUX_UNIT = u.erg / u.s / (u.cm**2) / DEFAULT_WAVE_UNIT

JWST_X1D = 'JWST x1d'
JWST_S2D = 'JWST s2d'
JWST_S3D = 'JWST s3d'

file_format_default = 'JWST x1d'
file_format_dict = {
    'x1d': JWST_X1D,
    's2d': JWST_S2D,
    's3d': JWST_S3D,
}


def get_fmt_from_path(path):
    end = path.stem.split('_')[-1]
    if not end in file_format_dict:
        return file_format_default
    return file_format_dict[end]


def dered(wave, z):
    return wave / (1 + z)


class Spectrum:
    def __init__(self, file, wave, spec, err):
        self.file = file
        self.wave = wave
        self.spec = spec
        self.err = err


    @classmethod
    def from_x1d_file(cls, fits_file, z=0.0, wave_unit=DEFAULT_WAVE_UNIT):
        # TODO: implement custom
        fmt = get_fmt_from_path(fits_file)
        spectrum = Spectrum1D.read(fits_file, format=fmt)
        wave = dered(spectrum.spectral_axis, z).to(wave_unit)
        return cls(fits_file, wave, spectrum.flux, spectrum.uncertainty)


    @classmethod
    def from_s2d_file(cls, fits_file, z=0.0, wave_unit=DEFAULT_WAVE_UNIT):
        # TODO: implement custom
        fmt = get_fmt_from_path(fits_file)
        spectrum = Spectrum1D.read(fits_file, format=fmt)
        wave = dered(spectrum.spectral_axis, z).to(wave_unit)
        return cls(fits_file, wave, spectrum.flux, spectrum.uncertainty)


    @classmethod
    def from_s3d_file(cls, fits_file, z=0.0, wave_unit=DEFAULT_WAVE_UNIT):
        hdu = fits.open(fits_file)
        sci_hdu = hdu['SCI']
        hdr = sci_hdu.header

        wave_values = np.arange(hdr['NAXIS3'])*hdr['CDELT3']+hdr['CRVAL3']
        wave_unit = hdr['CUNIT3']
        wave = u.Quantity(wave_values, unit=wave_unit)
        wave = dered(wave, z)

        flux_values = sci_hdu.data.T
        flux_unit = u.Unit(hdr['BUNIT'])
        spec = u.Quantity(flux_values, unit=flux_unit)

        try:
            if spec.unit.to_string().split(' / ')[1] == 'sr':
                pxar = hdr['PIXAR_SR'] * u.sr
                spec *= pxar
        except:
            pass

        spec = spec.to(DEFAULT_FLUX_UNIT, equivalencies=u.spectral_density(wave))

        # TODO: propogate uncertainty correctly
        if not hdu['ERR'].header['ERRTYPE'] == 'ERR':
            print('WARNING: Unexpected error type')
            err = None
        else:
            err_values = u.Quantity(hdu['ERR'].data.T, unit=u.Unit(hdu['ERR'].header['BUNIT']))
            err_values *= pxar
            err_values = err_values.to(DEFAULT_FLUX_UNIT, equivalencies=u.spectral_density(wave))
            plt.figure()
            plt.matshow(err_values[:,:,400])
            plt.colorbar(location='bottom')
            plt.savefig('testerr.pdf')
            plt.close()
            # err = StdDevUncertainty(err_values, unit=hdu['ERR'].header['BUNIT'])

        #     err * pxar

        err = None

        return cls(fits_file, wave, spec, err)


    @classmethod
    def from_file(cls, fits_file, z=0.0, wave_unit=DEFAULT_WAVE_UNIT):
        fmt = fits_file.stem.split('_')[-1]
        read_func_name = 'from_%s_file' % fmt
        if not hasattr(cls, read_func_name):
            raise Exception('Unsupported option file type: %s' % fmt)

        return getattr(cls, read_func_name)(fits_file, z=z, wave_unit=wave_unit)


    @classmethod
    def combine_spectra(cls, spectra, dispersion=None, plot=False):
        if not dispersion:
            dispersion = spectra[0].wave
        wave_unit = dispersion.unit
        flux_unit = spectra[0].spec.unit

        if plot:
            plt.figure()

        fluxes = np.zeros((len(spectra), len(dispersion)))
        errors = np.zeros_like(fluxes)
        for i, spectrum in enumerate(spectra):

            iwave_unit = spectrum.wave.unit
            iwave = spectrum.wave.to(wave_unit, equivalencies=u.spectral())

            iflux_unit = spectrum.spec.unit
            iflux = spectrum.spec.to(flux_unit, equivalencies=u.spectral_density(iwave))

            inter_flux_func = interpolate.interp1d(iwave.value, iflux.value, bounds_error=False)
            fluxes[i,:] = inter_flux_func(dispersion.value)*flux_unit

            if plot:
                plt.plot(iwave, iflux, label='%d_before'%i)
                plt.plot(dispersion, fluxes[i], label='%d_after'%i)

            # TODO: error propagation okay?
            ierror = spectrum.err.quantity.to(flux_unit, equivalencies=u.spectral_density(iwave))
            inter_err_func = interpolate.interp1d(iwave.value, ierror.value, bounds_error=False)
            errors[i,:] = inter_err_func(dispersion.value)*flux_unit

        mean_flux = np.nanmean(fluxes, axis=0)*flux_unit
        err = np.sqrt(1.0 / np.sum(1.0 / (errors**2), axis=0))*flux_unit

        if plot:
            plt.plot(dispersion, mean_flux, label='mean flux')
            plt.xlabel('Wavelength ({})'.format(dispersion.unit))
            plt.ylabel('Flux ({})'.format(mean_flux.unit))
            plt.grid(True)
            plt.legend()
            plt.show()

        return cls(None, dispersion, mean_flux, err)


    def to_astropy_spec(self):
        return Spectrum1D(spectral_axis=self.wave, flux=self.spec, uncertainty=self.err)


    def plot_mean_flux(self):
        mean_flux = self.spec.sum(axis=(0,1)) / (len(self.spec)*len(self.spec[0]))

        plt.figure(figsize=(16,16))
        plt.plot(self.wave, mean_flux, label='Full spectrum', color='white', linewidth=0.5)
        plt.xlabel('Wavelength ({})'.format(self.wave.unit))
        plt.ylabel('Flux ({})'.format(mean_flux.unit))
        plt.title(self.file.stem)
        plt.grid(True)
        plt.show()


    def plot_spaxel_flux(self, idx):
        plt.figure(figsize=(16,16))
        plt.plot(self.wave, self.spec[idx[0],idx[1]], label='Full spectrum', color='black', linewidth=3)
        plt.xlabel('Wavelength ({})'.format(self.wave.unit))
        plt.ylabel('Flux ({})'.format(self.spec.unit))
        plt.title(self.file.stem)
        plt.grid(True)
        plt.show()


    def plot_all_spaxels(self, outdir):
        for i, j in [(i,j) for i in range(self.spec.shape[0]) for j in range(self.spec.shape[1])]:
            plt.figure()
            plt.plot(self.wave, self.spec[i,j], label='(%d, %d)' % (i, j), color='black', linewidth=0.5)
            plt.xlabel('Wavelength ({})'.format(self.wave.unit))
            plt.ylabel('Flux ({})'.format(self.spec.unit))
            plt.title('%s (%d, %d)' % (self.file.stem, i, j))
            plt.savefig(outdir.joinpath('%s_%d_%d.png' % (self.file.stem, i, j)))
            plt.close()


    def plot(self, lines=None, out=None):
        plt.figure()
        plt.plot(self.wave, self.spec, label='Full spectrum', color='black', linewidth=1.0)
        plt.xlabel('Wavelength ({})'.format(self.wave.unit))
        plt.ylabel('Flux ({})'.format(self.spec.unit))

        if lines:
            for name, line in lines.items():
                plt.axvline(x=line.to(self.wave.unit).value, color='red', linestyle='--')

        plt.title(self.file.stem if self.file else '')
        plt.grid(True)

        if out:
            plt.savefig(out)
            plt.close()
        else:
            plt.show()


    def plot_plotly(self, lines=None, out=None):
        fig = plotly.subplots.make_subplots(rows=1, cols=1)
        fig.add_trace(plotly.graph_objects.Scatter(x=self.wave, y=self.spec))

        fig.update_layout(
            yaxis_title='Flux ({})'.format(self.spec.unit),
            xaxis_title='Wavelength ({})'.format(self.wave.unit),
            title=self.file.stem if self.file else '',
            hovermode='x',
            template='plotly_white'
        )

        fig.update_xaxes(
            range=(np.nanmin(self.wave.value), np.nanmax(self.wave.value)),
            constrain='domain'
        )
        fig.update_yaxes(
            range=(np.nanmin(self.spec.value)-.3, np.nanmax(self.spec.value)+.3),
            constrain='domain'
        )

        if lines:
            for name, line in lines.items():
                fig.add_vline(x=line.to(self.wave.unit).value, line_dash='dash', line_color='red', annotation_text=name)

        if out:
            fig.write_html(out, include_mathjax="cdn")
        else:
            fig.show()


    def plot_flux_map(self):
        plt.figure(figsize=(16,16))
        m = plt.imshow(np.nansum(self.spec, axis=2)/self.spec.shape[2], cmap='jet', norm=colors.LogNorm())
        plt.colorbar(m)
        plt.show()


    def plot_flux_map_plotly(self):
        # fig = px.imshow(np.nansum(self.spec, axis=2)/self.spec.shape[2])
        # fig.update_layout(dragmode='drawrect', newshape=dict(line_color='cyan'))
        # fig.show(config={'modeBarButtonsToAdd': ['drawcircle', 'drawrect', 'eraseshape', 'hoverClosestCartesian']})

        fig = px.imshow(self.spec, animation_frame=2, labels={'color':str(self.spec.unit), 'animation_frame':'slice'})
        fig.layout.updatemenus[0].buttons[0].args[1]['frame']['duration'] = 10
        fig.layout.updatemenus[0].buttons[0].args[1]['transition']['duration'] = 1
        fig.show()



def get_dered_spectrum(fits_file, z=0.0, wave_unit=DEFAULT_WAVE_UNIT):
    if z == 0.0:
        print('WARNING: de-redshifting spectrum with z = 0.0')

    if USE_SPECUTILS_READ:
        fmt = get_fmt_from_path(fits_file)
        spec = Spectrum1D.read(fits_file, format=fmt)
        wave = dered(spec.spectral_axis, z).to(wave_unit)
        spec_keys = ['uncertainty', 'mask', 'velocity_convention', 'rest_value', 'meta']
        return Spectrum1D(spectral_axis=wave, flux=spec.flux, **{key: getattr(obs_spec, key) for key in spec_keys})

    return Spectrum.from_file(fits_file, z=z, wave_unit=wave_unit).to_astropy_spec()



