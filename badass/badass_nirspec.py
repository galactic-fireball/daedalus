import astropy.io.fits as fits
import astropy.units as u
import pathlib
import numpy as np
import sys

ARGO = True # TODO: make generic toml option

if ARGO:
    # TODO: make generic toml option
    BADASS_DIR = pathlib.Path('/projects/ssatyapa/spectra/sdoan2/jwst/badass/')
    # BADASS_DIR = pathlib.Path(__file__).resolve().parent.parent.parent.joinpath('badass')
else:
    BADASS_DIR = pathlib.Path('/Users/sara/Dropbox/research/bgc/pubs/props/jwst_cycle2')

if BADASS_DIR.exists():
    sys.path.insert(0, str(BADASS_DIR))
    import badass_repo.badass as badass
    import badass_repo.badass_tools.badass_ifu as ifu

OPTIONS_DIR = pathlib.Path(__file__).resolve().parent.joinpath('options')

TRANSPOSE_FLUX = True # TODO: make generic toml option

# TODO: make generic toml option
# TARGET_FLUX_UNIT = u.erg / u.s / (u.cm**2) / u.um
TARGET_FLUX_UNIT = 1e-17 * u.erg / u.s / (u.cm**2) / u.AA
TARGET_WAVE_UNIT = u.AA


def prepare_cube(cube_fits, specres, z=0.0, plot=False):
    hdu = fits.open(cube_fits)

    dataext = 1
    CUNIT = 'CUNIT3'
    BUNIT = 'BUNIT'
    data = np.array((hdu[dataext].data).T, dtype='float32')
    header = hdu[dataext].header
    cunit = header[CUNIT]
    bunit = header[BUNIT]

    CRVAL = 'CRVAL3'
    CRPIX = 'CRPIX3'
    CDELT = 'CDELT3'
    nwave = data.shape[2]
    wave0 = header[CRVAL] - (header[CRPIX] - 1) * header[CDELT]
    wave = wave0 + np.arange(nwave)*header[CDELT]
    wave *= u.Unit(cunit)

    # TODO: check bunit
    pxar = header['PIXAR_SR'] * u.sr
    tf = (hdu['SCI'].data.T * u.Unit(bunit)) * pxar
    # tf = (hdu['SCI'].data.T * u.Unit(bunit))
    spec = tf.to(TARGET_FLUX_UNIT, equivalencies=u.spectral_density(wave)).value
    spec[np.isnan(spec)] = 0.0
    spec[spec < 0.0] = 0.0
    if TRANSPOSE_FLUX: spec = spec.T

    err = u.Quantity(hdu['ERR'].data.T, unit=u.Unit(hdu['ERR'].header['BUNIT']))
    err *= pxar # TODO: check bunit
    err = err.to(TARGET_FLUX_UNIT, equivalencies=u.spectral_density(wave)).value
    if TRANSPOSE_FLUX: err = err.T
    ivar = 1/(err**2)

    nx, ny, nz = header['NAXIS1'], header['NAXIS2'], header['NAXIS3']
    ra, dec = hdu[0].header['TARG_RA'], hdu[0].header['TARG_DEC']

    wave = wave.to(TARGET_WAVE_UNIT).value

    if isinstance(specres, dict):
        if len(specres['R']) == len(wave):
            specres = specres['R']
        else:
            specres = np.interp(wave, specres['wavelength'], specres['R'], left=specres['R'][0], right=specres['R'][-1])

    # TODO: in toml config
    aperture = [5, 45, 5, 45]

    wave,flux,ivar,mask,fwhm_res,binnum,\
    npixels,xpixbin,ypixbin,z,dataid,objname = ifu.prepare_ifu(str(cube_fits), z, format='user', ivar=ivar,
                                                                                  voronoi_binning=False,
                                                                                  fixed_binning=False,
                                                                                  fixed_bin_size=2,
                                                                                  use_and_mask=False,
                                                                                  nx=nx, ny=ny, nz=nz,
                                                                                  ra=ra, dec=dec,
                                                                                  wave=wave, flux=spec,
                                                                                  specres=specres,
                                                                                  aperture=aperture
                                                                                  )

    if plot:
        ifu.plot_ifu(cube_fits, wave, flux, ivar, mask, binnum, npixels, xpixbin, ypixbin, z, dataid)


def run_badass_spaxel(spaxel_dir, options_name):
    fits_file = list(spaxel_dir.glob('*.fits'))[0]
    print('Running spaxel fits: %s' % str(fits_file))

    run_dir = spaxel_dir.joinpath('spax_test') # TODO: fix
    options_file = OPTIONS_DIR.joinpath(options_name)
    badass.run_BADASS(fits_file, run_dir=run_dir, options_file=options_file, nprocesses=None, sdss_spec=False, ifu_spec=True)


def run_badass(spaxel_dir, options_name, output_name):
    options_file = OPTIONS_DIR.joinpath(options_name)
    if not options_file.exists():
        raise Exception('Could not find options file: %s' % str(options_file))

    # test_line_d = None
    # if test_line:
    #     test_line_d = {'bool':True,'line':'na_'+test_line}

    fits_files = []
    run_dirs = []
    spaxel_subdirs = list(spaxel_dir.glob('spaxel_*_*'))
    for spaxel in spaxel_subdirs:
        rd = spaxel.joinpath(output_name)
        if rd.exists():
            continue
        run_dirs.append(rd)
        rd.mkdir(exist_ok=True, parents=True)
        fits_files.append(spaxel.joinpath('%s.fits' % spaxel.name))
    if not fits_files:
        return

    print('Running badass on %d spaxels' % len(fits_files))
    badass.run_BADASS(fits_files, run_dir=run_dirs, options_file=options_file, nprocesses=48, sdss_spec=False, ifu_spec=True)


def run_cube_build(cube_file, output_name):
    _, _ = ifu.reconstruct_ifu(cube_file, output_name)
    ifu.plot_reconstructed_cube(cube_file.with_suffix('').joinpath('full_cube', output_name), animated=False)



# def prepare_cube(cube_fits, specres, z=0.0, aperture=None, plot=False):
#     hdu = fits.open(cube_fits)
#     sci_hdu = hdu['SCI']
#     hdr = sci_hdu.header
#     nx, ny, nz = hdu['SCI'].header['NAXIS1'], hdu['SCI'].header['NAXIS2'], hdu['SCI'].header['NAXIS3']
#     ra, dec = hdu[0].header['TARG_RA'], hdu[0].header['TARG_DEC']

#     wave_values = np.arange(hdr['NAXIS3'])*hdr['CDELT3']+hdr['CRVAL3']
#     wave_unit = hdr['CUNIT3']
#     wave_microns = u.Quantity(wave_values, unit=wave_unit)
#     wave = wave_microns.to(u.AA).value

#     flux_values = sci_hdu.data.T
#     flux_unit = u.Unit(hdr['BUNIT'])
#     flux = u.Quantity(flux_values, unit=flux_unit)

#     pxar = hdr['PIXAR_SR'] * u.sr
#     flux *= pxar

#     target_flux_unit = 1e-17 * u.erg / u.s / (u.cm**2) / u.AA
#     flux = flux.to(target_flux_unit, equivalencies=u.spectral_density(wave_microns)).value
#     flux = flux.T

#     err_values = u.Quantity(hdu['ERR'].data.T, unit=u.Unit(hdu['ERR'].header['BUNIT']))
#     err_values *= pxar
#     err_values = err_values.to(target_flux_unit, equivalencies=u.spectral_density(wave_microns)).value
#     err_values = err_values.T
#     ivar = 1/(err_values**2)

#     wave,flux,ivar,mask,fwhm_res,binnum,\
#     npixels,xpixbin,ypixbin,z,dataid,objname = ifu.prepare_ifu(str(cube_fits), z, format='user', ivar=ivar,
#                                                                                   voronoi_binning=False,
#                                                                                   fixed_binning=False,
#                                                                                   fixed_bin_size=2,
#                                                                                   use_and_mask=False,
#                                                                                   nx=nx, ny=ny, nz=nz,
#                                                                                   ra=ra, dec=dec,
#                                                                                   wave=wave, flux=flux,
#                                                                                   specres=specres
#                                                                                   )

#     if plot:
#         ifu.plot_ifu(cube_fits, wave, flux, ivar, mask, binnum, npixels, xpixbin, ypixbin, z, dataid)





