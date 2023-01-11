import argparse
from astropy.io import fits
import astropy.units as u
import importlib
import os
import matplotlib.pyplot as plt
import numpy as np
import pathlib
from shutil import copyfile
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent)) # for miri_consts

BADASS_DIR = pathlib.Path('/projects/ssatyapa/spectra/sdoan2/jwst/badass/')
# BADASS_DIR = pathlib.Path(__file__).resolve().parent.parent.parent.joinpath('badass')

if BADASS_DIR.exists():
    sys.path.insert(0, str(BADASS_DIR))
    import badass_repo.badass as badass
    import badass_repo.badass_tools.badass_ifu as ifu

import miri_consts as consts

OPTIONS_DIR = pathlib.Path(__file__).resolve().parent


def get_line(line_name):
    for line in consts.lines:
        if line.name == line_name:
            return line
    return None


def get_channel_resolution(channel, subarray):
    return consts.resolutions[channel][subarray]


def prepare_cube(cube_fits, specres, z=0.0, aperture=None, plot=False):
    hdu = fits.open(cube_fits)
    sci_hdu = hdu['SCI']
    hdr = sci_hdu.header
    nx, ny, nz = hdu['SCI'].header['NAXIS1'], hdu['SCI'].header['NAXIS2'], hdu['SCI'].header['NAXIS3']
    ra, dec = hdu[0].header['TARG_RA'], hdu[0].header['TARG_DEC']

    wave_values = np.arange(hdr['NAXIS3'])*hdr['CDELT3']+hdr['CRVAL3']
    wave_unit = hdr['CUNIT3']
    wave_microns = u.Quantity(wave_values, unit=wave_unit)
    wave = wave_microns.to(u.AA).value

    flux_values = sci_hdu.data.T
    flux_unit = u.Unit(hdr['BUNIT'])
    flux = u.Quantity(flux_values, unit=flux_unit)

    pxar = hdr['PIXAR_SR'] * u.sr
    flux *= pxar

    target_flux_unit = 1e-17 * u.erg / u.s / (u.cm**2) / u.AA
    flux = flux.to(target_flux_unit, equivalencies=u.spectral_density(wave_microns)).value
    flux = flux.T

    err_values = u.Quantity(hdu['ERR'].data.T, unit=u.Unit(hdu['ERR'].header['BUNIT']))
    err_values *= pxar
    err_values = err_values.to(target_flux_unit, equivalencies=u.spectral_density(wave_microns)).value
    err_values = err_values.T
    ivar = 1/(err_values**2)

    wave,flux,ivar,mask,fwhm_res,binnum,\
    npixels,xpixbin,ypixbin,z,dataid,objname = ifu.prepare_ifu(str(cube_fits), z, format='user', ivar=ivar,
                                                                                  voronoi_binning=False,
                                                                                  fixed_binning=False,
                                                                                  fixed_bin_size=2,
                                                                                  use_and_mask=False,
                                                                                  nx=nx, ny=ny, nz=nz,
                                                                                  ra=ra, dec=dec,
                                                                                  wave=wave, flux=flux,
                                                                                  specres=specres
                                                                                  )

    if plot:
        ifu.plot_ifu(cube_fits, wave, flux, ivar, mask, binnum, npixels, xpixbin, ypixbin, z, dataid)


def run_badass(spaxel_dir, options_name, output_name, test_line=None, fit_reg=None):
    options_file = OPTIONS_DIR.joinpath(options_name)
    if not options_file.exists():
        raise Exception('Could not find options file: %s' % str(options_file))

    test_line_d = None
    if test_line:
        test_line_d = {'bool':True,'line':'na_'+test_line}

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

    badass.run_BADASS(fits_files, run_dir=run_dirs, options_file=options_file, nprocesses=4, test_line=test_line_d, fit_reg=fit_reg, sdss_spec=False, ifu_spec=True)


def run_badass_extracted(fits_file, options_name, output_name, specres, z, test_line=None, fit_reg=None):
    options_file = OPTIONS_DIR.joinpath(options_name)
    if not options_file.exists():
        raise Exception('Could not find options file: %s' % str(options_file))

    run_dir = fits_file.parent.joinpath(fits_file.stem, output_name)
    run_dir.mkdir(exist_ok=True, parents=True)

    test_line_d = None
    if test_line:
        test_line_d = {'bool':True,'line':'na_'+test_line}

    hdu = fits.open(fits_file)
    # TODO: use hdu[1].header['WAVEUNIT']
    wave = hdu[1].data['wave']
    wave = (wave*u.um).to(u.AA).value
    spec = hdu[1].data['spec']
    spec = spec*1.e-17
    err = hdu[1].data['err']
    err = err*1.e-17
    fwhm_res = wave / specres

    badass.run_BADASS(fits_file, run_dir=run_dir, options_file=options_file, test_line=test_line_d, fit_reg=fit_reg, sdss_spec=False, ifu_spec=False, wave=wave, spec=spec, err=err, fwhm_res=fwhm_res, z=z)


def reconstruct_cube(cube_fits, out_label, line_test=False):
    ifu.reconstruct_ifu(cube_fits, out_label, line_test=line_test)


def plot_cube(cube_fits, out_label):
    cube_out_dir = cube_fits.with_suffix('').joinpath('full_cube', out_label)
    ifu.plot_reconstructed_cube(cube_out_dir, animated=False)


def plot_ratio(cube1, name1, cube2, name2, out_plot):
    partable = cube1.joinpath('log', 'cube_par_table.fits')
    if not partable.exists():
        raise Exception('Unable to find par table: %s' % str(partable))
    parhdu = fits.open(partable)
    if not name1 in parhdu:
        print('%s not found' % name1)
        return

    ox, oy = parhdu[0].header['ORIGINX'], parhdu[0].header['ORIGINY']
    data1 = parhdu[name1].data
    mask = data1 >= 0
    data1[mask] = np.nan
    data1 = 10**data1

    partable = cube2.joinpath('log', 'cube_par_table.fits')
    if not partable.exists():
        raise Exception('Unable to find par table: %s' % str(partable))
    parhdu = fits.open(partable)
    if not name2 in parhdu:
        print('%s not found' % name2)
        return

    data2 = parhdu[name2].data
    mask = data2 >= 0
    data2[mask] = np.nan
    data2 = 10**data2

    if data1.shape != data2.shape:
        print('Skipping %s vs %s' % (name1, name2))
        return

    data = data1 / data2

    fig, ax = plt.subplots()
    map_ = ax.imshow(data, origin='lower', cmap='jet',
              vmin=np.nanpercentile(data, 1),
              vmax=np.nanpercentile(data, 99),
              extent=[ox-.5, ox+data.shape[1]-.5, oy-.5, oy+data.shape[0]-.5])
    plt.colorbar(map_, ax=ax, label=out_plot.stem)
    ax.set_title(out_plot.stem)
    plt.savefig(out_plot, bbox_inches='tight', dpi=300)
    plt.close()






