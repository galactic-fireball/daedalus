from astropy.convolution import convolve, Tophat2DKernel
import astropy.io.fits as fits
import numpy as np

from jwst.residual_fringe import utils as rfutils


def product_to_path_name(product_name):
    return product_name.replace('_','-').replace('.','-')


def create_smoothed_cube(inst, cube_path):
    prod_name = product_to_path_name(inst.config.product_name)
    new_prod_name = prod_name + '-smoothed'
    out_path = cube_path.parent.joinpath(cube_path.name.replace(prod_name, new_prod_name))

    radius = 1.0
    kernel = Tophat2DKernel(radius)
    kernel.normalize()

    hdul = fits.open(cube_path)
    cube_flux = hdul['SCI'].data
    cube_err = hdul['ERR'].data
    cube_ivar = 1 / cube_err**2

    # cube shape: (nz, ny, nx)
    smoothed_flux = np.empty_like(cube_flux)
    smoothed_ivar = np.empty_like(cube_ivar)

    for i in range(cube_flux.shape[0]):
        smoothed_flux[i] = convolve(cube_flux[i], kernel, boundary='extend', nan_treatment='interpolate', preserve_nan=True)
        smoothed_ivar[i] = convolve(cube_ivar[i], kernel, boundary='extend', nan_treatment='interpolate', preserve_nan=True)

    smoothed_err = 1 / np.sqrt(smoothed_ivar)

    new_hdul = fits.HDUList([h.copy() for h in hdul])
    new_hdul['SCI'].data = smoothed_flux
    new_hdul['ERR'].data = smoothed_err
    new_hdul['SCI'].header['POSTPROC'] = ('Spatially smoothed with Tophat2DKernel, radius = %.02f'%radius)
    new_hdul.writeto(out_path, overwrite=True)
    new_hdul.close()

    hdul.close()

    return out_path


def create_1drfc_cube(cube_path, ch, out_path=None):
    if out_path is None:
        out_path = cube_path.parent.joinpath(cube_path.stem + '-1drfc.fits')

    hdul = fits.open(cube_path)
    cube_flux = hdul['SCI'].data

    header = hdul['SCI'].header
    nwave = hdul['SCI'].data.shape[0]
    wave0 = header['CRVAL3'] - (header['CRPIX3'] - 1) * header['CDELT3']
    wave = (wave0 + np.arange(nwave)*header['CDELT3'])

    rfc_flux = np.full_like(cube_flux, np.nan)
    for x in range(cube_flux.shape[1]):
        print('Row: %d / %d'%(x,cube_flux.shape[1]))
        for y in range(cube_flux.shape[2]):

            spax = cube_flux[:,x,y]
            if np.any(~np.isfinite(spax)):
                continue

            try:
                rfc_flux[:,x,y] = rfutils.fit_residual_fringes_1d(spax, wave, channel=ch)
            except:
                print('Skipping: (%d,%d)'%(x,y))

    new_hdul = fits.HDUList([h.copy() for h in hdul])
    new_hdul['SCI'].data = rfc_flux
    new_hdul['SCI'].header['POSTPROC'] = ('Performed spaxel-by-spaxel 1D residual fringe correction')
    new_hdul.writeto(out_path, overwrite=True)
    new_hdul.close()
    print('1D RFC cube created at %s'%str(out_path))

    hdul.close()


def rfcorr_post_process(inst):
    cube_name_fmt = product_to_path_name(inst.config.product_name)+'-level3'+'_ch%d-%s_s3d.fits'
    targ_out_dir = inst.config.output_dir

    for ch in range(1,5):
        for band in ['short', 'medium', 'long']:
            cube_path = targ_out_dir.joinpath(cube_name_fmt%(ch,band))
            smoothed_path = create_smoothed_cube(inst, cube_path)

            prod_name = product_to_path_name(inst.config.product_name)
            new_prod_name = prod_name + '-1drfc'
            out_path = cube_path.parent.joinpath(cube_path.name.replace(prod_name, new_prod_name))
            cube_1drfc_path = create_1drfc_cube(smoothed_path, ch, out_path=out_path)

            print('Completed 1D RFC for ch%d-%s: %s'%(ch,band,str(cube_1drfc_path)))

