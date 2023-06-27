import astropy.io.fits as fits
import astropy.units as u
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
import pathlib


PLOT_DIR = pathlib.Path(__file__).resolve().parent.joinpath('plots')
PLOT_DIR.mkdir(parents=True, exist_ok=True)
TARGET_FLUX_UNIT = u.erg / u.s / (u.cm**2) / u.um
PA_ALPHA_WAVE = 1.87510 * u.um
DEFAULT_NEIGHBOR_RAD = 1


def dered(wave, z):
    return wave / (1 + z)


def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx


def get_cube_data(cube_fits, z=0.0):
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

    tf = (hdu['SCI'].data.T * u.Unit(header['BUNIT']))
    spec = tf.to(TARGET_FLUX_UNIT, equivalencies=u.spectral_density(wave))

    return dered(wave, z), spec


def clip_spec(spec, fov=30):
    lenx, leny, _ = spec.shape
    clipxs = int(np.ceil((lenx - fov) / 2))
    clipxe = lenx - clipxs - fov
    clipys = int(np.ceil((leny - fov) / 2))
    clipye = leny - clipys - fov
    return spec[clipxs:-clipxe,clipys:-clipye,:]


def get_neighbor_med(spec, i, j, w, r):
    if np.isnan(spec[:,:,w]).all():
        if w+1 == spec.shape[2]:
            return 1e-20*TARGET_FLUX_UNIT # we've run out of options
        return get_neighbor_med(spec, i, j, w+1, r)

    maxi = spec.shape[0]-1
    maxj = spec.shape[1]-1

    old_val = spec[i,j,w]
    spec[i,j,w] = np.nan # set to NaN here so we can ignore it
    med_value = np.nanmedian(spec[max(0,i-r):min(maxi,i+r+1), max(0,j-r):min(maxj,j+r+1), w])
    spec[i,j,w] = old_val
    if np.isnan(med_value): # possible if the whole box is NaN
        return get_neighbor_med(spec, i, j, w, r+1)
    return med_value


def set_spec_invalid(spec):
    for i in range(spec.shape[0]):
        for j in range(spec.shape[1]):
            for w in range(spec.shape[2]):
                if np.isnan(spec[i,j,w]) or spec[i,j,w] < 0:
                    # print('Invalid pixel (%d, %d) @ %d' % (i, j, w))
                    spec[i,j,w] = get_neighbor_med(spec, i, j, w, DEFAULT_NEIGHBOR_RAD)


def flux_limit(wave, spec, line):
    LINE_PAD = 0.02*u.um
    for i in range(spec.shape[0]):
        for j in range(spec.shape[1]):
            _, min_idx = find_nearest(wave, line-LINE_PAD)
            _, max_idx = find_nearest(wave, line+LINE_PAD)
            max_value = np.nanmax(spec[i,j,min_idx:max_idx])
            # if i == 22 and j == 34:
            #   breakpoint()
            for w in range(spec.shape[2]):
                if spec[i,j,w] > max_value:
                    print('Hot pixel (%d, %d) @ %d' % (i, j, w))
                    spec[i,j,w] = get_neighbor_med(spec, i, j, w, DEFAULT_NEIGHBOR_RAD)


def post_process(wave, spec):
    set_spec_invalid(spec)
    flux_limit(wave, spec, PA_ALPHA_WAVE)
    return np.nan_to_num(spec, copy=False)


def plot_cube_sum(wave, spec, outsuffix=''):
    fnu_sum = np.sum(spec, axis=(0, 1))
    cube_sum = np.sum(spec, axis=2)

    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    ax1.plot(wave, fnu_sum)
    # ax1.set_xlim(0.95, 1.5)
    ax1.set_title('Spaxel sums')
    ax1.set_xlabel('Wavelength (um)')
    ax1.set_ylabel('Flux')

    ax2.imshow(cube_sum, norm=LogNorm())
    ax2.set_title('Slice sums')

    plt.savefig(PLOT_DIR.joinpath('cube_sum%s.png'%outsuffix))


def plot_circular_aperture(wave, spec, center_xy=None, outsuffix=''):
    r_pix = 5.92

    if not center_xy:
        center_xy = [spec.shape[0]/2, spec.shape[1]/2]

    aperture = CircularAperture(center_xy, r=r_pix)

    cylinder_sum = []
    for slice2d in data:
        phot_table = aperture_photometry(slice2d, aperture)
        cylinder_sum.append(phot_table['aperture_sum'][0])


def plot_cone_aperture():
    lambda0 = wavelength[0]
    print('Reference wavelength:', lambda0)

    cone_sum = []
    idx = -1
    for (slice2d, wave) in zip(data, wavelength):
        idx = idx + 1
        r_cone = r_pix * wave / lambda0
        aperture_cone = CircularAperture(center_xy, r=r_cone)
        phot_table = aperture_photometry(slice2d, aperture_cone, wcs=w.celestial, method='exact')
        cone_sum.append(phot_table['aperture_sum'][0])


def test_spaxel_spec():
    test_cube = pathlib.Path(__file__).parent.joinpath('data_sets/J1601/nirspec_ifu/pipeline_1.10.2/badass/jw01983-o001_t001_nirspec_g235m-f170lp_s3d.fits').resolve()
    z = 0.03055706

    xsub = 13
    ysub = 12

    test_spaxels = [(28,31), (28,30), (28,29), (28,27), (28,26), (30, 25)]
    # test_spaxels = [(25, 35), (22, 34), (23,18)]
    # test_spaxels = [(22, 34), (23,18)]
    nspax = len(test_spaxels)

    wave, spec = get_cube_data(test_cube, z=z)
    spec = clip_spec(spec)
    fig, axs = plt.subplots(nspax, 2, figsize=(20,10))

    def plot_spaxels(col):
        for i in range(nspax):
            spax = test_spaxels[i]
            axs[i,col].plot(wave, spec[spax[0]-xsub,spax[1]-ysub,:])
            axs[i,col].set_title('Spaxel (%d, %d)' % (spax[0], spax[1]))

    plot_spaxels(0)
    post_process(wave, spec)
    plot_spaxels(1)
    fig.savefig(PLOT_DIR.joinpath('test_spaxels.png'))

    plot_cube_sum(wave, spec)


def main():
    test_spaxel_spec()


if __name__ == '__main__':
    main()
