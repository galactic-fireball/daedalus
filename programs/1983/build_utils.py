import astropy.io.fits as fits
import astropy.units as u
from datetime import datetime
import json
import jwst
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
import os
import pathlib

from synphot import SourceSpectrum, units
from synphot.models import Empirical1D

INST_DATA_DIR = pathlib.Path(__file__).resolve().parent.parent.parent.joinpath('instruments', 'data')
WEBBPSF_DIR = INST_DATA_DIR.joinpath('webbpsf', 'webbpsf-data')
os.environ['WEBBPSF_PATH'] = str(WEBBPSF_DIR)

import webbpsf


PIPELINE_DIR = pathlib.Path(__file__).resolve().parent.joinpath('data_sets', 'J1601', 'nirspec_ifu', 'pipeline_1.10.2', 'pipeline')
# PIPELINE_DIR = pathlib.Path(__file__).resolve().parent.joinpath('../1335/data_sets/XID2028/nirspec_ifu/pipeline_1.10.2/pipeline')

TEST_CUBE = PIPELINE_DIR.parent.joinpath('badass', 'jw01983-o001_t001_nirspec_g235m-f170lp_s3d.fits')
PLOT_DIR = TEST_CUBE.parent.joinpath('plots')
PLOT_DIR.mkdir(parents=True, exist_ok=True)
REDSHIFT = 0.03055706


PROGRAM = 1983
OBS_ID = 1
TARGET_ID = 1
PATT_TYPE = '2-point-nod|4-point-nod|along-slit-nod'
GRATING = 'g395m'
FILTER = 'f290lp'


TARGET_FLUX_UNIT = u.erg / u.s / (u.cm**2) / u.um


def check_spec2_asns():
    asns = PIPELINE_DIR.glob('*_spec2_*')
    for asn in asns:
        print(asn.name)
        d = json.load(asn.open('r'))
        # breakpoint()
        print('\tasn id: %s' % d['asn_id'])
        print('\ttarget: %s' % d['target'])
        print('\tproduct: %s' % d['products'][0]['name'])
        cons = d['constraints']
        detector = cons.split('[\'detector\']')[1][12:16]
        print('\tdetector: %s' % detector)
        filt = cons.split('[\'filter\']')[1][12:18]
        print('\tfilt: %s' % filt)
        grating = cons.split('[\'grating\']')[1][12:17]
        print('\tgrating: %s' % grating)
        is_imprt = cons.split('[\'is_imprt\']')[1][11:12]
        print('\tis_imprt: %s' % is_imprt)

        if grating == 'g395m':
            print(d)
        print('\n')


CONSTRAINTS_FMT = '''
DMSAttrConstraint({'name': 'program', 'sources': ['program'], 'value': '%d'})
DMSAttrConstraint({'name': 'instrument', 'sources': ['instrume'], 'value': 'nirspec'})
Constraint_Target({'name': 'target', 'sources': ['targetid'], 'value': '%d'})
DMSAttrConstraint({'name': 'exp_type', 'sources': ['exp_type'], 'value': 'nrs_ifu'})
DMSAttrConstraint({'name': 'patttype', 'sources': ['patttype'], 'value': '%s'})
DMSAttrConstraint({'name': 'opt_elem', 'sources': ['grating'], 'value': '%s'})
Constraint_Obsnum({'name': 'obs_num', 'sources': ['obs_num'], 'value': None})
Constraint_TargetAcq({'name': 'target_acq', 'value': 'target_acquisition'})
DMSAttrConstraint({'name': 'acq_obsnum', 'sources': ['obs_num'], 'value': <function AsnMixin_Science.__init__.<locals>.<lambda> at 0x2b0dbeae3c10>})
DMSAttrConstraint({'name': 'asn_candidate', 'sources': ['asn_candidate'], 'value': \"\\\\('o%03d',\\\\ 'observation'\\\\)\"})
'''

def build_spec3_asn():
    date = datetime.utcnow().strftime('%Y%m%dt%H%M%S')

    asn_dict = {}
    asn_dict['asn_type'] = 'spec3'
    asn_dict['asn_rule'] = 'candidate_Asn_Lv3NRSIFU'
    asn_dict['code_version'] = jwst.__version__
    asn_dict['version_id'] = date
    asn_dict['asn_pool'] = 'jw01983_20230602t094127_pool.csv'
    asn_dict['degraded_status'] = 'No known degraded exposures in association.'
    asn_dict['program'] = '%05d'%PROGRAM
    asn_dict['constraints'] = CONSTRAINTS_FMT % (PROGRAM, TARGET_ID, PATT_TYPE, GRATING, OBS_ID)
    asn_dict['asn_id'] = 'o%03d'%OBS_ID
    asn_dict['target'] = 't%03d'%TARGET_ID

    name = 'jw%05d-o%03d_t%03d_nirspec_%s' % (PROGRAM, OBS_ID, TARGET_ID, GRATING)
    members = []
    cal_files = PIPELINE_DIR.glob('jw01983001001_05101_*_cal.fits')
    for cal in cal_files:
        mem_dict = {}
        mem_dict['expname'] = cal.name
        mem_dict['exptype'] = 'science'
        mem_dict['exposerr'] = 'null'
        mem_dict['asn_candidate'] = '[(\'o%03d\', \'observation\')]' % OBS_ID
        members.append(mem_dict)

    asn_dict['products'] = [{'name':name, 'members':members}]

    out_file = PIPELINE_DIR.joinpath('jw%05d-o%03d_%s_spec3_00001_asn.json'%(PROGRAM, OBS_ID, date))
    json.dump(asn_dict, out_file.open('w'), indent=4)



lines = '''
SiX_1,1.43008,um,351,Cann+18
FeV_1h,1.47688,um,,cloudy
FeIV_1d,1.50629,um,,cloudy
FeIV_1e,1.507,um,,cloudy
FeVI_1j,1.5138,um,,cloudy
FeIV_1f,1.54122,um,,cloudy
FeIV_1g,1.54322,um,,cloudy
FeIV_1h,1.54429,um,,cloudy
FeV_1i,1.62393,um,,cloudy
FeIV_1i,1.66239,um,,cloudy
FeV_1j,1.81848,um,,cloudy
PaA_1,1.87510,um,14,Satyapal+20
FeIV_1j,1.9074,um,,cloudy
FeV_1k,1.91721,um,,cloudy
FeV_1l,1.921,um,,cloudy
SiXI_1,1.93446,um,523,Cann+18
BrD_1,1.94454,um,14,Satyapal+20
SiVI_1,1.96247,um,167,Cann+18
FeV_2,2.03151,um,,cloudy
FeV_2a,2.03959,um,,cloudy
AlIX_2,2.04444,um,285,cloudy
FeV_2b,2.04669,um,,cloudy
FeIV_2,2.08829,um,,cloudy
FeIV_2a,2.12139,um,,cloudy
Br_G,2.16551,um,14,Satyapal+20
FeIV_2b,2.18716,um,,cloudy
CaVIII_2,2.32117,um,127,Cann+18
SiVII_2,2.48071,um,205,cloudy
SiIX_2,2.580,um,304,Satyapal+20
FeVII_2,2.62872,um,,cloudy
FeIV_2c,2.71341,um,,cloudy
FeIV_2d,2.71569,um,,cloudy
FeIV_2e,2.77325,um,,cloudy
FeIV_2f,2.77563,um,,cloudy
FeIV_2g,2.80554,um,31,Satyapal+20
FeIV_2h,2.83081,um,,cloudy
FeIV_2i,2.83329,um,,cloudy
FeIV_2j,2.83562,um,31,Satyapal+20
FeIV_2k,2.86447,um,31,Satyapal+20
FeV_2c,2.89265,um,,cloudy
FeIV_2l,2.90104,um,,cloudy
AlV_2,2.9045,um,120,Cann+18
FeIV_2m,2.90821,um,,cloudy
FeIV_2n,2.97281,um,,cloudy
FeV_3,3.02517,um,,cloudy
MgVIII,3.027950,um,225,Satyapal+20
FeIV_3,3.03352,um,,cloudy
FeV_3a,3.15073,um,,cloudy
FeIV_3a,3.16981,um,,cloudy
KVII_3,3.18966,um,,cloudy
CaIV_3,3.2061,um,51,Cann+18
FeIV_3b,3.21932,um,,cloudy
'''


lines = '''
FeV_3b,3.27093,um,,cloudy
FeV_3c,3.36066,um,,cloudy
FeVII_3,3.38363,um,,cloudy
FeIV_3c,3.38592,um,,cloudy
FeIV_3d,3.39109,um,,cloudy
FeV_3d,3.51633,um,,cloudy
FeIV_3e,3.58762,um,,cloudy
FeIV_3f,3.58878,um,,cloudy
AlVI_3,3.65932,um,154,cloudy
AlVIII_3,3.69722,um,,cloudy
SIX_3,3.75414,um,,cloudy
FeIV_3g,3.8082,um,,cloudy
SiIX_3,3.9282,um,304,cloudy
MgV_3,3.9656,um,,cloudy
DI_4,4.01974,um,14,Satyapal+20
BrA_4,4.05113,um,14,Satyapal+20
FeVI_4,4.1431,um,,cloudy
CaV_4,4.15739,um,,cloudy
FeVI_4a,4.40526,um,,cloudy
MgIV_4,4.48712,um,80,Cann+18
ArVI_4,4.528,um,75,Cann+18
KIII_4,4.61683,um,32,Satyapal+20
NaVII_4,4.68257,um,,cloudy
FeVI_4b,4.68828,um,,cloudy
FeIV_4,4.72617,um,,cloudy
FeII_4,4.889,um,7.9,Pereira-Santaella+22
FeV_4,4.89702,um,,cloudy
ArV_4,4.92778,um,,cloudy
H2S8_5,5.053,um,,Pereira-Santaella+22
H2O3,5.17942,um,15.37,Kameron
FeV_5,5.23472,um,,cloudy
FeII_5,5.33881,um,8,Satyapal+20
FeIV_5,5.39606,um,,cloudy
FeVIII_5,5.447,um,125,Pereira-Santaella+22
H2O9,5.448,um,15.37,Kameron
MgVII_5,5.50177,um,187,Cann+18
H2S7_5,5.511,um,15.37,Pereira-Santaella+22;Kameron
KVI_5,5.57324,um,,cloudy
MgV_5,5.607,um,109,Cann+18
FeV_5a,5.71894,um,,cloudy
AlVIII_5,5.82933,um,,cloudy
KIV_5,5.9803,um,,cloudy
'''


def create_line_list():

    for line in lines.split('\n'):
        if line == '':
            continue

        name, wave, _, _, _ = line.split(',')
        print('\'na_%s\': {\'center\':%.3f, \'line_type\':\'user\', \'line_profile\':\'gaussian\'},' % (name, (wave*u.um).to(u.AA)))





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

    return dered(wave, z), np.nan_to_num(spec)



def get_webbpsf2():

    wave, spec = get_cube_data(TEST_CUBE)

    nrs = webbpsf.NIRSpec()
    nrs.image_mask = None # No MSA for IFU mode
    # nl = np.linspace(2.87e-6, 5.27e-6, 6)
    nl = wave.to(u.m).value

    # Calculate PSF datacube
    breakpoint()
    cube = nrs.calc_datacube(wavelengths=nl[100:120], fov_pixels=27, oversample=4)

    # Display the contents of the data cube
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(10,7))
    for iy in range(2):
        for ix in range(3):
            ax=axes[iy,ix]
            i = iy*3+ix
            wl = cube[0].header['WAVELN{:02d}'.format(i)]

            # Note that when displaying datacubes, you have to set the "cube_slice" parameter
            webbpsf.display_psf(cube, ax=ax, cube_slice=i,
                                title="NIRSpec, $\lambda$ = {:.3f} $\mu$m".format(wl*1e6),
                                vmax=.2, vmin=1e-4, ext=1, colorbar=False)
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
    plt.show()
    breakpoint()


def clip_spec(spec, fov=30):
    lenx, leny, _ = spec.shape
    clipxs = int(np.ceil((lenx - fov) / 2))
    clipxe = lenx - clipxs - fov
    clipys = int(np.ceil((leny - fov) / 2))
    clipye = leny - clipys - fov
    return spec[clipxs:-clipxe,clipys:-clipye,:]



def get_webbpsf():
    wave, spec = get_cube_data(TEST_CUBE)

    nrs = webbpsf.NIRSpec()
    nrs.image_mask = 'IFU'
    # nrs.filter = 'IFU'
    waves = wave.to(u.m).value

    # wave_a = wave.to(u.AA).value

    # wave = [1000, 2000, 3000, 4000, 5000]  # Angstrom
    # flux = [1e-17, -2.3e-18, 1.8e-17, 4.5e-17, 9e-18] * units.FLAM
    # sp = SourceSpectrum(Empirical1D, points=wave, lookup_table=flux)

    # sp = SourceSpectrum(Empirical1D, points=wave_a[100:104], lookup_table=spec[:,:,100:104])
    # breakpoint()

    # breakpoint()
    # outfits = nrs.calc_psf(source={'wavelengths':waves[100:104]})
    # breakpoint()

    psfcube = nrs.calc_datacube(waves, fov_pixels=30, oversample=4, add_distortion=True)
    outfile = TEST_CUBE.parent.joinpath('webbpsf_ifucube.fits')
    psfcube.writeto(outfile)




def plot_cube_sum(wave, spec, outsuffix=''):
    fnu_sum = np.sum(spec, axis=(0, 1))
    cube_sum = np.sum(spec, axis=0)
    breakpoint()

    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5)) 
    ax1.plot(wave, fnu_sum) 
    # ax1.set_xlim(0.95, 1.5)
    ax1.set_title('Spaxel sums')
    ax1.set_xlabel('Wavelength (um)')  
    ax1.set_ylabel('Flux')

    ax2.imshow(cube_sum, norm=LogNorm())
    ax2.set_title('Slice sums')

    plt.savefig(PLOT_DIR.joinpath('cube_sum%s.png'%outsuffix))


def test_psf():
    wave, spec = get_cube_data(TEST_CUBE)
    spec = clip_spec(spec)
    plot_cube_sum(wave, spec)


def main():
    test_psf()






if __name__ == '__main__':
    main()
