import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np
import pathlib
from photutils.aperture import CircularAperture, RectangularAperture
from photutils.segmentation import detect_sources
from photutils.segmentation import SourceCatalog
from skimage.draw import disk
import sys

from jwst.datamodels import dqflags

INSTRUMENTS_DIR = pathlib.Path(__file__).parent
sys.path.insert(0, str(INSTRUMENTS_DIR))
import nsclean as nc # TODO: check for failure and mark rate clean unavailable

# TODO: make toml options

# Snowball detection options
THRESH_SIGMA = 3.0
BKG = 0.0
STD_DEV = 0.00007
NPIX_DET = 50 # The minimum number of connected pixels to detect a snowball

# Thresholds for determining if a snowball candidate is real
ECC_THRESH = 0.6 # eccentricity
BOX_RATIO_THRESH = 1.5 # bounding box width/height ratio threshold

HALO_RAD_FACT = 2 # Factor to increase the radius of the snowball to cut
SAT_RADIUS = 3 # Additional saturation radius to cut

PLOT = True

# Only a slightly modified version of Chris Willot's dosnowballflags.py:
# https://github.com/chriswillott/jwst/blob/master/dosnowballflags.py

def flag_snowballs(ramp_file):

    plot_dir = ramp_file.parent.joinpath('snowball_plots')
    plot_dir.mkdir(parents=True, exist_ok=True)

    hdu = fits.open(ramp_file)

    groupdq = hdu['GROUPDQ'].data
    header = hdu[0].header
    nints = header['NINTS']
    ngroups = header['NGROUPS']

    # iterate over integrations
    for integ in range(nints):
        snow_cnt = 0
        # Skip first group because no jumps there
        for group in range(1, ngroups):

            # Find pixels in this group with jump detected and/or saturation detected
            jumps = np.zeros((2048,2048), dtype='uint8')
            sat = np.zeros((2048,2048), dtype='uint8')
            jumps_sat = np.zeros((2048,2048), dtype='uint8')
            cur_gdq = np.squeeze(groupdq[integ,group,:,:])

            i_yy,i_xx = np.where(np.bitwise_and(cur_gdq, dqflags.group['JUMP_DET']) != 0)
            jumps[i_yy,i_xx] = 1
            jumps_sat[i_yy,i_xx] = 1

            i_yy,i_xx = np.where(np.bitwise_and(cur_gdq, dqflags.group['SATURATED']) != 0)
            sat[i_yy,i_xx] = 1
            jumps_sat[i_yy,i_xx] = 1

            det_threshold = BKG + (THRESH_SIGMA * STD_DEV)
            # Run initial find on jumps or saturated because some short ramps do not have a jump in the regions that saturate
            detects = detect_sources(jumps_sat, det_threshold, npixels=NPIX_DET)
            det_image = detects.data.astype(np.uint32)

            if np.max(det_image) <= 0:
                continue

            if PLOT:
                plt.figure()
                plt.imshow(det_image.data)

            det_tbl = SourceCatalog(jumps, detects).to_table()

            ctsnowballs = det_tbl['xcentroid'][:].size
            # iterate over each possible snowball
            for k in range(ctsnowballs): # TODO: can just iterate over each row of the table ?

                # If low eccentricity proceed, otherwise remove source from segmentation image
                # Use both eccentricity and segmentation box axis since not always consistent, e.g. for merged jumps
                segboxaxisratio = np.abs((det_tbl['bbox_xmax'][k]-det_tbl['bbox_xmin'][k])/(det_tbl['bbox_ymax'][k]-det_tbl['bbox_ymin'][k]))
                if segboxaxisratio < 1.0:
                    # TODO: segboxaxisratio = 1 / segboxaxisratio ?
                    segboxaxisratio = np.abs((det_tbl['bbox_ymax'][k]-det_tbl['bbox_ymin'][k])/(det_tbl['bbox_xmax'][k]-det_tbl['bbox_xmin'][k]))

                if ((det_tbl['eccentricity'][k] >= ECC_THRESH) | (segboxaxisratio >= BOX_RATIO_THRESH)):
                    continue

                if PLOT:
                    xcent = (det_tbl['bbox_xmax'][k]-det_tbl['bbox_xmin'][k])/2+det_tbl['bbox_xmin'][k]
                    ycent = (det_tbl['bbox_ymax'][k]-det_tbl['bbox_ymin'][k])/2+det_tbl['bbox_ymin'][k]
                    width = (det_tbl['bbox_xmax'][k]-det_tbl['bbox_xmin'][k])
                    height = (det_tbl['bbox_ymax'][k]-det_tbl['bbox_ymin'][k])
                    ap = RectangularAperture((xcent,ycent), width, height)
                    ap.plot()

                xmin, xmax, ymin, ymax = int(det_tbl['bbox_xmin'][k]), int(det_tbl['bbox_xmax'][k]), int(det_tbl['bbox_ymin'][k]), int(det_tbl['bbox_ymax'][k])
                jumpscutout = jumps[ymin:ymax,xmin:xmax]
                satcutout = sat[ymin:ymax,xmin:xmax]
                jumps_sat_cutout = jumps_sat[ymin:ymax,xmin:xmax]

                # Triple box size for increased area to flag further out
                bigoffsetx = int((det_tbl['bbox_xmax'][k]-int(det_tbl['bbox_xmin'][k])))
                bigoffsety = int((det_tbl['bbox_ymax'][k]-int(det_tbl['bbox_ymin'][k])))

                bigsizex = bigoffsetx * 3
                bigsizey = bigoffsety * 3

                jumpsbigcutout = np.zeros((bigsizey,bigsizex), dtype=np.uint8)
                jumpsbigcutout[bigoffsety:(bigoffsety+jumpscutout.shape[0]), bigoffsetx:(bigoffsetx+jumpscutout.shape[1])] = jumpscutout

                satbigcutout = np.zeros((bigsizey,bigsizex), dtype=np.uint8)
                satbigcutout[bigoffsety:(bigoffsety+jumpscutout.shape[0]), bigoffsetx:(bigoffsetx+jumpscutout.shape[1])] = satcutout

                # For jumps assume round and use all jump or saturated pixels to get area
                num_jumps_sat = jumps_sat_cutout[np.where(jumps_sat_cutout > 0)].size
                radiusjumporsat = np.sqrt(num_jumps_sat/np.pi)
                radius = int(HALO_RAD_FACT*radiusjumporsat)
                rr, cc = disk((bigsizey/2-0.5,bigsizex/2-0.5), radius) # TODO: why -0.5?
                jumpsbigcutout[rr, cc] = 4 # TODO: arbitrary?

                #For saturation assume round and use saturated pixels to get area
                numsat = satcutout[np.where(satcutout>0)].size
                radiussat = np.sqrt(numsat/np.pi)
                radius = int(radiussat+SAT_RADIUS)
                rr, cc = disk((bigsizey/2-0.5,bigsizex/2-0.5), radius) # TODO: why -0.5?
                satbigcutout[rr, cc] = 2 # TODO: arbitrary?

                xlo = xmin-bigoffsetx
                ylo = ymin-bigoffsety

                #Update pixels in GROUPDQ array for halo
                i_yy,i_xx = np.where(jumpsbigcutout > 0)
                i_yy += ylo
                i_xx += xlo
                for l in range(len(i_xx)): # TODO: better loop init
                    if ((i_xx[l] > 3) & (i_xx[l] < 2044) & (i_yy[l] > 3) & (i_yy[l] < 2044)): # TODO: why?
                        groupdq[integ,group,i_yy[l],i_xx[l]] = np.bitwise_or(groupdq[integ,group,i_yy[l],i_xx[l]], dqflags.group['JUMP_DET'])
                        # if the snowball happened in the third group, flag the second group similarly in case first effects happened there and not enough data for good ramps up to there anyway
                        if group == 2:
                            groupdq[integ,group-1,i_yy[l],i_xx[l]] = np.bitwise_or(groupdq[integ,group-1,i_yy[l],i_xx[l]], dqflags.group['JUMP_DET'])

                #Update pixels in GROUPDQ array for saturated core
                i_yy,i_xx = np.where(satbigcutout > 0)
                i_yy += ylo
                i_xx += xlo
                for l in range(len(i_xx)): # TODO: better loop init
                    if ((i_xx[l] > 3) & (i_xx[l] < 2044) & (i_yy[l] > 3) & (i_yy[l] < 2044)): # TODO: why?
                        groupdq[integ,group,i_yy[l],i_xx[l]] = np.bitwise_or(groupdq[integ,group,i_yy[l],i_xx[l]], dqflags.group['SATURATED'])
                        # if the snowball happened in the third group, flag the second group similarly in case first effects happened there and not enough data for good ramps up to there anyway
                        if group == 2:
                            groupdq[integ,group-1,i_yy[l],i_xx[l]] = np.bitwise_or(groupdq[integ,group-1,i_yy[l],i_xx[l]], dqflags.group['SATURATED'])

                snow_cnt += 1

            if PLOT:
                plt.tight_layout()
                plt.savefig(plot_dir.joinpath('%s_%d_%d.png'%(ramp_file.stem, integ, group)))
                plt.close()

        # Any pixel flagged as saturated in a group must be flagged as saturated in all subsequent groups  
        for group in range(ngroups):
            # Find pixels in this group with saturation detected
            sat = np.zeros((2048,2048), dtype='uint8')
            cur_gdq = np.squeeze(groupdq[integ,group,:,:])
            i_yy,i_xx = np.where(np.bitwise_and(cur_gdq, dqflags.group['SATURATED']) != 0)
            for l in range(len(i_xx)):
                groupdq[integ,group:,i_yy[l],i_xx[l]] = np.bitwise_or(groupdq[integ,group:,i_yy[l],i_xx[l]], dqflags.group['SATURATED'])

        header['HISTORY'] = 'Corrected {} snowballs in integration {}'.format(snow_cnt,(integ+1))

    hdu['GROUPDQ'].data = groupdq 
    header.set('SNOWCORR', 'COMPLETE', 'snowball flagging applied')
    hdu[0].header = header
    out_file = ramp_file.parent.joinpath(ramp_file.stem+'_snow.fits')
    hdu.writeto(out_file, overwrite=True)
    hdu.close()
    return out_file


def run_nsclean(rate_file):
    unclean_file = rate_file.parent.joinpath(rate_file.stem + '_unclean.fits')
    rate_file.rename(unclean_file)

    hdu = fits.open(unclean_file)
    data = hdu['SCI'].data

    detector = 'NRS%d' % int(unclean_file.stem.split('nrs')[1][0])
    mask_file = INSTRUMENTS_DIR.joinpath('data', 'nsclean', '%s_ifu_mask_thorough.fits'%detector.lower())
    if not mask_file.exists():
        raise Exception('Unable to find NSClean mask file: %s' % str(mask_file))

    mhdu = fits.open(mask_file)
    mask = np.array(mhdu[0].data, dtype=np.bool_)
    mhdu.close()

    cleaner = nc.NSClean(detector, mask)
    clean_data = cleaner.clean(data)

    hdu['SCI'].data = clean_data
    hdu[0].header.set('NSCLEAN', 'COMPLETE', 'Processed by NSClean v' + nc.__version__)
    hdu.writeto(rate_file, overwrite=True)

