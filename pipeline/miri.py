from datetime import datetime
import json
import multiprocessing as mp
import os
import pathlib
import sys

import ssl
# Needed to work around ssl certificate verification during crds downloads
ssl._create_default_https_context = ssl._create_unverified_context

USE_CRDS_OPS = True
# Needs to be set before crds/jwst imports
if USE_CRDS_OPS:
	os.environ['CRDS_PATH'] = str(pathlib.Path(__file__).resolve().parent.joinpath('crds_cache', 'ops'))
	os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu'
else:
	os.environ['CRDS_PATH'] = str(pathlib.Path(__file__).resolve().parent.joinpath('crds_cache', 'pub'))
	os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds-pub.stsci.edu'

import crds
import jwst
from jwst.associations import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from jwst.associations.lib.rules_level3 import Asn_Lv3SpectralSource
from jwst.pipeline.calwebb_detector1 import Detector1Pipeline
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline
from jwst.pipeline.calwebb_spec3 import Spec3Pipeline
from jwst.residual_fringe import ResidualFringeStep

print('Using jwst pipeline version: %s' % jwst.__version__)


RESIDUAL_FRINGE_CORRECT = False
USE_MULTIPROCESS = True
NPROCESSES = 4
INCLUDE_BK = True



def run_stage1_single(ufile, output_dir):
	print('Processing: {}'.format(str(ufile)))
	if not output_dir.joinpath(ufile.name.replace('uncal', 'rate')).exists():
		detector1 = Detector1Pipeline()
		detector1.output_dir = str(output_dir)
		detector1.output_file = str(output_dir.joinpath(ufile.stem))
		detector1.run(ufile)


def run_stage1_all(uncal_dir, output_dir):
	uncal_files = uncal_dir.glob('*_uncal.fits')

	if not USE_MULTIPROCESS:
		for ufile in uncal_files:
			run_stage1_single(ufile, output_dir)
		return

	args = [(ufile, output_dir) for ufile in uncal_files]
	pool = mp.Pool(processes=NPROCESSES, maxtasksperchild=1)
	pool.starmap(run_stage1_single, args, chunksize=1)
	pool.close()
	pool.join()



def run_stage2_single(rfile, output_dir, skip_cubes):

	if not output_dir.joinpath(rfile.name.replace('rate', 'cal')).exists():
		spec2 = Spec2Pipeline()
		spec2.output_dir = str(output_dir)

		print('Processing: {}'.format(str(rfile)))
		spec2.output_file = str(rfile.with_suffix(''))

		if (skip_cubes):
			spec2.cube_build.skip = True
			spec2.extract_1d.skip = True

		spec2.run(rfile)
		spec2.suffix = 'spec2'

	if RESIDUAL_FRINGE_CORRECT:
		# TODO: if skip_cubes == False, does this need to be on the x1d files?
		rfs_in = output_dir.joinpath(rfile.name.replace('rate', 'cal'))
		if not output_dir.joinpath(rfs_in.name.replace('cal', 'residual_fringe')).exists():
			print('Performing residual fringe correction on: {}'.format(str(rfs_in)))
			rfringe_step = ResidualFringeStep()
			rfringe_step.output_file = str(rfs_in.with_suffix(''))
			rfringe_step.output_dir = str(output_dir)
			rfringe_step.run(rfs_in)


def run_stage2_all(input_dir, output_dir, skip_cubes=False):
	rate_files = input_dir.glob('*_rate.fits')

	if not USE_MULTIPROCESS:
		for rfile in rate_files:
			run_stage2_single(rfile, output_dir, skip_cubes)
		return

	args = [(rfile, output_dir, skip_cubes) for rfile in rate_files]
	pool = mp.Pool(processes=NPROCESSES, maxtasksperchild=1)
	pool.starmap(run_stage2_single, args, chunksize=1)
	pool.close()
	pool.join()



def run_stage3(sci_input, bk_input, product_name, output_dir):
	if RESIDUAL_FRINGE_CORRECT:
		sci_files = sci_input.glob('*_residual_fringe.fits')
	else:
		print(str(sci_input))
		sci_files = sci_input.glob('*_cal.fits')
	sci_files = [str(f) for f in sci_files]
	print(sci_files)

	asn = asn_from_list.asn_from_list(sci_files, rule=DMS_Level3_Base, product_name=product_name)

	if INCLUDE_BK:
		bk_files = bk_input.glob('*_x1d.fits')
		for bkfile in bk_files:
			asn['products'][0]['members'].append({'expname': str(bkfile), 'exptype': 'background'})

	asnfile = output_dir.joinpath('%s_asn.json' % product_name)
	with open(asnfile, 'w') as afile:
		afile.write(asn.dump()[1])

	crds_config = Spec3Pipeline.get_config_from_reference(str(asnfile))
	spec3 = Spec3Pipeline.from_config_section(crds_config)

	spec3.output_dir = str(output_dir)
	spec3.save_results = True
	spec3.run(asnfile)
