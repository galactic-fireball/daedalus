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

print('Using jwst pipeline version: %s' % jwst.__version__)


USE_MULTIPROCESS = True
NPROCESSES = 4



def run_stage1_single(ufile, output_dir):
	print('Processing: {}'.format(str(ufile)))
	if output_dir.joinpath(ufile.name.replace('uncal', 'rate')).exists():
		return None

	detector1 = Detector1Pipeline()
	detector1.output_dir = str(output_dir)
	detector1.output_file = str(output_dir.joinpath(ufile.stem))
	detector1.run(ufile)
	return None


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


def update_asn(asn_temp):
	data = json.load(open(asn_temp))

	name = data['products'][0]['name']
	data['code_version'] = jwst.__version__
	data['version_id'] = datetime.utcnow().strftime('%Y%m%dt%H%M%S')
	data['asn_pool'] = ''

	asn_file = asn_temp.parent.joinpath('%s_updated.json' % asn_temp.stem)
	json.dump(data, open(asn_file, 'w'))

	return asn_file, name


def run_stage2_single(asn_temp, output_dir):
	print('Processing: {}'.format(str(asn_temp)))
	asn_file, name = update_asn(asn_temp)

	if output_dir.joinpath('%s_cal.fits'%name).exists():
		return None

	crds_config = Spec2Pipeline.get_config_from_reference(str(asn_file))
	spec2 = Spec2Pipeline.from_config_section(crds_config)

	spec2.output_dir = str(output_dir)
	spec2.save_results = True
	spec2.run(asn_file)
	return None



def run_stage2_all(indir, output_dir):
	asn_files = indir.glob('*_spec2_*_asn.json')

	if not USE_MULTIPROCESS:
		for asn_file in asn_files:
			run_stage2_single(asn_file, output_dir)
		return

	args = [(asn_file, output_dir) for asn_file in asn_files]
	pool = mp.Pool(processes=NPROCESSES, maxtasksperchild=1)
	pool.starmap(run_stage2_single, args, chunksize=1)
	pool.close()
	pool.join()


def run_stage3_single(asn_temp, out_dir):
	print('Processing: {}'.format(str(asn_temp)))
	asn_file, name = update_asn(asn_temp)

	if len(list(out_dir.glob('%s*_s3d.fits'%name))) == 1:
		return None

	crds_config = Spec3Pipeline.get_config_from_reference(str(asn_file))
	spec3 = Spec3Pipeline.from_config_section(crds_config)

	spec3.output_dir = str(out_dir)
	spec3.save_results = True
	spec3.run(asn_file)
	return None


def run_stage3_all(indir, output_dir):
	asn_files = indir.glob('*_spec3_*_asn.json')

	if not USE_MULTIPROCESS:
		for asn_file in asn_files:
			run_stage3_single(asn_file, output_dir)
		return

	args = [(asn_file, output_dir) for asn_file in asn_files]
	# pool = mp.Pool(processes=NPROCESSES, maxtasksperchild=1)
	# pool.starmap(run_stage3_single, args, chunksize=1)
	pool = mp.Pool()
	pool.starmap(run_stage3_single, args)
	pool.close()
	pool.join()
