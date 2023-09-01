import argparse
from astroquery.mast import Mast, Observations
from astropy.table import unique, vstack
import os
import pandas as pd
import pathlib
import sys

INSTRUMENTS = {'miri': 'MIRI', 'miri_ifu':'MIRI/IFU', 'nirspec': 'NIRSPEC', 'nirspec_ifu': 'NIRSPEC/IFU'}
PRODUCTS = {'asn': 'ASN', 'cal': 'CAL', 'uncal': 'UNCAL'}

# ['intentType', 'obs_collection', 'provenance_name', 'instrument_name', 'project',
# 'filters', 'wavelength_region', 'target_name', 'target_classification', 'obs_id', 's_ra', 's_dec',
# 'dataproduct_type', 'proposal_pi', 'calib_level', 't_min', 't_max', 't_exptime', 'em_min', 'em_max',
# 'obs_title', 't_obs_release', 'proposal_id', 'proposal_type', 'sequence_number', 's_region', 'jpegURL',
# 'dataURL', 'dataRights', 'mtFlag', 'srcDen', 'obsid', 'objID']

# array(['MIRI/IMAGE', 'NIRISS/WFSS', 'NIRCAM/IMAGE', 'NIRSPEC/MSA',
#        'NIRCAM/GRISM', 'MIRI/IFU', 'MIRI/SLIT', 'NIRSPEC', 'NIRCAM/CORON',
#        'NIRISS/IMAGE', 'NIRSPEC/IFU', 'MIRI/TARGACQ', 'NIRSPEC/SLIT',
#        'FGS/FGS2', 'NIRISS/AMI', 'MIRI/CORON', 'FGS/FGS1', 'NIRISS/SOSS',
#        'MIRI/SLITLESS', 'NIRISS', 'MIRI', 'NIRCAM', 'NIRCAM/TARGACQ',
#        'NIRSPEC/IMAGE', 'FGS'], dtype=object)

logged_in = False
def login(mast_api_token=None):
	global logged_in
	if logged_in:
		return

	auth_token = os.environ.get('MAST_API_TOKEN', mast_api_token)
	if auth_token == '':
		raise Exception('Must provide a MAST API token')

	Observations.login(auth_token)
	logged_in = True


def get_full_jwst_table(cache_file=None):
    if cache_file and cache_file.exists():
        return pd.read_csv(cache_file)

    df = Observations.query_criteria(obs_collection=['JWST']).to_pandas()
    if cache_file:
        df.to_csv(cache_file, index=False)
    return df


def get_all_service_data(service, cache_file=None):
	if cache_file and cache_file.exists():
		return pd.read_csv(cache_file)

	parameters = {'columns': '*', 'filters': []}
	response = Mast.service_request_async(service, parameters)
	results = response[0].json()['data']
	df = pd.DataFrame(results)

	if cache_file:
		df.to_csv(cache_file, index=False)

	return df


def get_all_miri_data(cache_file=None):
	return get_all_service_data('Mast.Jwst.Filtered.Miri', cache_file=cache_file)


def get_all_nirspec_data(cache_file=None):
	return get_all_service_data('Mast.Jwst.Filtered.Nirspec', cache_file=cache_file)


def get_jwst_data(args, expand=False, filt=None):
	args['obs_collection'] = ['JWST']
	obs_list = Observations.query_criteria(**args)
	if expand:
		obs_list = [Observations.get_product_list(obs) for obs in obs_list]
		obs_list = unique(vstack(obs_list), keys='productFilename')

	if filt:
		obs_list = Observations.filter_products(obs_list, **filt)

	return obs_list

def get_program_data(program_id, instrument_name):
	return get_jwst_data({'proposal_id':program_id, 'instrument_name':instrument_name}, expand=True)
	# obs_list = Observations.query_criteria(obs_collection=['JWST'], proposal_id=program_id, instrument_name=instrument_name)
	# all_products = [Observations.get_product_list(obs) for obs in obs_list]
	# return unique(vstack(all_products), keys='productFilename')


def get_instrument_data(instrument_name):
	return Observations.query_criteria(obs_collection=['JWST'], instrument_name=instrument_name)


def get_data_products(prog_id, inst, calib_level, product_type):
	prog_data = get_program_data(prog_id, inst)
	return Observations.filter_products(prog_data, calib_level=[calib_level], productSubGroupDescription=product_type)


def download_file(file_name, dest=None):
	uri = f'mast:jwst/product/%s' % file_name
	result = Observations.download_file(uri, base_url='https://mast.stsci.edu/api/v0.1/Download/file')
	if result[0] == 'ERROR':
		raise RuntimeError('Error retrieving file: ' + result[1])

	if dest:
		pathlib.Path(file_name).rename(dest)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--download', help='file to download', type=str, default=None)
	parser.add_argument('--query', help='query JWST', action='store_true', default=False)
	parser.add_argument('--all', help='expand JWST query to all observations', action='store_true', default=False)
	parser.add_argument('--instrument', help='intrument to query', choices=INSTRUMENTS.keys(), default=None)
	parser.add_argument('--level', help='calibration level', type=int, choices=[1,2,3], default=None)
	parser.add_argument('--program', help='program id', type=str, default=None)
	parser.add_argument('--product', help='product type', choices=PRODUCTS.keys(), default=None)
	parser.add_argument('--public', help='get publicly available only', action='store_true', default=False)
	parser.add_argument('--out', help='output file', type=str, default=None)
	args = parser.parse_args()

	if args.download:
		download_file(args.download)
		return

	if args.query:
		query = {}
		if args.instrument:
			query['instrument_name'] = INSTRUMENTS[args.instrument]
		if args.program:
			query['proposal_id'] = args.program

		filt = {}
		if args.level:
			filt['calib_level'] = [args.level]
		if args.product:
			filt['productSubGroupDescription'] = PRODUCTS[args.product]
		if args.public:
			filt['dataRights'] = 'PUBLIC'

		obs_list = get_jwst_data(query, expand=args.all, filt=filt)
		if args.out:
			outfile = pathlib.Path(args.out).resolve()
			obs_list.to_pandas().to_csv(outfile, index=False)
		return


if __name__ == '__main__':
	main()
