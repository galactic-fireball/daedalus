import argparse
from astroquery.mast import Observations
from astropy.table import unique, vstack
import os
import pandas as pd
import pathlib
import sys

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


def get_proposal_data(proposal_id, instrument_name):
	return Observations.query_criteria(obs_collection=['JWST'], proposal_id=proposal_id, instrument_name=instrument_name)


def get_data_products(prog_id, inst, calib_level, product_type):
	obs_list = get_proposal_data(prog_id, inst)
	all_products = [Observations.get_product_list(obs) for obs in obs_list]
	file_names = unique(vstack(all_products), keys='productFilename')
	return Observations.filter_products(file_names, calib_level=[calib_level], productSubGroupDescription=product_type)


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
	args = parser.parse_args()

	if args.download:
		download_file(args.download)
		return


if __name__ == '__main__':
	main()
