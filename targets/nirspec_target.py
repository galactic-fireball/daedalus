from astroquery.mast import Observations
import pandas as pd
import pathlib
import sys

JWST_UTILS_DIR = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(JWST_UTILS_DIR))

import mast.mast as mast
import pipeline.nirspec as nirspec
# import badass.badass_nirspec as badass_nirspec
import pipeline.crds_utils as crds_utils

class NirspecProgram:

    download_funcs = {
        'all': 'download_all',
        'uncal': 'download_uncal',
        '1d': 'download_mast_final',
        'asn2': 'download_asn2',
    }

    def __init__(self, attr_dict):
        self.__dict__.update(attr_dict)

        # for attr in ['program_id', 'redshift', 'sci_obs', 'target_name', 'data_dir', 'product_name']:
        for attr in ['program_id', 'redshift', 'sci_obs', 'data_dir', 'product_name']:
            if not hasattr(self, attr) or getattr(self, attr) is None:
                raise Exception('Program input missing expected value: {attr}'.format(attr=attr))

        self.data_dir.mkdir(exist_ok=True, parents=True)
        if not hasattr(self, 'mast_cube_prefix'):
            self.mast_cube_prefix = None


    def download_all(self):
        cached_file = self.data_dir.joinpath('nirspec_{pid}.csv'.format(pid=self.program_id))
        if cached_file.exists():
            df = pd.read_csv(cached_file)
        else:
            df = mast.get_program_data(str(self.program_id), 'NIRSPEC')
            # df.to_csv(cached_file, index=False)

        dl_files = []

        # uncal
        uncal_df = Observations.filter_products(df, calib_level=[1], productSubGroupDescription='UNCAL').to_pandas()
        obs_fmt = 'jw%05d%03d' % (self.program_id, self.sci_obs)
        obs_df = uncal_df[uncal_df.obs_id.str.contains(obs_fmt)]

        # dl_data_dir = self.data_dir.joinpath('input', 'stage1', 'sci')
        dl_data_dir = self.data_dir.joinpath('pipeline_data')
        dl_data_dir.mkdir(exist_ok=True, parents=True)

        for file_name in obs_df.productFilename.values:
            dl_files.append((file_name, dl_data_dir.joinpath(file_name)))

        # asn2
        asn2_df = Observations.filter_products(df, calib_level=[2], productSubGroupDescription='ASN').to_pandas()
        obs_fmt = 'jw%05d%03d' % (self.program_id, self.sci_obs)
        obs_df = asn2_df[(asn2_df.obs_id.str.contains(obs_fmt)) & (asn2_df.productFilename.str.contains('spec'))]

        # dl_data_dir = self.data_dir.joinpath('input', 'stage2', 'sci')
        # dl_data_dir.mkdir(exist_ok=True, parents=True)

        for file_name in obs_df.productFilename.values:
            dl_files.append((file_name, dl_data_dir.joinpath(file_name)))

        # msa
        msa_df = Observations.filter_products(df, productSubGroupDescription='MSA').to_pandas()
        obs_fmt = 'jw%05d%03d' % (self.program_id, self.sci_obs)
        obs_df = msa_df[msa_df.obs_id.str.contains(obs_fmt)]

        # dl_data_dir = self.data_dir.joinpath('input', 'stage2', 'sci')
        # dl_data_dir.mkdir(exist_ok=True, parents=True)

        for file_name in obs_df.productFilename.values:
            dl_files.append((file_name, dl_data_dir.joinpath(file_name)))

        # asn3
        asn3_df = Observations.filter_products(df, calib_level=[3], productSubGroupDescription='ASN').to_pandas()
        obs_fmt = 'jw%05d-o%03d' % (self.program_id, self.sci_obs)
        obs_df = asn3_df[asn3_df.obs_id.str.contains(obs_fmt)]

        # dl_data_dir = self.data_dir.joinpath('input', 'stage3', 'sci')
        # dl_data_dir.mkdir(exist_ok=True, parents=True)

        for file_name in obs_df.productFilename.values:
            dl_files.append((file_name, dl_data_dir.joinpath(file_name)))

        for file_name, dest in dl_files:
            mast.download_file(file_name, dest=dest)


    def download_uncal(self):
        cached_file = self.data_dir.joinpath('nirspec_{pid}_uncal.csv'.format(pid=self.program_id))
        if cached_file.exists():
            df = pd.read_csv(cached_file)
        else:
            df = mast.get_data_products(str(self.program_id), 'NIRSPEC', 1, 'UNCAL').to_pandas()
            df.to_csv(cached_file, index=False)

        obs_fmt = 'jw%05d%03d' % (self.program_id, self.sci_obs)
        obs_df = df[df.obs_id.str.contains(obs_fmt)]

        dl_data_dir = self.data_dir.joinpath('input', 'stage1', 'sci')
        dl_data_dir.mkdir(exist_ok=True, parents=True)

        for file_name in obs_df.productFilename.values:
            mast.download_file(file_name, dest=dl_data_dir.joinpath(file_name))


    def download_asn2(self):
        cached_file = self.data_dir.joinpath('nirspec_{pid}_spec2_asn.csv'.format(pid=self.program_id))
        if cached_file.exists():
            df = pd.read_csv(cached_file)
        else:
            df = mast.get_data_products(str(self.program_id), 'NIRSPEC', 2, 'ASN').to_pandas()
            df.to_csv(cached_file, index=False)

        obs_fmt = 'jw%05d%03d' % (self.program_id, self.sci_obs)
        obs_df = df[(df.obs_id.str.contains(obs_fmt)) & (df.productFilename.str.contains('spec'))]

        dl_data_dir = self.data_dir.joinpath('input', 'stage2', 'sci')
        dl_data_dir.mkdir(exist_ok=True, parents=True)

        for file_name in obs_df.productFilename.values:
            mast.download_file(file_name, dest=dl_data_dir.joinpath(file_name))


    def download_mast_final(self):
        cached_file = self.data_dir.joinpath('nirspec_{pid}_x1d.csv'.format(pid=self.program_id))
        if cached_file.exists():
            df = pd.read_csv(cached_file)
        else:
            df = mast.get_data_products(str(self.program_id), 'NIRSPEC', 3, 'X1D').to_pandas()
            df.to_csv(cached_file, index=False)

        # obs_id : jw02736-o008_s06355_nirspec_f170lp-g235m
        sci_fmt = 'jw%05d-o%03d_s%05d' % (self.program_id, self.sci_obs, self.target_name)
        sci_df = df[df.obs_id.str.contains(sci_fmt)]

        dl_data_dir = self.data_dir.joinpath('output', 'stage3', 'sci')
        dl_data_dir.mkdir(exist_ok=True, parents=True)

        for file_name in sci_df.productFilename.values:
            mast.download_file(file_name, dest=dl_data_dir.joinpath(file_name))


    def run_pipeline(self):
        crds_utils.cache_crds('nirspec')

        # STAGE 1
        # input_dir = self.data_dir.joinpath('input', 'stage1', 'sci')
        # output_dir = self.data_dir.joinpath('output', 'stage1', 'sci')
        # output_dir.mkdir(exist_ok=True, parents=True)
        input_dir = self.data_dir.joinpath('pipeline_data')
        output_dir = input_dir
        nirspec.run_stage1_all(input_dir, output_dir)

        # STAGE 2
        # input_dir = self.data_dir.joinpath('input', 'stage2', 'sci')
        # # Move the stage 1 output to stage 2 input
        # input_dir.mkdir(exist_ok=True, parents=True)
        # for f in output_dir.glob('*'):
        #   f.rename(input_dir.joinpath(f.name))

        # output_dir = self.data_dir.joinpath('output', 'stage2', 'sci')
        # output_dir.mkdir(exist_ok=True, parents=True)
        nirspec.run_stage2_all(input_dir, output_dir)

        # STAGE 3
        # input_dir = self.data_dir.joinpath('input', 'stage3', 'sci')
        # Move the stage 2 output to stage 3 input
        # input_dir.mkdir(exist_ok=True, parents=True)
        # for f in output_dir.glob('*'):
        #     f.rename(input_dir.joinpath(f.name))

        # output_dir = self.data_dir.joinpath('output', 'stage3', 'sci')
        # output_dir.mkdir(exist_ok=True, parents=True)
        nirspec.run_stage3_all(input_dir, output_dir)

        print('Pipeline for {pname} complete!'.format(pname=self.product_name))
        self.pipeline_out_dir = output_dir


    def init_for_badass(self):
        if not hasattr(self, 'pipeline_out_dir'):
            raise Exception('Set pipeline output directory!')

        self.badass_dir = self.data_dir.joinpath('badass')
        self.badass_dir.mkdir(exist_ok=True, parents=True)

        cubes = list(self.pipeline_out_dir.glob('*-shortmediumlong_s3d.fits'))
        print('Found {n} cubes'.format(n=len(cubes)))
        for cube in cubes:
            cube.rename(self.badass_dir.joinpath(cube.name))

        # grab the list this way since we might not have actually moved anything
        cubes = self.badass_dir.glob('*_s3d.fits')
        for cube in cubes:
            spaxel_dir = cube.with_suffix('')
            if not spaxel_dir.exists():
                channel = int(re.findall(r'_ch\d+-shortmediumlong_s3d', cube.name)[0].split('-')[0][3:])
                specres = badass_miri.get_channel_resolution(channel, 'short')
                badass_miri.prepare_cube(cube, specres, z=self.redshift)







