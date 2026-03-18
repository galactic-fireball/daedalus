from astropy.io import fits
import copy
import json
import multiprocessing as mp
import pandas as pd
import pathlib
from typing import NamedTuple

import daedalus.mast.mast as mast
from daedalus.instruments.instrument_common import Instrument

from jwst.associations import asn_from_list as afl
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase
from jwst.associations import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline
from jwst.pipeline.calwebb_spec3 import Spec3Pipeline


def file_to_exp(file_name):
    if isinstance(file_name, str):
        return file_name.rsplit('_',1)[0]
    if isinstance(file_name, pathlib.Path):
        return file_name.stem.rsplit('_',1)[0]
    return None


def exp_to_file(exp, ftype, ext='.fits'):
    return exp + '_' + ftype + ext


def file_to_file(file_name, ftype, ext1='.fits', ext2='.fits'):
    if isinstance(file_name, str):
        cur_ftype = file_name.rsplit('_',1)[1].replace(ext1,'')
        return file_name.replace(cur_ftype, ftype).replace(ext1, ext2)
    if isinstance(file_name, pathlib.Path):
        cur_ftype = file_name.stem.rsplit('_',1)[1]
        return file_name.parent.joinpath(file_name.name.replace(cur_ftype, ftype).replace(ext1, ext2))
    return None


class ChBand(NamedTuple):
    ch: str
    band: str


class MIRI_IFU(Instrument):

    EXP_LISTS_FILE = 'exp_lists.json'


    def __init__(self, config):
        super().__init__(config)
        self.sci_obs = self.config.target.instrument.sci_obs
        self.bkgd_obs = self.config.target.instrument.bkgd_obs

        # store list of exposure names for science and background
        self.sci_exps = []
        self.bkgd_exps = []

        # this will be initialized later so we only have to find
        # the ch and band of each file once
        self.file_ch_band = {}


    def fill_exp_lists(self):
        if len(self.sci_exps) > 0:
            return

        exp_file = self.config.uncal_dir.joinpath(self.EXP_LISTS_FILE)
        if not exp_file.exists():
            raise Exception('exposure file not found: %s'%str(exp_file))

        with open(exp_file, 'r') as f:
            exp_lists = json.load(f)

        self.sci_exps = exp_lists['sci']
        self.bkgd_exps = exp_lists['bkgd']


    def run_download(self, opts):
        if not opts.mast_token is None:
            mast.login(opts.mast_token)

        program_id = self.config.target.program
        df = mast.get_data_products(str(program_id), 'MIRI/IFU', 1, 'UNCAL').to_pandas()

        # Get only IFU data
        df = df[(df.obs_id.str.contains('mirifulong')) | (df.obs_id.str.contains('mirifushort'))]

        # Download science files
        sci_fmt = 'jw%05d%03d' % (program_id, self.sci_obs)
        sci_df = df[df.obs_id.str.contains(sci_fmt)]

        for file_name in sci_df.productFilename.values:
            self.sci_exps.append(file_to_exp(file_name))
            dest = self.config.uncal_dir.joinpath(file_name)
            if dest.exists():
                continue
            print('downloading %s' % file_name)
            mast.download_file(file_name, dest=dest)

        # Download background files
        bkgd_fmt = 'jw%05d%03d' % (program_id, self.bkgd_obs)
        bkgd_df = df[df.obs_id.str.contains(bkgd_fmt)]

        for file_name in bkgd_df.productFilename.values:
            self.bkgd_exps.append(file_to_exp(file_name))
            dest = self.config.uncal_dir.joinpath(file_name)
            if dest.exists():
                continue
            print('downloading %s' % file_name)
            mast.download_file(file_name, dest=dest)

        exp_file = self.config.uncal_dir.joinpath(self.EXP_LISTS_FILE)
        with open(exp_file, 'w') as f:
            json.dump({'sci':self.sci_exps,'bkgd':self.bkgd_exps}, f, indent=2)


    # TODO: fix and add as download option
    def download_mast_final():
        cached_file = self.data_dir.joinpath('miri_{pid}_s3d.csv'.format(pid=self.program_id))
        if cached_file.exists():
            df = pd.read_csv(cached_file)
        else:
            df = mast.get_data_products(str(self.program_id), 'MIRI', 3, 'S3D').to_pandas()
            df.to_csv(cached_file, index=False)

        sci_fmt = 'jw%05d-o%03d' % (self.program_id, self.sci_obs)
        sci_df = df[df.obs_id.str.contains(sci_fmt)]

        self.mast_dir = self.data_dir.joinpath('mast')
        self.mast_dir.mkdir(parents=True, exist_ok=True)

        for file_name in sci_df.productFilename.values:
            if not self.mast_cube_prefix:
                self.mast_cube_prefix = file_name.name.split('-shortmediumlong_s3d.fits')[0].split('_ch')[0]
            mast.download_file(file_name, dest=self.mast_dir.joinpath(file_name))


    def init_ch_band_dict(self, files):
        for file in files:
            exp_name = file_to_exp(file)
            with fits.open(file) as hdu:
                hdr = hdu[0].header
                self.file_ch_band[exp_name] = ChBand(ch=hdr['CHANNEL'], band=hdr['BAND'])


    # The following functions for stage 2 and 3 were modified from David Law and Kirsten Larson's jupyter notebook:
    #   https://github.com/spacetelescope/jwst-pipeline-notebooks/blob/main/notebooks/MIRI/MRS/JWPipeNB-MIRI-MRS.ipynb
    def create_stage2_asn(self, rfile, bkgd_files, all_files):
        asn = afl.asn_from_list([str(rfile),], rule=DMSLevel2bBase, product_name=self.config.product_name+'_level2')

        target_ch_band = self.file_ch_band.get(file_to_exp(rfile), None)
        if target_ch_band is None:
            raise Exception('File %s ch and band not found'%str(rfile))

        for file in bkgd_files:
            ch_band = self.file_ch_band.get(file_to_exp(file), None)
            if ch_band is None:
                raise Exception('File %s ch and band not found'%str(file))

            if (ch_band.ch == target_ch_band.ch) and (ch_band.band == target_ch_band.band):
                asn['products'][0]['members'].append({'expname':str(file), 'exptype':'background'})

        for file in all_files:
            ch_band = self.file_ch_band.get(file_to_exp(file), None)
            if ch_band is None:
                raise Exception('File %s ch and band not found'%str(file))

            if (ch_band.ch == target_ch_band.ch):
                asn['products'][0]['members'].append({'expname':str(file), 'exptype':'selfcal'})

        asn_file = file_to_file(rfile, 'asn', ext2='.json')
        _,asn_serialized = asn.dump()
        with open(asn_file, 'w') as f:
            f.write(asn_serialized)
        return asn_file


    def run_stage2_single(self, infile, opts, bkgd_files, all_files):
        asn_file = str(self.create_stage2_asn(infile, bkgd_files, all_files))

        build_args = {
            'output_dir': str(self.config.stage3_dir),
            'save_results': True,
            'steps': opts.stage3.steps,
        }

        config, _ = Spec2Pipeline.build_config(asn_file, **build_args)
        print('Running Spec 2 Pipeline with config: \n'+json.dumps(config.dict(), indent=4))
        Spec2Pipeline.call(asn_file, **build_args)


    def run_stage2(self, opts):
        sci_rate_files = [self.config.stage2_dir.joinpath(exp_to_file(exp,'rate')) for exp in self.sci_exps]
        bkgd_rate_files = [self.config.stage2_dir.joinpath(exp_to_file(exp,'rate')) for exp in self.bkgd_exps]

        all_rate_files = sci_rate_files + bkgd_rate_files
        self.init_ch_band_dict(all_rate_files)

        if not opts.multiprocess:
            for rfile in sci_rate_files:
                self.run_stage2_single(rfile, opts, bkgd_rate_files, all_rate_files)
            for rfile in bkgd_rate_files:
                self.run_stage2_single(rfile, opts, [], all_rate_files)
            return

        proc_args = [(rfile, opts, bkgd_rate_files, all_rate_files) for rfile in sci_rate_files]
        pool = mp.Pool(processes=opts.nprocesses, maxtasksperchild=1)
        pool.starmap(self.run_stage2_single, proc_args, chunksize=1)
        pool.close()
        pool.join()

        proc_args = [(rfile, opts, [], all_rate_files) for rfile in bkgd_rate_files]
        pool = mp.Pool(processes=opts.nprocesses, maxtasksperchild=1)
        pool.starmap(self.run_stage2_single, proc_args, chunksize=1)
        pool.close()
        pool.join()


    def run_stage3(self, opts):
        sci_files = [self.config.stage3_dir.joinpath(exp_to_file(exp,'cal')) for exp in self.sci_exps]
        bkgd_files = [self.config.stage3_dir.joinpath(exp_to_file(exp,'x1d')) for exp in self.bkgd_exps]

        asn = afl.asn_from_list(sci_files, rule=DMS_Level3_Base, product_name=self.config.product_name+'_level3')

        for file in bkgd_files:
            asn['products'][0]['members'].append({'expname':str(file), 'exptype':'background'})

        asn_file = self.config.output_dir.joinpath('l3asn.json')
        _,asn_serialized = asn.dump()
        with open(asn_file, 'w') as f:
            f.write(asn_serialized)
        asn_file = str(asn_file)

        build_args = {
            'output_dir': str(self.config.output_dir),
            'save_results': True,
            'steps': opts.stage3.steps,
        }

        config, _ = Spec3Pipeline.build_config(asn_file, **build_args)
        print('Running Spec 3 Pipeline with config: \n'+json.dumps(config.dict(), indent=4))
        Spec3Pipeline.call(asn_file, **build_args)

