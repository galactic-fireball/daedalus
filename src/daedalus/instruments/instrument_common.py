from abc import ABC, abstractmethod
import json
import multiprocessing as mp
import pathlib

from daedalus.instruments.utilities import flag_snowballs, run_nsclean

from jwst.pipeline.calwebb_detector1 import Detector1Pipeline

INSTRUMENTS_DIR = pathlib.Path(__file__).parent

SUPPORTED_INSTRUMENTS = ['miri', 'nirspec_ifu', 'nirspec_mos']

POST_JUMP_STEPS = ['ramp_fit', 'gain_scale'] # as of 1.12.2


class Instrument(ABC):

    ACTIONS = {
        'download': 'run_download',
        'pipeline': 'run_pipeline',
    }

    def __init__(self, config):
        self.config = config

        # output_dir = <user_out>/<program_id>/<target_name>/<instrument_name>/<product_name>/
        self.config.output_dir = self.config.output_dir.joinpath(
            str(self.config.target.program),
            self.config.target.name,
            self.config.target.instrument.name,
            self.config.product_name
        )
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        print('Output directory: %s'%str(self.config.output_dir))

        if self.config.uncal_dir is None:
            self.config.uncal_dir = self.config.output_dir.joinpath('uncal')
        self.config.uncal_dir.mkdir(parents=True, exist_ok=True)

        if self.config.stage2_dir is None:
            self.config.stage2_dir = self.config.output_dir.joinpath('stage2')
        self.config.stage2_dir.mkdir(parents=True, exist_ok=True)

        if self.config.stage3_dir is None:
            self.config.stage3_dir = self.config.output_dir.joinpath('stage3')
        self.config.stage3_dir.mkdir(parents=True, exist_ok=True)


    def run_actions(self):
        res = []
        for action in self.config.actions:
            if not action.name in Instrument.ACTIONS:
                print('Specified action [%s] unsupported'%action.name)
            print('Running action [%s]'%action.name)
            res.append(getattr(self, Instrument.ACTIONS[action.name])(action))
        return res


    def run_pipeline(self, opts):
        print('run_pipeline')

        if not opts.stage1.skip:
            self.run_stage1(opts)

        if not opts.stage2.skip:
            self.run_stage2(opts)

        if not opts.stage3.skip:
            self.run_stage3(opts)

        return True


    @abstractmethod
    def run_download(self, opts):
        pass


    def run_stage1_single(self, ufile, opts):
        print('Processing: {}'.format(str(ufile)))

        out_file = self.config.stage2_dir.joinpath(ufile.name.replace('uncal', 'rate'))
        if out_file.exists() and not opts.overwrite:
            return

        build_args = {
            'output_dir': str(out_file.parent),
            'output_file': str(out_file),
            'steps': opts.stage1.steps,
        }

        config, _ = Detector1Pipeline.build_config(ufile, **build_args)
        print('Running Detector 1 Pipeline with config: \n'+json.dumps(config.dict(), indent=4))
        Detector1Pipeline.call(ufile, **build_args)


    def run_stage1(self, opts):

        # TODO: move all option parsing up to Instrument.__init__
        if opts.stage1.flag_snowballs:
            if not 'jump' in opts.stage1.steps:
                opts.stage1.steps['jump'] = {}
            if self.config.target.instrument in ['nirspec_ifu', 'nirspec_mos',]:
                # As of JWST 1.20.2 / CRDS 13.1.10, this is on by default for NIR observations, but make it explicit
                opts.stage1.steps['jump']['expand_large_events'] = True
            elif self.config.target.instrument in ['miri']:
                opts.stage1.steps['jump']['find_showers'] = True

        # As of JWST 1.20.2 / CRDS 13.1.10, this step is skipped by default, turn on here
        if opts.stage1.clean_noise:
            if not 'clean_flicker_noise' in opts.stage1.steps:
                opts.stage1.steps['clean_flicker_noise'] = {}
            opts.stage1.steps['clean_flicker_noise']['skip'] = False

        uncal_files = self.config.uncal_dir.glob('*_uncal.fits')
        if not opts.multiprocess:
            for ufile in uncal_files:
                self.run_stage1_single(ufile, opts)
            return

        proc_args = [(ufile, opts) for ufile in uncal_files]
        pool = mp.Pool(processes=opts.nprocesses, maxtasksperchild=1)
        pool.starmap(self.run_stage1_single, proc_args, chunksize=1)
        pool.close()
        pool.join()


    @abstractmethod
    def run_stage2(self, opts):
        pass


    @abstractmethod
    def run_stage3(self, opts):
        pass


