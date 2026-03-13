import multiprocessing as mp
import pathlib

from daedalus.instruments.utilities import flag_snowballs, run_nsclean

from jwst.pipeline.calwebb_detector1 import Detector1Pipeline

INSTRUMENTS_DIR = pathlib.Path(__file__).parent

SUPPORTED_INSTRUMENTS = ['miri', 'nirspec_ifu', 'nirspec_mos']

POST_JUMP_STEPS = ['ramp_fit', 'gain_scale'] # as of 1.12.2


class Instrument:

    ACTIONS = {
        'download': 'run_download',
        'pipeline': 'run_pipeline',
    }

    def __init__(self, config):
        self.config = config
        self.config.output_dir = self.config.output_dir.joinpath(str(self.config.target.program), self.config.target.name, self.config.product_name)
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        print('Output directory: %s'%str(self.config.output_dir))
        breakpoint()

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
        self.run_stage1_all(opts)
        return True


    def run_download(self, opts):
        print('run_download')
        return True


    @staticmethod
    def create_stage1_detector(output_file, step_opts):
        detector1 = Detector1Pipeline()
        detector1.output_dir = str(output_file.parent)
        detector1.output_file = str(output_file)

        for step, options in step_opts.items():
            for opt, val in options.items():
                setattr(getattr(detector1, step), opt, val)

        return detector1


    def run_pipeline_old(self, context, args):
        self.run_stage1_all(context, args)
        self.run_stage2_all(context, args)
        self.run_stage3_all(context, args)
        print('Pipeline for {pname} complete!'.format(pname=context.config.product_name))


    def run_stage1_sb_flagging(self, ufile, output_dir, context, args):
        stage_args = args.get('stage1', {})
        step_opts = stage_args.get('steps', {})

        detector1 = create_stage1_detector(output_dir.joinpath(ufile.stem), step_opts)
        detector1.save_calibrated_ramp = True

        # First run everything up to and including 'jump' step
        for name in detector1.step_defs.keys():
            if name in POST_JUMP_STEPS:
                getattr(detector1, name).skip = True
        detector1.run(ufile)

        ramp_file = output_dir.joinpath(ufile.name.replace('uncal', 'ramp'))
        if not ramp_file.exists():
            raise Exception('Failed to create ramp file: %s' % str(ramp_file))

        print('Flagging snowballs')
        sb_file = flag_snowballs(ramp_file)

        # Then run everything after 'jump' step with the snowball flagged ramp file
        detector1 = create_stage1_detector(output_dir.joinpath(ufile.stem), step_opts)
        for name in detector1.step_defs.keys():
            if not name in POST_JUMP_STEPS:
                getattr(detector1, name).skip = True
        detector1.run(sb_file)


    def run_stage1_single(self, ufile, opts):
        print('Processing: {}'.format(str(ufile)))

        out_file = self.config.stage2_dir.joinpath(ufile.name.replace('uncal', 'rate'))
        if out_file.exists() and not opts.overwrite:
            return None

        detector1 = Instrument.create_stage1_detector(self.config.stage2_dir.joinpath(ufile.stem), opts.stage1.steps)
        detector1.run(ufile)
        return None

        if opts.stage1.flag_snowballs:
            self.run_stage1_sb_flagging(ufile, output_dir, context, args)
        else:
            step_opts = stage_args.get('steps', {})
            detector1 = create_stage1_detector(output_dir.joinpath(ufile.stem), step_opts)
            detector1.run(ufile)

        if stage_args.get('clean_rates', False):
            print('Running NSClean on %s' % str(out_file))
            run_nsclean(out_file)

        return None


    def run_stage1_all(self, opts):
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
