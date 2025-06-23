import multiprocessing as mp
import pathlib
from prodict import Prodict

from instruments.utilities import flag_snowballs, run_nsclean

from jwst.pipeline.calwebb_detector1 import Detector1Pipeline

INSTRUMENTS_DIR = pathlib.Path(__file__).parent

POST_JUMP_STEPS = ['ramp_fit', 'gain_scale'] # as of 1.12.2


def create_stage1_detector(output_file, step_opts):
    detector1 = Detector1Pipeline()
    detector1.output_dir = str(output_file.parent)
    detector1.output_file = str(output_file)

    for step, options in step_opts.items():
        for opt, val in options.items():
            setattr(getattr(detector1, step), opt, val)

    return detector1


class Instrument(Prodict):

    def post_setup(self, context):
        pass


    def run_pipeline(self, context, args):
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


    def run_stage1_single(self, ufile, output_dir, context, args):
        print('Processing: {}'.format(str(ufile)))

        stage_args = args.get('stage1', {})
        overwrite = stage_args.get('overwrite', False)

        out_file = output_dir.joinpath(ufile.name.replace('uncal', 'rate'))
        if out_file.exists() and not overwrite:
            return None

        if stage_args.get('flag_snowballs', False):
            self.run_stage1_sb_flagging(ufile, output_dir, context, args)
        else:
            step_opts = stage_args.get('steps', {})
            detector1 = create_stage1_detector(output_dir.joinpath(ufile.stem), step_opts)
            detector1.run(ufile)

        if stage_args.get('clean_rates', False):
            print('Running NSClean on %s' % str(out_file))
            run_nsclean(out_file)

        return None


    def run_stage1_all(self, context, args, indir=None, outdir=None):
        stage1_opts = args.get('stage1', {})
        input_dir = stage1_opts.get('input_dir', indir if indir else context.pipeline_dir)
        output_dir = stage1_opts.get('output_dir', outdir if outdir else context.pipeline_dir)
        multiprocess = args.get('multiprocess', False)
        nprocesses = args.get('nprocesses', 1)

        uncal_files = input_dir.glob('*_uncal.fits')

        if not multiprocess:
            for ufile in uncal_files:
                self.run_stage1_single(ufile, output_dir, context, args)
            return

        proc_args = [(ufile, output_dir, context, args) for ufile in uncal_files]
        pool = mp.Pool(processes=nprocesses, maxtasksperchild=1)
        pool.starmap(self.run_stage1_single, proc_args, chunksize=1)
        pool.close()
        pool.join()
