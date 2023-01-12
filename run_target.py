import pathlib
import sys
import toml

JWST_UTILS_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(JWST_UTILS_DIR))

from targets.target_common import Target
from targets.miri_target import MiriTarget
from targets.nirspec_target import NirspecTarget

# See jwst-utils/badass for possible options files
BADASS_OPTIONS_FILE = 'ir_options1.py'


INSTRUMENTS = {
	'miri': MiriTarget,
	'nirspec': NirspecTarget,
}


def main():
	if len(sys.argv) < 2:
		raise Exception('Config file expected')

	config_file = pathlib.Path(sys.argv[1]).resolve()
	if not config_file.exists():
		raise Exception('Config file %s not found' % str(config_file))

	config = toml.load(config_file)
	Target.from_config(config, INSTRUMENTS).run()











	# if args.download:
	# 	getattr(program, program.download_funcs[args.download])()

	# if args.pipeline:
	# 	program.run_pipeline()

	# if args.cube:
	# 	if not args.pipeline:
	# 		program.pipeline_out_dir = program.data_dir.joinpath('output', 'stage3', 'sci')
	# 	program.init_for_badass()

	# if args.line:
	# 	if not args.cube:
	# 		program.badass_dir = program.data_dir.joinpath('badass')
	# 	if args.extract:
	# 		program.run_badass_line_region_extracted(args.line, BADASS_OPTIONS_FILE)
	# 	else:
	# 		program.run_badass_line_region(args.line, BADASS_OPTIONS_FILE)
	# 		program.run_cube_reconstruct(args.line)
	# elif args.extract:
	# 	raise Exception('Extract 1D spectra not yet implemented!')

	# if args.ratio:
	# 	program.create_all_ratio_maps()


if __name__ == '__main__':
	main()
