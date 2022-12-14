import argparse
import pathlib
import sys
import toml

from miri_program import MiriProgram
from nirspec_program import NirspecProgram

CUR_DIR = pathlib.Path(__file__).resolve().parent

# See jwst-utils/badass for possible options files
BADASS_OPTIONS_FILE = 'ir_options1.py'


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--config', help='Program/target config file', type=str, required=True)
	parser.add_argument('--download', help='Download uncal or cube', type=str, default=None)
	parser.add_argument('--pipeline', help='Run through pipeline', action='store_true', default=False)
	parser.add_argument('--line', help='Line region for BADASS', type=str, default=None)
	args = parser.parse_args()

	config_file = pathlib.Path(args.config).resolve()
	if not config_file.exists():
		raise Exception('Config file %s not found' % str(config_file))

	config = toml.load(config_file)

	if 'data_dir' in config:
		config['data_dir'] = CUR_DIR.joinpath(config['data_dir'])

	instrument = config['instrument']
	if instrument == 'miri':
		program = MiriProgram(config)
	elif instrument == 'nirspec':
		program = NirspecProgram(config)
	else:
		raise Exception('Unknown instrument: %s' % instrument)

	if args.download:
		getattr(program, program.download_funcs[args.download])()

	if args.pipeline:
		program.run_pipeline()

	if args.line:
		if not args.pipeline:
			program.pipeline_out_dir = program.data_dir.joinpath('output', 'stage3', 'sci')
		program.init_for_badass()
		program.run_badass_line_region(args.line, BADASS_OPTIONS_FILE)
		program.run_cube_reconstruct(args.line)


if __name__ == '__main__':
	main()
