import argparse
import pathlib
import sys

from miri_program import MiriProgram

DATA_SET = 'pipeline_1.8.2'
# DATA_SET = 'pipeline_1.6.2'
# DATA_SET = 'mast_direct'

CUR_DIR = pathlib.Path(__file__).resolve().parent

PROGRAMS = {
	'1328': {
		'program_id': 1328,
		'redshift': 0.01627, # NGC 7469
		'mast_cube_prefix': 'jw01328-o015_t014_miri',
		'sci_obs': 15, # For NGC 7469, see https://www.stsci.edu/jwst/phase2-public/1328.pdf
		'bkgd_obs': 16,
		'data_dir': CUR_DIR.joinpath('1328', 'data_sets', DATA_SET),
		'product_name': 'ngc7469_%s'%DATA_SET
	},
	'2732': {
		'program_id': 2732,
		'redshift': 0.022823, # NGC 7319
		'sci_obs': 4, # For NGC 7319, see see https://www.stsci.edu/jwst/phase2-public/2732.pdf
		'bkgd_obs': 5,
		'data_dir': CUR_DIR.joinpath('2732', 'data_sets', DATA_SET),
		'product_name': 'ngc7319_%s'%DATA_SET
	}
}

# See jwst-utils/badass for possible options files
BADASS_OPTIONS_FILE = 'ir_options1.py'


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--program', help='Target JWST program', type=str, choices=PROGRAMS.keys(), required=True)
	parser.add_argument('--download', help='Download uncal or cube', type=str, choices=MiriProgram.download_funcs.keys(), default=None)
	parser.add_argument('--pipeline', help='Run through pipeline', action='store_true', default=False)
	parser.add_argument('--line', help='Line region for BADASS', type=str, default=None)
	args = parser.parse_args()

	program = MiriProgram(PROGRAMS[args.program])

	if args.download:
		getattr(program, MiriProgram.download_funcs[args.download])()

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
