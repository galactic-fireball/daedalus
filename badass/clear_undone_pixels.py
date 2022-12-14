import argparse
import pathlib
import shutil


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--spaxels', help='Spaxel dir', type=str, required=True)
	parser.add_argument('--out', help='BADASS output label', type=str, default=None)
	parser.add_argument('--all', help='All output labels', action='store_true', default=False)
	args = parser.parse_args()

	spaxel_dir = pathlib.Path(args.spaxels).resolve()
	if not spaxel_dir.exists():
		raise Exception('Could not find spaxel directory: %s' % str(spaxel_dir))

	if args.out:
		bout = args.out
	elif args.all:
		bout = '*'
	else:
		raise Exception('Not output label provided')

	out_dirs = [d for d in spaxel_dir.glob('spaxel_*/%s' % bout) if d.is_dir()]
	for odir in out_dirs:
		if (not odir.joinpath('log').exists()) or (not odir.joinpath('log', 'par_table.fits').exists()):
			shutil.rmtree(str(odir))


if __name__ == '__main__':
	main()
