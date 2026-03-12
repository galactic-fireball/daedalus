import pathlib
import sys
import toml

from daedalus.models import Config

def main():
	if len(sys.argv) < 2:
		print('Config file required, see usage')
		return

	config_file = pathlib.Path(sys.argv[1]).resolve()
	if not config_file.exists():
		print('Config file %s not found'%str(config_file))
		return

	config = Config(**toml.load(config_file))


if __name__ == '__main__':
	main()
