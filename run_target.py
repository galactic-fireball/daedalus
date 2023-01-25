import pathlib
import sys
import toml

JWST_UTILS_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(JWST_UTILS_DIR))

from instruments.instrument_common import Instrument
from instruments.miri_ifu import MIRI_IFU
from instruments.nirspec_ifu import NIRSpec_IFU
from instruments.nirspec_mos import NIRSpec_MOS

# See jwst-utils/badass for possible options files
BADASS_OPTIONS_FILE = 'ir_options1.py'


INSTRUMENTS = {
	'miri': MIRI_IFU,
	'nirspec_ifu': NIRSpec_IFU,
	'nirspec_mos': NIRSpec_MOS,
}


def main():
	if len(sys.argv) < 2:
		raise Exception('Config file expected')

	config_file = pathlib.Path(sys.argv[1]).resolve()
	if not config_file.exists():
		raise Exception('Config file %s not found' % str(config_file))

	config = toml.load(config_file)
	instrument = Instrument.from_config(config, INSTRUMENTS)

	# TODO: workaround for now, need a better way to handle this
	if len(sys.argv) > 2:
		instrument.badass['line'] = sys.argv[2]

	instrument.run()


if __name__ == '__main__':
	main()
