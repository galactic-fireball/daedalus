import os
import pathlib
from prodict import Prodict
import sys
import toml

import ssl
# Needed to work around ssl certificate verification during crds downloads
ssl._create_default_https_context = ssl._create_unverified_context

JWST_UTILS_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(JWST_UTILS_DIR))

DEFAULT_CONFIG = JWST_UTILS_DIR.joinpath('default.toml')
PROGRAMS_DIR = pathlib.Path(__file__).parent.joinpath('programs')

def configure_crds(cache=None, use_ops=True):
    if use_ops:
        os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu'
        if cache:
            os.environ['CRDS_PATH'] = str(pathlib.Path(cache).joinpath('ops'))
            os.environ['CRDS_CONFIG_URI'] = str(pathlib.Path(cache).joinpath('ops', 'config', 'jwst'))
    else:
        os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds-pub.stsci.edu'
        if cache:
            os.environ['CRDS_PATH'] = str(pathlib.Path(cache).joinpath('pub'))
            os.environ['CRDS_CONFIG_URI'] = str(pathlib.Path(cache).joinpath('pub', 'config', 'jwst'))


def get_config():
    if DEFAULT_CONFIG.exists():
        config = toml.load(DEFAULT_CONFIG)
    else:
        config = {}

    if len(sys.argv) < 2:
        return config

    config_file = pathlib.Path(sys.argv[1]).resolve()
    if not config_file.exists():
        raise Exception('Config file %s not found' % str(config_file))

    config.update(toml.load(config_file))
    return Prodict.from_dict(config)


config = get_config()
if 'CRDS_CACHE' in config or 'USE_CRDS_OPS' in config:
    configure_crds(cache=config.get('CRDS_CACHE', None), use_ops=config.get('USE_CRDS_OPS', True))


import instruments.instrument_common as ic

from instruments.miri_ifu import MIRI_IFU
from instruments.nirspec_ifu import NIRSpec_IFU
from instruments.nirspec_mos import NIRSpec_MOS

INSTRUMENTS = {
    'miri': MIRI_IFU,
    'nirspec_ifu': NIRSpec_IFU,
    'nirspec_mos': NIRSpec_MOS,
}

import actions
ACTIONS = actions.get_actions()


class RunContext():
    def __init__(self, config):
        self.config = config

        self.target = None # TODO: set default variables
        if 'target' in self.config:
            self.target = Prodict.from_dict(self.config.target)

        self.instrument = None # TODO: set default variables
        if 'instrument' in self.config:
            inst_name = self.config.instrument.name
            if not inst_name in INSTRUMENTS:
                raise Exception('Unregistered instrument: %s' % inst_name)
            self.instrument = INSTRUMENTS[inst_name].from_dict(self.config.instrument)

        self.actions = []
        for action in self.config.action:
            aname = action['name']
            if not aname in ACTIONS:
                raise Exception('Unregistered action: %s' % action)
            self.actions.append({'args':action,'fnc':ACTIONS[aname]})

        self.data_dir = PROGRAMS_DIR.joinpath(str(self.target.program_id), 'data_sets', self.target.name, self.instrument.name, self.config.product_name)
        self.data_dir.mkdir(parents=True, exist_ok=True)

        # subdirectories: mkdir these in specific actions as needed
        self.pipeline_dir = self.data_dir.joinpath('pipeline')
        self.badass_dir = self.data_dir.joinpath('badass')

        self.instrument.post_setup(self)


    def run(self):
        results = []
        for action in self.actions:
            results.append(action['fnc'](self, action['args']))
        return results


def main():
    RunContext(config).run()

if __name__ == '__main__':
    main()
