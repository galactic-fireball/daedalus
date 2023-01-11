import argparse
import os
import pathlib
import sys

import ssl
# Needed to work around ssl certificate verification during crds downloads
ssl._create_default_https_context = ssl._create_unverified_context

USE_CRDS_OPS = True
# Needs to be set before crds/jwst imports
if USE_CRDS_OPS:
    os.environ['CRDS_PATH'] = str(pathlib.Path(__file__).resolve().parent.joinpath('crds_cache', 'ops'))
    os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu'
else:
    os.environ['CRDS_PATH'] = str(pathlib.Path(__file__).resolve().parent.joinpath('crds_cache', 'pub'))
    os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds-pub.stsci.edu'

import crds
from crds.client.api import cache_references
from crds.core.exceptions import IrrelevantReferenceTypeError

supported_instruments = ['miri', 'nirspec']

# Adapted from spacetelescope/jwst/jwst/datamodels/tests/test_schema_against_crds.py
def cache_crds(instrument):
    context = crds.get_default_context('jwst')
    crds.api.dump_mappings(context)
    pmap = crds.get_cached_mapping(context)
    imap = pmap.get_imap(instrument)

    reftypes = imap.get_filekinds()
    _ = [reftypes.remove(name) for name in reftypes[::-1] if name.startswith('pars-')]

    for reftype in reftypes:
        try:
            r = imap.get_rmap(reftype)
            for f in r.reference_names():
                if 'fits' in f or 'asdf' in f:
                    refs = cache_references(context, {reftype: f})
        except IrrelevantReferenceTypeError as e:
            pass


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--instrument', help='Instrument to cache crds files for', type=str, choices=supported_instruments, required=True)
    args = parser.parse_args()

    cache_crds(args.instrument)


if __name__ == '__main__':
    main()

