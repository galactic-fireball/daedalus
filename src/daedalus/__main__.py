import os
import pathlib
import requests
import sys
import toml

SC_URL = 'https://jwst-crds.stsci.edu/unchecked_get/config/jwst/server_config'
PMAP_URL = 'https://jwst-crds.stsci.edu/unchecked_get/mappings/jwst/%s'
REPO_DIR = pathlib.Path(__file__).resolve().parent.parent.parent


instrument_map = {
    'nirspec_ifu': 'nirspec',
    'nirspec_mos': 'nirspec',
}


def run(config):
    instrument = config.target.instrument.create_instrument(config)
    res = instrument.run_actions()
    return all(res)


def main():

    import ssl
    # Needed to work around ssl certificate verification during crds downloads
    ssl._create_default_https_context = ssl._create_unverified_context

    if len(sys.argv) < 2:
        print('Config file required, see usage')
        return


    config_file = pathlib.Path(sys.argv[1]).resolve()
    if not config_file.exists():
        print('Config file %s not found'%str(config_file))
        return


    # intercept config and configure crds before anything else is imported
    toml_config = toml.load(config_file)
    crds_cache = toml_config.get('crds_cache', None)
    os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu'

    if crds_cache:
        crds_cache = pathlib.Path(crds_cache)
        if not crds_cache.is_absolute():
            crds_cache = REPO_DIR.joinpath(crds_cache)

        print('Setting CRDS cache: %s'%str(crds_cache))
        os.environ['CRDS_PATH'] = str(crds_cache.joinpath('ops'))
        uri = crds_cache.joinpath('ops', 'config', 'jwst')
        uri.mkdir(parents=True, exist_ok=True)
        os.environ['CRDS_CONFIG_URI'] = str(uri)

        # if the cache is new, we might need the server_config file
        sc_file = crds_cache.joinpath('ops', 'config', 'jwst', 'server_config')
        if not sc_file.exists():
            res = requests.get(SC_URL)
            if not res.ok:
                print('Unable to download server_config. Exiting...')
                return

            sc_file.parent.mkdir(parents=True, exist_ok=True)
            with open(sc_file, 'wb') as scf:
                scf.write(res.content)

        import crds
        pmap = crds.get_default_context('jwst')
        ctx_map_file = crds_cache.joinpath('ops', 'mappings', 'jwst', pmap)
        if not ctx_map_file.exists():
            url = PMAP_URL%pmap
            res = requests.get(url)
            if not res.ok:
                print('Unable to download pmap file: %s. Exiting...'%pmap)
                return

            ctx_map_file.parent.mkdir(parents=True, exist_ok=True)
            with open(ctx_map_file, 'wb') as cmf:
                cmf.write(res.content)

        # FileNotFoundError: [Errno 2] No such file or directory: '/projects/ssatyapa/spectra/jwst/crds_cache/ops/mappings/jwst/jwst_miri_0517.imap'

        # make sure we have all the instrument references mapped
        # instrument = toml_config.get('target', {}).get('instrument', {}).get('name', None)
        # if instrument is None:
        #     print('Need instrument! Exiting...')
        #     return

        # instrument = instrument_map.get(instrument, instrument)
        # mappings = crds.api.get_mapping_names(pmap)
        # imaps = [m for m in mappings if '.imap' in m]
        # target_maps = [m for m in imaps if instrument in m]

        # for mfile in mappings:
        #     local_file = crds_cache.joinpath('ops', 'mappings', 'jwst', mfile)
        #     if local_file.exists():
        #         continue

        #     url = PMAP_URL%mfile
        #     res = requests.get(url)
        #     if not res.ok:
        #         print('Unable to download imap file: %s. Exiting...'%mfile)

        #     with open(local_file, 'wb') as lmf:
        #         lmf.write(res.content)


    from daedalus.models import Config
    config = Config(**toml_config)
    success = run(config)


if __name__ == '__main__':
    main()
