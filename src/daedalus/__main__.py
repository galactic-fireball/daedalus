import os
import pathlib
import requests
import sys
import toml

SC_URL = 'https://jwst-crds.stsci.edu/unchecked_get/config/jwst/server_config'
REPO_DIR = pathlib.Path(__file__).resolve().parent.parent.parent

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

            with open(sc_file, 'wb') as scf:
                scf.write(res.content)


    from daedalus.models import Config
    config = Config(**toml_config)
    success = run(config)


if __name__ == '__main__':
    main()
