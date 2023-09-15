# jwst-utils

## Actions

### Download

#### Options

* `mast_token`: MAST token needed to login to MAST to access proprietary data.
* `output_dir`: Directory to download JWST uncalibrated data. Defaults to pipeline directory.

### Pipeline

#### Options

* `multiprocess`: Run the JWST Pipeline with multiprocessing where available. Default: False.
* `nprocesses`: Number of processes to use while multiprocessing. Default: 1.
* `stage1`: dict containing JWST Pipeline Stage 1 options
	- `input_dir`: Directory to find JWST Pipeline Stage 1 uncalibrated data. Defaults to pipeline directory.
	- `output_dir`: Directory to output JWST Pipeline Stage 1 products. Defaults to pipeline directory.
	- `steps`: dict of JWST Pipeline Stage 1 steps and their corresponding options