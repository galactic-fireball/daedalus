#!/bin/sh

rsync -zv run_target.py hopper:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/
rsync -rzv --exclude '__pycache__' badass hopper:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/
rsync -rzv --exclude '__pycache__' --exclude 'jwst_all.csv' mast hopper:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/
rsync -rzv --exclude '__pycache__' --exclude 'crds_cache' pipeline hopper:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/
rsync -rzv --exclude '__pycache__' --exclude 'webbpsf-data' instruments hopper:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/
rsync -rzv --exclude '__pycache__' --exclude 'data_sets' --exclude 'archive' --exclude 'test_data' programs hopper:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/
rsync -rzv --exclude '__pycache__' utils hopper:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/
rsync -rzv --exclude '__pycache__' --exclude 'webbpsf' data hopper:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/
