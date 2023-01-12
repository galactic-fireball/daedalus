#!/bin/sh

rsync -zv run_target.py argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/
rsync -rzv --exclude '__pycache__' badass argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/
rsync -rzv --exclude '__pycache__' mast argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/
rsync -rzv --exclude '__pycache__' --exclude 'crds_cache' pipeline argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/
rsync -rzv --exclude '__pycache__' spec argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/
rsync -rzv --exclude '__pycache__' targets argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/
rsync -rzv --exclude '__pycache__' --exclude 'data_sets' --exclude 'archive' programs argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/
