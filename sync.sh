#!/bin/sh

rsync -zv badass/badass_miri.py argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/badass/badass_miri.py
rsync -zv badass/miri_consts.py argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/badass/miri_consts.py
rsync -zv badass/ir_options1.py argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/badass/ir_options1.py
rsync -zv badass/clear_undone_pixels.py argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/badass/clear_undone_pixels.py

rsync -zv mast/mast.py argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/mast/mast.py

rsync -zv pipeline/miri.py argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/pipeline/miri.py
rsync -zv pipeline/nirspec.py argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/pipeline/nirspec.py
rsync -zv pipeline/crds_utils.py argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/pipeline/crds_utils.py
rsync -zv pipeline/run_crds_cache.slurm argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/pipeline/run_crds_cache.slurm
rsync -zv pipeline/plot_specs.py argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/pipeline/plot_specs.py
rsync -zv pipeline/spec_common.py argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/pipeline/spec_common.py

rsync -zv programs/miri_program.py argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/programs/miri_program.py
rsync -zv programs/nirspec_program.py argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/programs/nirspec_program.py
rsync -zv programs/run_program.py argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/programs/run_program.py

rsync -zv programs/1328/run_all_lines.sh argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/programs/1328/run_all_lines.sh
rsync -zv programs/1328/ngc_7469_pipeline_1.8.2.toml argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/programs/1328/ngc_7469_pipeline_1.8.2.toml
rsync -zv programs/1328/ngc_7469_NIRSPECIFU_pipeline_1.8.5.toml argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/programs/1328/ngc_7469_NIRSPECIFU_pipeline_1.8.5.toml
rsync -zv programs/1328/run_NGC7469.slurm argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/programs/1328/run_NGC7469.slurm
rsync -zv programs/1328/run_NGC7469_NSIFU.slurm argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/programs/1328/run_NGC7469_NSIFU.slurm
rsync -zv programs/1328/VV114_pipeline_1.8.4.toml argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/programs/1328/VV114_pipeline_1.8.4.toml
rsync -zv programs/1328/run_VV114.slurm argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/programs/1328/run_VV114.slurm
rsync -zv programs/1328/IIZw096_pipeline_1.8.4.toml argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/programs/1328/IIZw096_pipeline_1.8.4.toml
rsync -zv programs/1328/run_IIZw096.slurm argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/programs/1328/run_IIZw096.slurm

rsync -zv programs/2732/run_all_lines.sh argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/programs/2732/run_all_lines.sh
rsync -zv programs/2732/ngc_7319_pipeline_1.8.2.toml argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/programs/2732/ngc_7319_pipeline_1.8.2.toml
rsync -zv programs/2732/run_NGC7319.slurm argo:/projects/ssatyapa/spectra/sdoan2/jwst/jwst-utils/programs/2732/run_NGC7319.slurm
