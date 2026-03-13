## `product_name`
*Type:* `str`<br/>
*Default:* `None [Required Field]`<br/>
*Description:* Unique product name for this build. It can be anything, but will be used to refer to this build.

## `target`
*Type:* `TargetConfig`<br/>
*Default:* `None [Required Field]`<br/>
*Description:* Configuration for observed target of this build.

## `target.name`
*Type:* `str`<br/>
*Default:* `None [Required Field]`<br/>
*Description:* Name of target. This does not have to be any particularly referenced name, it is for your benefit.

## `target.program`
*Type:* `int`<br/>
*Default:* `None [Required Field]`<br/>
*Description:* The program ID of the target

## `target.instrument`
*Type:* `NIRSpecIFUInstrumentConfig | MIRIInstrumentConfig`<br/>
*Default:* `None [Required Field]`<br/>
*Description:* Configuration of the instrument and observation of the target for this build.

## `target.NIRSpecIFUInstrumentConfig.name`
*Type:* `nirspec_ifu`<br/>
*Default:* `nirspec_ifu`<br/>
*Options:* `nirspec_ifu`

## `target.NIRSpecIFUInstrumentConfig.obs`
*Type:* `int`<br/>
*Default:* `None [Required Field]`<br/>
*Description:* The observation number for the target found in the program PDF.

## `target.MIRIInstrumentConfig.name`
*Type:* `miri`<br/>
*Default:* `miri`<br/>
*Options:* `miri`

## `target.MIRIInstrumentConfig.sci_obs`
*Type:* `int`<br/>
*Default:* `None [Required Field]`<br/>
*Description:* The science observation number for the target found in the program PDF.

## `target.MIRIInstrumentConfig.bkgd_obs`
*Type:* `int`<br/>
*Default:* `None [Required Field]`<br/>
*Description:* The background observation number for the target found in the program PDF.

## `actions`
*Type:* `list`<br/>
*Default:* `None [Required Field]`<br/>
*Description:* The actions to be performed during this build.

## `DownloadConfig.name`
*Type:* `download`<br/>
*Default:* `download`<br/>
*Options:* `download`

## `DownloadConfig.overwrite`
*Type:* `bool`<br/>
*Default:* `False`<br/>
*Description:* Overwrite existing files.

## `DownloadConfig.mast_token`
*Type:* `str | NoneType`<br/>
*Default:* `None`<br/>
*Description:* MAST API token. Only needed to download proprietary data that you have access to.

## `PipelineConfig.name`
*Type:* `pipeline`<br/>
*Default:* `pipeline`<br/>
*Options:* `pipeline`

## `PipelineConfig.overwrite`
*Type:* `bool`<br/>
*Default:* `False`<br/>
*Description:* Overwrite existing files.

## `PipelineConfig.multiprocess`
*Type:* `bool`<br/>
*Default:* `False`<br/>
*Description:* Switch to run the pipeline build in multiprocess mode.

## `PipelineConfig.nprocesses`
*Type:* `int`<br/>
*Default:* `1`<br/>
*Description:* The number of processes to use when running in multiprocess mode.

## `PipelineConfig.stage1`
*Type:* `PipelineStage1Options`<br/>
*Default:* `stage=1 steps={} flag_snowballs=False clean_rates=False`<br/>
*Description:* Options for Stage 1 of the pipeline build.

## `PipelineConfig.stage1.stage`
*Type:* `1 | 2 | 3`<br/>
*Default:* `1`<br/>
*Options:* `1`, `2`, `3`

## `PipelineConfig.stage1.steps`
*Type:* `dict`<br/>
*Default:* `None [Required Field]`<br/>
*Description:* Additional options to send to specific steps in this stage.

## `PipelineConfig.stage1.flag_snowballs`
*Type:* `bool`<br/>
*Default:* `False`<br/>
*Description:* TODO

## `PipelineConfig.stage1.clean_rates`
*Type:* `bool`<br/>
*Default:* `False`<br/>
*Description:* TODO

## `PipelineConfig.stage2`
*Type:* `PipelineStage2Options`<br/>
*Default:* `stage=2 steps={}`<br/>
*Description:* Options for Stage 2 of the pipeline build.

## `PipelineConfig.stage2.stage`
*Type:* `2`<br/>
*Default:* `2`<br/>
*Options:* `2`

## `PipelineConfig.stage2.steps`
*Type:* `dict`<br/>
*Default:* `None [Required Field]`<br/>
*Description:* Additional options to send to specific steps in this stage.

## `PipelineConfig.stage3`
*Type:* `PipelineStage3Options`<br/>
*Default:* `stage=3 steps={} outlier_detection=False`<br/>
*Description:* Options for Stage 3 of the pipeline build.

## `PipelineConfig.stage3.stage`
*Type:* `3`<br/>
*Default:* `3`<br/>
*Options:* `3`

## `PipelineConfig.stage3.steps`
*Type:* `dict`<br/>
*Default:* `None [Required Field]`<br/>
*Description:* Additional options to send to specific steps in this stage.

## `PipelineConfig.stage3.outlier_detection`
*Type:* `bool`<br/>
*Default:* `False`<br/>
*Description:* TODO

## `uncal_dir`
*Type:* `Path | str`<br/>
*Default:* `None`<br/>
*Description:* Directory of uncalibrated data. If one is not provided, Daedalus will create one for you and download the necessary files.

## `stage2_dir`
*Type:* `Path | str`<br/>
*Default:* `None`<br/>
*Description:* Directory of data processed by stage 1, ready for stage 2. If one is not provided, Daedalus will create one for you and download or build the necessary files.

## `stage3_dir`
*Type:* `Path | str`<br/>
*Default:* `None`<br/>
*Description:* Directory of data processed by stage 2, ready for stage 3. If one is not provided, Daedalus will create one for you and download or build the necessary files.

## `output_dir`
*Type:* `Path | str`<br/>
*Default:* `None [Required Field]`<br/>
*Description:* Daedalus build output directory. Note that the final output directory will be: <output_dir>/<program_id>/<target_name>/<product_name>

## `crds_cache`
*Type:* `Path | str`<br/>
*Default:* `None`<br/>
*Description:* CRDS cache directory. Will default to some directory as determined by the CRDS module. Warning! These directories can get up to 10s of GBs.

