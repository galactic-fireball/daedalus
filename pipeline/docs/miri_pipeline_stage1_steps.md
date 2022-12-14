# MIRI Pipeline Stage 1

## Group Scale

[Space Telescope Description](https://github.com/spacetelescope/jwst/blob/master/docs/jwst/group_scale/description.rst)

Rescales group data to account for on-board frame averaging that did not use `FRMDIVSR = NFRAMES`. All groups in the exposure are rescaled by `FRMDIVSR/NFRAMES`.


## DQ Init

[Space Telescope Description](https://github.com/spacetelescope/jwst/blob/master/docs/jwst/dq_init/description.rst)

Initialize the Data Quality extension from the mask reference file. The `dq_init` step initializes the pixeldq attribute of the input datamodel using the MASK reference file. For some FGS `exp_types`, initialize the dq attribute of the input model instead. The dq attribute of the MASK model is bitwise OR'd with the pixeldq (or dq) attribute of the input model.


## Saturation

[Space Telescope Description](https://github.com/spacetelescope/jwst/blob/master/docs/jwst/saturation/description.rst)

Flags saturated pixels.
