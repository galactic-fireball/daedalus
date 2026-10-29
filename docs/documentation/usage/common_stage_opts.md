---
title: Common Stage Options
---

Last updated for JWST version 3.0.0


# Stage 1

## `clean_flicker_noise`

Applies 1/f noise correction. **Skipped by default.**

To turn on:
```
[[actions]]
name = 'pipeline'
stage1.clean_noise = true
```

---

# Stage 2

## `clean_flicker_noise`

Applies 1/f noise correction. **Skipped by default.**

To turn on:
```
[[actions]]
name = 'pipeline'
stage2.clean_noise = true
```


## `residual_fringe`

Applies residual fringe correction. For MIRI MRS only. **Skipped by default.**

To turn on:
```
[[actions]]
name = 'pipeline'
stage2.rf_corr = true
```


## `adaptive_trace_model`

**Skipped by default.**

To turn on:
```
[[actions]]
name = 'pipeline'
stage2.atm = true
```

---

# Stage 3

## `adaptive_trace_model`

**Skipped by default.**

To turn on:
```
[[actions]]
name = 'pipeline'
stage3.atm = true
```


## `cube_build` pixel scale

Pixel scale of final cube, in arcseconds. Default is 0.1"

To modify:
```
[[actions]]
name = 'pipeline'
stage3.cube_scale = 0.05
```
