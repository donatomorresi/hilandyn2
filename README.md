
<h1>
High-dimensional detection of Landscape Dynamics 2 (HILANDYN2)
<img src="man/figures/hex_logo_hilandyn2.png" align="right" width="120" alt="hilandyn2 hex logo" />
</h1>

<!-- badges: start -->

[![GitHub release][release-badge]][latest-release]
[![R build status][build-badge]][build-status]

[release-badge]:
  https://img.shields.io/github/v/release/donatomorresi/hilandyn2
[latest-release]: https://github.com/donatomorresi/hilandyn2/releases/latest
[build-badge]:
  https://github.com/donatomorresi/hilandyn2/actions/workflows/R-CMD-check.yaml/badge.svg
[build-status]: https://github.com/donatomorresi/hilandyn2/actions/workflows/R-CMD-check.yaml
<!-- badges: end -->

`hilandyn2` is an *R* package that implements the High-dimensional detection of 
Landscape Dynamics 2 framework (HILANDYN2; Morresi *et al.* submitted) for 
mapping stand-replacing and non-stand replacing forest disturbance by temporally 
segmenting high-dimensional optical time series. 
Compared to HILANDYN (Morresi *et al.* 2024), HILANDYN2 supports intra-annual 
satellite time series and introduces improved gap-filling and noise detection. 
Image composites (Level 3) from both Landsat and Sentinel-2 imagery are supported. 
Compatible Level 3 data can be obtained using 
[FORCE](https://github.com/davidfrantz/force) or `hilandyn2_composite()`.

High-dimensional optical time series include information from the
spatial and spectral dimensions. Time series matrices at the pixel level are 
created by extracting spectral variables time series from pixels within a 
spatial kernel. Spectral variables can include multiple original bands and 
spectral indices. `hilandyn2` uses a modified version of the High-dimensional 
Trend Segmentation (HiTS) procedure proposed by Maeng (2019). The HiTS procedure 
aims to detect changepoints in a piecewise linear signal where the number and 
location of changepoints are unknown. Changes can occur in the intercept, slope 
or both of linear trends. The HiTS procedure implemented in `hilandyn2` uses an 
early stopping criterion based on a modified BIC to find the optimal time series 
segmentation threshold. 

`hilandyn2_composite()` provides pixel-wise image compositing methods based on 
the medoid (Flood, 2013) and geometric median (Roberts *et al.*, 2017).
Observation-based weights can be used to reduce intra-annual spectral 
variability. Gap-filling, noise filtering, and seasonal adjustment are performed 
before temporal segmentation. Seasonal adjustment leverages the maximal overlap 
discrete wavelet packet transform (MODWPT; Percival & Walden, 2000). 
      

## Installation

You can install the latest code from the default branch of `hilandyn2` on
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("donatomorresi/hilandyn2")
```

## Reusing output files

With `overwrite = FALSE`, `hilandyn2_map()` reuses readable, non-empty output
rasters by filename. `hilandyn2_composite()` reuses existing products by filename
and regenerates zero-byte files. Neither function compares the inputs or
processing settings with those used to create an existing output.

Use a separate `out_path` (mapping) or `out_file` root (compositing) for each
analysis configuration. When changing inputs, variables, dates, thresholds,
masks, binning, weighting, or output options in the same directory, use
`overwrite = TRUE` to regenerate the results.

### References

Morresi, D., Maeng, H., Marzano, R., Lingua, E., Motta, R., Garbarino, M., 2024. 
High-dimensional detection of Landscape Dynamics: a Landsat time series-based 
algorithm for forest disturbance mapping and beyond. 
GIScience Remote Sens. 61. https://doi.org/10.1080/15481603.2024.2365001

Maeng, H., 2019. Adaptive Multiscale Approaches to Regression and Trend 
Segmentation. The London School of Economics and Political Science (LSE).

Flood, N., 2013. Seasonal composite Landsat TM/ETM+ Images using the medoid 
(a multi-dimensional median). Remote Sens. 5, 6481–6500. 
https://doi.org/10.3390/rs5126481

Roberts, D., Mueller, N., McIntyre, A., 2017. High-Dimensional Pixel Composites 
From Earth Observation Time Series. IEEE Trans. Geosci. Remote Sens. 
55, 6254–6264. https://doi.org/10.1109/TGRS.2017.2723896

Percival, D.B., Walden, A.T., 2000. Wavelet Methods for Time Series Analysis. 
Cambridge University Press. https://doi.org/10.1017/cbo9780511841040
