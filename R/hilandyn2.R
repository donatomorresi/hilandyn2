#' hilandyn2: Implementation of the High-dimensional detection of Landscape
#' Dynamics (HILANDYN) 2 algorithm for mapping forest disturbance dynamics
#' through the temporal segmentation of high-dimensional intra-annual time
#' series.
#'
#' HILANDYN employs the High-dimensional Trend Segmentation (HiTS) procedure
#' proposed by Maeng (2019) to detect changepoints in linear trends from Landsat
#' time series including information in the spatial and spectral dimensions.
#' Optical Sentinel-2 and Landsat data used by the compositing and on-the-fly
#' spectral-index modules must be FORCE Analysis Ready Data from a
#' [FORCE data cube](https://github.com/davidfrantz/force).
#' The HiTS procedure aims to detect changepoints in a piecewise linear signal
#' where their number and location are unknown. Changes can occur in the
#' intercept, slope or both of linear trends.
#' To start with, see the function \code{hilandyn2_map}.
#'
#' @author Donato Morresi \email{donato.morresi@@gmail.com}, Hyeyoung Maeng
#'   \email{hyeyoung.maeng@@durham.ac.uk}
#' @references Morresi, D., Maeng, H., Marzano, R., Lingua, E., Motta, R.,
#'   Garbarino, M., 2024. High-dimensional detection of Landscape Dynamics: a
#'   Landsat time series-based algorithm for forest disturbance mapping and
#'   beyond. GIScience Remote Sens. 61.
#'   https://doi.org/10.1080/15481603.2024.2365001
#' @seealso \code{\link{hilandyn2_map}}, \code{\link{hilandyn2_win}}
#' @name hilandyn2
#' @useDynLib hilandyn2
NULL
