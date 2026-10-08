#' High-dimensional detection of Landscape Dynamics 2
#'
#' This is the main function of the package.
#'
#' \code{hilandyn2_map()} produces maps relative to landscape dynamics through
#' the segmentation of high-dimensional satellite time series into linear
#' trends.
#' These latter include information from the spatial and spectral domains and
#' are analysed using a modified version of the High-dimensional Trend
#' Segmentation (HiTS) procedure proposed by
#' \insertCite{maeng2019adaptive;textual}{hilandyn2}.
#' The HiTS procedure aims to detect changepoints in a piecewise linear signal
#' where their number and location are unknown. Changes can occur in the
#' intercept, slope or both of linear trends.
#' High-dimensional time series can include single or multiple spectral
#' bands/indices, hereafter referred to as bands.
#' Impulsive noise, i.e. outliers in the time series caused by unmasked clouds,
#' cloud shadows, and other types of noise, are removed using a robust baseline
#' as reference. The \code{noise_zscore} parameter controls the noise filter.
#' Sentinel-2 or Landsat reflectance must come from a
#' [FORCE data cube](https://github.com/davidfrantz/force). The selected folder
#' must be the data-cube root containing `X####_Y####` tile folders with one
#' FORCE Level-3 `BAP` (including PAC-generated BAP), `MED`, or `GEO` multiband
#' composite per time step. If more than one composite product is present,
#' select it with `composite_product`. The `include_adjacent` argument controls
#' whether observations from neighbouring tiles contribute to the analysis at
#' tile boundaries.
#'
#' @param FORCE_datacube Character string. Path to the tiled FORCE data-cube
#'   root.
#' @param ts_ids numeric vector. Unique identifiers of six-digits associated
#'   with each time step in the input time series. Values should contain the
#'   four-digit year followed by a two-digit value associated with the month.
#'   See the example for more details.
#' @param satellite_family Character string selecting either `"Sentinel-2"` or
#'   `"Landsat"` files. This is required because both families may coexist in a
#'   FORCE data cube.
#' @param composite_product Optional character string selecting `"BAP"` (Best
#'   Available Pixel), `"MED"` (medoid), or `"GEO"` (geometric median)
#'   composites. `"BAP"` includes both static and phenology-adaptive FORCE
#'   composites, which share the same BAP product code. 
#'   It is inferred when only one product is present.
#' @param tiles Optional character vector of target FORCE tile names. By
#'   default, every available tile is processed.
#' @param roi_vec Character string. Path to a polygon vector file, such as a
#'   shapefile, defining the region of interest for masking and cropping.
#' @param raster_mask Optional path to a single-layer raster mask. Cells with
#'   `NA` are excluded and cells with any non-`NA` value are retained. Unlike
#'   `roi_vec`, it does not crop the output extent.
#'   If the mask does not fully cover the input rasters, a warning is issued
#'   and its extent is extended for processing with value `1`, retaining cells
#'   outside the original mask extent. Existing `NA` cells remain excluded,
#'   and the original mask file is unchanged. The mask must have the same CRS,
#'   resolution, and origin as the FORCE composites.
#' @param spectral_bands Optional character vector selecting FORCE reflectance
#'   bands for the change analysis. Sentinel-2 names are BLUE, GREEN, RED, 
#'   REDEDGE1, REDEDGE2, REDEDGE3, NIR1, NIR2, SWIR1, and SWIR2; Landsat names 
#'   are BLUE, GREEN, RED, NIR, SWIR1, and SWIR2. `NULL` selects every band and
#'   `character(0)` uses reflectance only as the source for on-the-fly indices.
#' @param spectral_indices Optional character vector naming spectral indices to
#'   compute from the selected FORCE reflectance product during the analysis.
#'   Sentinel-2 supports NDVI, NDRE1, NDRE2,
#'   NDRE3, NDMI, NBR, MSI, TCB, TCG, TCW, TCA, CRSWIR, CRRE, DRS, IRECI, and
#'   NMDI. Landsat supports NDVI, NDMI, NBR, MSI, TCB, TCG, TCW, TCA, DRS, and
#'   NMDI.
#' @param rescale_reflectance Logical; multiply FORCE composite reflectance by
#'   0.0001 before computing indices. This should normally be `TRUE` for FORCE's
#'   integer-scaled Level-3 composites.
#' @param cng_dir Optional numeric vector overriding the automatically inferred
#'   direction (either +1 or -1) of disturbance-related change for the selected
#'   spectral variables. When indices are requested, it may describe either
#'   only `spectral_bands` (index directions are appended automatically) or all
#'   selected bands and indices.
#' @param include_adjacent Logical; use all available touching tiles to populate
#'   the spatial-window extension. If `FALSE`, cells outside each target tile
#'   are filled with `NA`.
#' @param win_side integer. Width (in cells) of the spatial kernel used to
#'   extract input data from rasters. Must be an odd number.
#' @param cell_weights logical. Enable or disable the use of cell-based weights
#'   within the spatial kernel. A similarity measure based on the Spectral Angle
#'   Mapper method is computed for every neighboring cell with respect to the
#'   focal one.
#' @param ts_len_min Minimum number of observations in the input time series.
#' @param gap_len_max Maximum number of consecutive missing observations
#'   allowed.
#' @param noise_zscore Z-score threshold used by the impulsive-noise filter.
#' @param wavelet_filters Character vector of candidate wavelet filters.
#' @param wavelet_level Optional wavelet decomposition level.
#' @param wavelet_filter_penalty Penalty applied during wavelet-filter
#'   selection.
#' @param acf_min Minimum autocorrelation required for seasonal adjustment.
#' @param sd_method character string. Statistics used for estimating the
#'   standard deviation of each time series at the beginning of the HiTS
#'   procedure. It can be either "median_abs_dev" (Median Absolute Deviation;
#'   Hampel, 1974), "mean_abs_dev" (Mean Absolute Deviation; Pinsky and
#'   Klawansky, 2023), "biweight_mid" (Biweight midvariance; Lax, 1985),
#'   "std_dev" (Standard deviation).
#' @param th_min,th_max Lower and optional upper bounds of the threshold search.
#' @param dBIC_min Minimum BIC improvement required to retain a model.
#' @param seg_len_min Optional minimum segment length.
#' @param cpt_max Maximum number of expected changepoints.
#' @param out_path character. The path where output rasters (in \emph{.tif}
#'   format) will be saved. Folders are created recursively if they do not
#'   exist. If \code{NULL} the output is a \code{SpatRaster} object created by
#'   the \pkg{terra} package.
#' @param overwrite logical. Controls whether an existing raster with the same
#'   name is overwritten. When `FALSE`, a readable, non-empty raster is
#'   returned without reprocessing; an unreadable or empty file
#'   produces an error.
#'   Reuse is based on the output filename, without comparing input data,
#'   selected variables, time IDs, thresholds, masks, ROI, or output options.
#'   Use a separate \code{out_path} for each analysis configuration, or set
#'   \code{overwrite = TRUE} after changing the inputs or settings.
#' @param rmse_out,len_seg_out,slo_seg_out,bln_mat_out
#'   Logical switches controlling optional output layers.
#' @param zsc_mat_out,mag_mat_out,sgn_mat_out,est_mat_out
#'   Logical switches controlling optional output layers.
#' @param gdal_options GDAL creation options for output files.
#' @param num_cores Positive integer number of target tiles processed
#'   concurrently.
#'   Parallel processing requires a non-`NULL` `out_path`.
#' @param num_threads Positive integer number of computational threads per tile.
#'   Up to `num_cores * num_threads` computational threads can be active.
#'   Builds without OpenMP support process each tile's C++ kernel serially.
#' @param num_copies Positive integer passed to [terra::writeStart()] when
#'   estimating output block sizes. Higher values request smaller chunks and
#'   reduce memory use. The default is `16L`.
#' @param temp_dir Optional directory for temporary processing files. Defaults
#'   to the current \pkg{terra} temporary directory. Temporary tile directories
#'   are removed after successful processing when `out_path` is supplied.
#' @param debug Logical; save the input data associated with a processing error
#'   to a CSV file in the working directory.
#'   Failed cells retain missing outputs. A summary warning is emitted after
#'   each affected chunk finishes, including the failure count and first error.
#'
#' @return For one selected tile, a \code{SpatRaster}; for multiple tiles, a named
#'   list of \code{SpatRaster} objects. A non-\code{NULL} \code{out_path} is
#'   required for multiple or parallel tile processing. Files are stored beneath
#'   \verb{<out_path>/<tile>/}. Each raster contains the outputs below, with
#'   separate layers for each band and/or time step where applicable. Disabled
#'   optional outputs are omitted. Unprocessable windows contain missing values.
#'   \item{D_MAX_TS}{Time identifier of the disturbance selected by the largest
#'   Euclidean norm of focal-band changes standardized by fitted-series standard
#'   deviations.}
#'   \item{D_MAX_MG}{Median absolute relative magnitude (percent), across finite
#'   band/cell values, of the event selected by \code{D_MAX_TS}. This need not
#'   be the largest relative magnitude among disturbance events.}
#'   \item{D_FRS_TS, D_LST_TS}{Time identifiers of the first and last disturbance
#'   events, respectively.}
#'   \item{D_FRS_MG, D_LST_MG}{Median absolute relative magnitudes (percent) of
#'   the first and last disturbance events, respectively.}
#'   \item{D_NUM}{Number of disturbance events; missing when none are detected.}
#'   \item{G_MAX_TS}{Time identifier of the greening event selected using the
#'   same standardized focal-band norm as \code{D_MAX_TS}.}
#'   \item{G_MAX_MG}{Median absolute relative magnitude (percent) of the event
#'   selected by \code{G_MAX_TS}.}
#'   \item{G_FRS_TS, G_LST_TS}{Time identifiers of the first and last greening
#'   events, respectively.}
#'   \item{G_FRS_MG, G_LST_MG}{Median absolute relative magnitudes (percent) of
#'   the first and last greening events, respectively.}
#'   \item{G_NUM}{Number of greening events; missing when none are detected.}
#'   \item{S_ORIG}{Seasonal autocorrelation score before seasonal adjustment.}
#'   \item{S_ADJ}{Seasonal autocorrelation score of the selected adjustment
#'   candidate, even if that candidate is rejected. Equals \code{S_ORIG} when
#'   adjustment is skipped because the original score is below \code{acf_min}.
#'   Both seasonal scores are missing when \code{wavelet_level = 0}.}
#'   \item{TH_CONST}{Final dimension-normalised threshold constant \eqn{C}.
#'   The hard threshold is \eqn{\lambda=C\sqrt{N}}, where \eqn{N} is the
#'   number of rows in the analysed time-series matrix.}
#'   \item{CPT_ID}{Change type at each time step: 101 (disturbance), 102
#'   (greening), or 9 (other change). Missing outside detected changepoints.}
#'   \item{REL_MAG}{Median absolute relative magnitude (percent) across finite
#'   band/cell values at each changepoint; missing at other time steps.}
#'   \item{RMSE}{Optional root mean square error for each focal-cell band,
#'   comparing the preprocessed input with the fitted signal.}
#'   \item{LEN}{Optional segment length, in time steps, at each segment start;
#'   missing at other time steps.}
#'   \item{SLO}{Optional slope for each focal-cell band at segment starts;
#'   missing at other time steps.}
#'   \item{BLN}{Optional robust noise-filter baseline for each focal-cell band
#'   and time step. Missing where no baseline was computed.}
#'   \item{ZSC}{Optional absolute noise-filter z-score for each focal-cell band
#'   and time step. Missing where no score was computed.}
#'   \item{MAG}{Optional signed magnitude for each focal-cell band at
#'   changepoints: fitted pre-change value minus fitted value at the changepoint.
#'   The pre-change value is taken from the same within-year position in the
#'   previous year when available, otherwise from the preceding time step.}
#'   \item{SGN}{Optional preprocessed input signal for each focal-cell band and
#'   time step, after enabled filtering, gap filling, and seasonal adjustment.}
#'   \item{EST}{Optional fitted signal for each focal-cell band and time step.}
#'
#' @author Donato Morresi, \email{donato.morresi@@gmail.com}
#'
#' @seealso \code{\link{hilandyn2_win}}
#'
#' @references
#' \insertRef{Morresi2024}{hilandyn2}
#'
#'
#' @examples
#' \dontrun{
#' # Use a FORCE Level-3 Sentinel-2 data cube containing GEO composites.
#' # Replace the paths and tile ID with those of your data cube.
#' # This example assumes four composites per year (June to September),
#' # ordered as within-year positions 01 to 04, from 2016 through 2025.
#' time_steps <- as.numeric(paste0(rep(2016:2025, each = 4),
#'                                 sprintf("%02d", rep(1:4, times = 10))))
#' rsout <- hilandyn2_map(
#'   FORCE_datacube = "path/to/Sentinel2_FORCE_datacube",
#'   satellite_family = "Sentinel-2",
#'   composite_product = "GEO",
#'   tiles = "X0000_Y0000",
#'   spectral_indices = c("NBR", "MSI", "CRSWIR", "CRRE", "NMDI"),
#'   ts_ids = time_steps,
#'   win_side = 3L,
#'   wavelet_level = NULL, # automatic seasonal adjustment
#'   out_path = "path/to/hilandyn2_output",
#'   est_mat_out = FALSE
#' )
#' 
#' # Plot disturbance magnitude and its corresponding time identifier.
#' terra::plot(rsout[["D_MAX_MG"]])
#' terra::plot(rsout[["D_MAX_TS"]])
#' }
#'
#' @export
#' @import Rcpp
#' @import RcppArmadillo
#' @import terra
#' @importFrom lubridate seconds_to_period
#' @importFrom Rdpack reprompt

hilandyn2_map <- function(
  FORCE_datacube,
  ts_ids,
  satellite_family = c("Sentinel-2", "Landsat"),
  composite_product = NULL,
  tiles = NULL,
  roi_vec = NULL,
  raster_mask = NULL,
  spectral_bands = NULL,
  spectral_indices = NULL,
  rescale_reflectance = TRUE,
  cng_dir = NULL,
  include_adjacent = TRUE,
  win_side = 3L,
  cell_weights = TRUE,
  ts_len_min = 16L,
  gap_len_max = 2L,
  noise_zscore = 3,
  wavelet_filters = c("LA4", "LA6", "LA8"),
  wavelet_level = NULL,
  wavelet_filter_penalty = 0.1,
  acf_min = 0.2,
  sd_method = "mean_abs_dev",
  th_min = 1,
  th_max = NULL,
  dBIC_min = 2,
  seg_len_min = NULL,
  cpt_max = 3,
  out_path = NULL,
  overwrite = FALSE,
  rmse_out = TRUE,
  len_seg_out = FALSE,
  slo_seg_out = FALSE,
  bln_mat_out = FALSE,
  zsc_mat_out = FALSE,
  mag_mat_out = FALSE,
  sgn_mat_out = FALSE,
  est_mat_out = TRUE,
  gdal_options = c("COMPRESS=ZSTD"),
  num_cores = 1L,
  num_threads = 1L,
  num_copies = 16L,
  temp_dir = NULL,
  debug = FALSE
) {

  map_arguments <- list(
    roi_vec = roi_vec,
    raster_mask = raster_mask,
    rmse_out = rmse_out,
    len_seg_out = len_seg_out,
    slo_seg_out = slo_seg_out,
    bln_mat_out = bln_mat_out,
    zsc_mat_out = zsc_mat_out,
    mag_mat_out = mag_mat_out,
    sgn_mat_out = sgn_mat_out,
    est_mat_out = est_mat_out,
    cell_weights = cell_weights,
    ts_len_min = ts_len_min,
    gap_len_max = gap_len_max,
    sd_method = sd_method,
    cng_dir = cng_dir,
    noise_zscore = noise_zscore,
    wavelet_filters = wavelet_filters,
    wavelet_level = wavelet_level,
    wavelet_filter_penalty = wavelet_filter_penalty,
    acf_min = acf_min,
    th_min = th_min,
    th_max = th_max,
    dBIC_min = dBIC_min,
    seg_len_min = seg_len_min,
    cpt_max = cpt_max,
    debug = debug,
    gdal_options = gdal_options,
    rescale_reflectance = rescale_reflectance
  )
  .hilandyn2_map_tiled(
    FORCE_datacube = FORCE_datacube,
    satellite_family = satellite_family,
    composite_product = composite_product,
    out_path = out_path,
    temp_dir = temp_dir,
    ts_ids = ts_ids,
    tiles = tiles,
    spectral_bands = spectral_bands,
    spectral_indices = spectral_indices,
    include_adjacent = include_adjacent,
    num_cores = num_cores,
    win_side = win_side,
    num_threads = num_threads,
    num_copies = num_copies,
    overwrite = overwrite,
    map_arguments = map_arguments
  )
}

.hilandyn2_map_tile <- function(
  FORCE_datacube,
  satellite_family,
  composite_product,
  out_path = NULL,
  roi_vec = NULL,
  raster_mask = NULL,
  pad_focal_edges = FALSE,
  overwrite = FALSE,
  rmse_out = TRUE,
  len_seg_out = FALSE,
  slo_seg_out = FALSE,
  bln_mat_out = FALSE,
  zsc_mat_out = FALSE,
  mag_mat_out = FALSE,
  sgn_mat_out = FALSE,
  est_mat_out = TRUE,
  cell_weights = TRUE,
  win_side = 3L,
  ts_ids,
  ts_len_min = 16L,
  gap_len_max = 2L,
  sd_method = "mean_abs_dev",
  cng_dir = NULL,
  noise_zscore = 3,
  wavelet_filters = c("LA4", "LA6", "LA8"),
  wavelet_level = NULL,
  wavelet_filter_penalty = 0.1,
  acf_min = 0.2,
  th_min = 1,
  th_max = NULL,
  dBIC_min = 2,
  seg_len_min = NULL,
  cpt_max = 3,
  num_threads = 1L,
  num_copies = 16L,
  debug = FALSE,
  gdal_options = c("COMPRESS=ZSTD"),
  spectral_bands = NULL,
  spectral_indices = NULL,
  rescale_reflectance = TRUE,
  analysis_spec = NULL
) {

  num_copies <- as.integer(num_copies)
  if (length(num_copies) != 1L || is.na(num_copies) || num_copies < 1L) {
    stop("num_copies must be one positive integer", call. = FALSE)
  }
  satellite_family <- match.arg(
    satellite_family,
    choices = c(
      "Sentinel-2",
      "Landsat"
    )
  )

  if (win_side%%2 == 0) {
    stop("window sides must be odd")
  }
  if (!is.logical(pad_focal_edges) ||
      length(pad_focal_edges) != 1L || is.na(pad_focal_edges)) {
    stop(
      "pad_focal_edges must be TRUE or FALSE",
      call. = FALSE
    )
  }

  if (is.null(ts_ids) || length(ts_ids) < ts_len_min) {
    stop(
      "ts_ids must contain at least ",
      ts_len_min,
      " time steps"
    )
  }

  nt <- length(ts_ids)

  if (is.null(analysis_spec)) {
    analysis_spec <- .force_analysis_spec(
      satellite_family,
      spectral_bands,
      spectral_indices
    )
  }

  message("\nReading FORCE reflectance data ...")
  sr_source <- .force_map_raster(
    FORCE_datacube,
    satellite_family,
    composite_product,
    nt
  )
  force_spec <- .force_time_series_spec(
    sr_source,
    nt,
    analysis_spec$family
  )
  selected_names <- analysis_spec$band_names
  selected_positions <- analysis_spec$band_positions
  index_names <- analysis_spec$index_names
  index_ids <- analysis_spec$index_ids

  canonical_groups <- force_spec$groups[, force_spec$positions, drop = FALSE]
  
  if (length(index_names)) {
    rs <- .select_grouped_layers(
      sr_source,
      force_spec$groups,
      force_spec$positions,
      force_spec$band_names
    )
  } else {
    rs <- do.call(
      c, lapply(
        seq_len(nt),
        function(i) {
          sr_source[[canonical_groups[i, selected_positions]]]
        }
      )
    )
    names(rs) <- rep(
      selected_names,
      times = nt
    )
  }

  first_step_names <- analysis_spec$variable_names
  external_nb <- analysis_spec$external_variables
  computed_nb <- analysis_spec$computed_variables

  if (!is.null(roi_vec)) {
    roi_vec <- vect(roi_vec)

    if (crs(roi_vec) != crs(rs)) {
      roi_vec <- project(roi_vec, crs(rs))
    }

    if (relate(
      ext(rs),
      roi_vec,
      "covers"
    )) {
      rs <- crop(rs, roi_vec)
    } else {
      stop("input rasters do not cover roi_vec completely")
    }

    message("\nMasking input rasters using vector ROI ...")
    rs <- mask(rs, roi_vec)
  }

  pad_processing_edges <- pad_focal_edges ||
    !is.null(roi_vec)
  raster_mask_spec <- if (is.null(raster_mask)) {
    NULL
  } else {
    .map_raster_mask_spec(
      terra::rast(raster_mask),
      rs
    )
  }

  nb <- length(first_step_names)
  n_lyrs <- nb * nt

  cng_dir <- .resolve_change_directions(
    first_step_names,
    cng_dir,
    external_nb,
    computed_nb
  )

  # Parameters for the focal computation
  nc <- win_side * win_side
  hw <- (win_side - 1)/2
  nr <- nb * nc
  cng_dir <- rep(cng_dir, each = nc)

  # Output layer names
  valn <- c(
    "D_MAX_TS",
    "D_MAX_MG",
    "D_FRS_TS",
    "D_FRS_MG",
    "D_LST_TS",
    "D_LST_MG",
    "D_NUM",
    "G_MAX_TS",
    "G_MAX_MG",
    "G_FRS_TS",
    "G_FRS_MG",
    "G_LST_TS",
    "G_LST_MG",
    "G_NUM",
    "S_ORIG",
    "S_ADJ",
    "TH_CONST"
  )

  vec1n <- c()  # No values included by default
  if (rmse_out)
    vec1n <- c("RMSE")

  vec2n <- c("CPT_ID", "REL_MAG")
  if (len_seg_out)
    vec2n <- c(vec2n, "LEN")

  matn <- c()  # No values included by default
  if (slo_seg_out)
    matn <- c(matn, "SLO")
  if (bln_mat_out)
    matn <- c(matn, "BLN")
  if (zsc_mat_out)
    matn <- c(matn, "ZSC")
  if (mag_mat_out)
    matn <- c(matn, "MAG")
  if (sgn_mat_out)
    matn <- c(matn, "SGN")
  if (est_mat_out)
    matn <- c(matn, "EST")

  # Set output names and number of output layers
  out_nm <- valn
  if (length(vec1n) > 0)
    out_nm <- c(
      out_nm, paste0(
        rep(vec1n, each = nb),
        "_",
        seq(1, nb)
      )
    )
  if (length(vec2n) > 0)
    out_nm <- c(
      out_nm, paste0(
        rep(vec2n, each = nt),
        "_",
        ts_ids
      )
    )
  if (length(matn) > 0)
    out_nm <- c(
      out_nm, paste0(
        rep(
          matn, each = nb * nt
        ),
        "_",
        seq(1, nb),
        "_",
        rep(ts_ids, each = nb)
      )
    )

  nout <- length(out_nm)

  # Assign values to NULL parameters
  if (is.null(th_max))
    th_max <- -1
  if (is.null(seg_len_min))
    seg_len_min <- -1
  if (is.null(wavelet_level))
    wavelet_level <- -1

  # Set output parameters
  if (is.null(out_path)) {
    filename = ""
    wopt = list()
  } else {
    if (!dir.exists(out_path))
      dir.create(
        out_path,
        recursive = TRUE
      )
    filename = .map_output_filename(
      out_path,
      nb,
      nc,
      ts_ids
    )
    wopt = list(
      filetype = "GTiff",
      datatype = "FLT4S",
      NAflag = -9999,
      gdal = gdal_options
    )
  }

  if (pad_processing_edges) {
    out <- rast(rs, nlyrs = nout)
  } else {
    rs_extent <- ext(rs)
    rs_resolution <- res(rs)
    output_extent <- ext(
      xmin(rs_extent) + hw * rs_resolution[1],
      xmax(rs_extent) - hw * rs_resolution[1],
      ymin(rs_extent) + hw * rs_resolution[2],
      ymax(rs_extent) - hw * rs_resolution[2]
    )
    out <- rast(
      nrows = nrow(rs) - 2 * hw,
      ncols = ncol(rs) - 2 * hw,
      nlyrs = nout,
      crs = crs(rs),
      extent = output_extent
    )
  }

  set.names(out, out_nm)
  n_cols <- if (pad_processing_edges) 
    ncol(rs) + 2L * hw else ncol(rs)
  readStart(rs)
  on.exit(
    readStop(rs),
    add = TRUE
  )
  if (!is.null(raster_mask_spec)) {
    readStart(raster_mask_spec$raster)
    on.exit(
      readStop(raster_mask_spec$raster),
      add = TRUE
    )
  }

  b <- writeStart(
    out,
    filename,
    overwrite,
    n = num_copies,
    wopt = wopt
  )
  output_open <- TRUE
  on.exit({
    if (output_open) try(writeStop(out), silent = TRUE)
  }, add = TRUE)
  tic <- proc.time()

  for (i in seq_len(b$n)) {
    message(
      "\nProcessing chunk #",
      i,
      "/",
      b$n
    )

    str_row <- b$row[i]

    if (pad_processing_edges) {
      padded_chunk <- .map_read_padded_chunk(
        rs,
        str_row,
        b$nrows[i],
        hw
      )
      input_values <- padded_chunk$values
      n_rows <- padded_chunk$n_rows
    } else {
      n_rows <- b$nrows[i] +
        2 * hw
      input_values <- rs@pntr$readValues(
        str_row - 1,
        n_rows,
        0,
        n_cols
      )
    }
    mask_values <- NULL
    if (!is.null(raster_mask_spec)) {
      mask_values <- .map_read_mask_chunk(
        raster_mask_spec,
        str_row,
        n_rows,
        n_cols,
        half_window = if (pad_processing_edges) hw else 0L
      )
    }
    if (length(index_names)) {
      input_values <- .force_map_chunk_values(
        input_values,
        n_rows * n_cols,
        terra::nlyr(rs),
        nt,
        selected_positions,
        index_ids,
        force_spec$id,
        rescale_reflectance,
        num_threads,
        mask_values
      )
    } else if (!is.null(mask_values)) {
      input_values <- mask_chunk_values_cpp(
        input_values,
        mask_values,
        n_rows * n_cols,
        terra::nlyr(rs)
      )
    }

    m <- hilandyn2_map_cpp(
      input_values,
      c(
        n_rows,
        n_cols,
        n_lyrs
      ),
      win_side,
      nb,
      nc,
      nt,
      nr,
      ts_ids,
      cell_weights,
      rmse_out,
      len_seg_out,
      slo_seg_out,
      bln_mat_out,
      zsc_mat_out,
      mag_mat_out,
      sgn_mat_out,
      est_mat_out,
      nout,
      ts_len_min,
      gap_len_max,
      sd_method,
      cng_dir,
      noise_zscore,
      wavelet_filters,
      wavelet_level,
      wavelet_filter_penalty,
      acf_min,
      th_min,
      th_max,
      dBIC_min,
      seg_len_min,
      cpt_max,
      num_threads,
      debug
    )

    writeValues(
      out,
      m,
      str_row,
      b$nrows[i]
    )
  }

  toc <- proc.time()
  elapsed <- round(seconds_to_period((toc - tic)[3]))
  message(
    "\nProcessing took ", elapsed
  )

  out <- writeStop(out)
  output_open <- FALSE
  gc(verbose = FALSE)
  return(invisible(out))
}

.map_read_padded_chunk <- function(
  x,
  output_row,
  output_rows,
  half_window,
  source_rows = terra::nrow(x),
  source_columns = terra::ncol(x),
  row_offset = 0L,
  column_offset = 0L
) {
  layers <- terra::nlyr(x)
  padded_rows <- output_rows + 2L * half_window
  padded_columns <- source_columns + 2L * half_window

  desired_first_row <- output_row - half_window
  desired_last_row <- output_row + output_rows - 1L +
    half_window
  first_source_row <- max(1L, desired_first_row)
  last_source_row <- min(source_rows, desired_last_row)
  rows_to_read <- last_source_row - first_source_row +
    1L

  source_values <- x@pntr$readValues(
    row_offset + first_source_row - 1L,
    rows_to_read,
    column_offset,
    source_columns
  )
  padded_values <- matrix(
    NA_real_,
    nrow = padded_rows * padded_columns,
    ncol = layers
  )

  destination_rows <- first_source_row - desired_first_row +
    seq_len(rows_to_read)
  destination_cells <- unlist(
    lapply(
      destination_rows, function(row) {
        (row - 1L) * padded_columns + half_window +
          seq_len(source_columns)
      }
    ),
    use.names = FALSE
  )
  padded_values[destination_cells, ] <- source_values
  dim(padded_values) <- NULL

  list(
    values = padded_values,
    n_rows = padded_rows,
    n_columns = padded_columns
  )
}

.map_raster_mask_spec <- function(raster_mask, x) {
  if (terra::nlyr(raster_mask) != 1L) {
    stop("raster_mask must contain exactly one layer", call. = FALSE)
  }
  same_origin <- isTRUE(
    all.equal(
      terra::origin(raster_mask),
      terra::origin(x),
      tolerance = .Machine$double.eps^0.5
    )
  )
  same_grid <- terra::compareGeom(
    raster_mask,
    x,
    crs = TRUE,
    ext = FALSE,
    rowcol = FALSE,
    res = TRUE,
    stopOnError = FALSE
  )
  if (!same_grid || !same_origin) {
    stop(
      "raster_mask and FORCE composites must have the same CRS, ",
      "resolution, and origin",
      call. = FALSE
    )
  }
  resolution <- terra::res(x)
  row_offset <- round(
    (terra::ymax(raster_mask) -
      terra::ymax(x))/resolution[2L]
  )
  column_offset <- round(
    (terra::xmin(x) -
      terra::xmin(raster_mask))/resolution[1L]
  )
  mask_rows <- terra::nrow(raster_mask)
  mask_columns <- terra::ncol(raster_mask)
  if (row_offset < 0L || column_offset < 0L || row_offset + terra::nrow(x) > mask_rows || column_offset + terra::ncol(x) > mask_columns) {
    warning(
      "raster_mask does not fully cover the input rasters; ",
      "its extent is extended for processing with value 1 (unmasked) ",
      "outside the original mask extent",
      call. = FALSE
    )
  }

  list(
    raster = raster_mask,
    row_offset = as.integer(row_offset),
    column_offset = as.integer(column_offset),
    n_rows = mask_rows,
    n_columns = mask_columns
  )
}

.map_read_mask_chunk <- function(
  specification,
  row,
  nrows,
  ncols,
  half_window = 0L
) {

  mask_row <- specification$row_offset +
    row - 1L - half_window
  mask_column <- specification$column_offset -
    half_window
  first_row <- max(0L, mask_row)
  last_row <- min(
    specification$n_rows, mask_row +
      nrows
  )
  first_column <- max(0L, mask_column)
  last_column <- min(
    specification$n_columns, mask_column +
      ncols
  )
  rows_to_read <- last_row - first_row
  columns_to_read <- last_column - first_column

  if (rows_to_read <= 0L || columns_to_read <= 0L) {
    return(rep.int(1, nrows * ncols))
  }
  source_values <- specification$raster@pntr$readValues(
    first_row,
    rows_to_read,
    first_column,
    columns_to_read
  )
  if (rows_to_read == nrows && columns_to_read == ncols) {
    return(source_values)
  }

  values <- matrix(
    1,
    nrow = ncols,
    ncol = nrows
  )
  values[first_column - mask_column +
    seq_len(columns_to_read),
    first_row - mask_row + seq_len(rows_to_read)] <- source_values
  dim(values) <- NULL
  values
}
