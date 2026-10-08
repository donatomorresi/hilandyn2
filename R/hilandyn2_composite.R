#' Produce FORCE optical image composites
#'
#' Creates temporal composites from FORCE Level-2 bottom-of-atmosphere (BOA)
#' reflectance using either the medoid or geometric median of the valid
#' observations in each bin. Input rasters must come from a
#' [FORCE data cube](https://github.com/davidfrantz/force). Sentinel-2 BOA must
#' contain its ten land-surface bands; Landsat 4/5/7/8/9 BOA must contain its
#' six legacy bands. `satellite_family` selects either the `SEN` or `LND` FORCE
#' acquisitions, allowing both families to coexist in a data cube.
#' Acquisition dates are read from the leading `YYYYMMDD` component of each BOA
#' filename. Only acquisitions within `start_date` and `end_date` are included.
#' Bins without valid observations remain empty (`NA`). Each bin target date
#' is its lower temporal midpoint. Consequently, a complete 30-day
#' calendar-month bin is dated on the 15th. 
#' Products follow the FORCE pattern `YYYYMMDD_LEVEL3_SENSOR_PRODUCT.tif` 
#' and are organized in FORCE tile directories. `INF` is always derived from the
#' selected medoid, including when the reflectance product is GEO. Consequently,
#' MED and GEO can coexist in one output data cube and share one
#' method-independent
#' `INF` file. Stable sensor IDs are SEN2A--SEN2D = 1--4 and LND04, LND05,
#' LND07, LND08, LND09 = 1--5.
#'
#' @param FORCE_datacube Path to a FORCE data cube containing
#'   `datacube-definition.prj` and its `X####_Y####` tile directories. A path
#'   directly to one tile remains supported. Each tile must contain matching
#'   `_BOA` and `_QAI` files. BOA bands must use FORCE order:
#'   BLUE, GREEN, RED, REDEDGE1, REDEDGE2, REDEDGE3, NIR1, NIR2, SWIR1, SWIR2
#'   for Sentinel-2; or BLUE, GREEN, RED, NIR, SWIR1, SWIR2 for Landsat. Each
#'   basename must follow the FORCE convention beginning with `YYYYMMDD` and
#'   ending in `_BOA.tif`; its QAI companion must have the same prefix and end
#'   in
#'   `_QAI.tif` (GeoTIFF and VRT extensions are accepted).
#' @param satellite_family Satellite family to composite, either
#'   `"Sentinel-2"` (the default) or `"Landsat"`. FORCE files are selected using
#'   mission tokens beginning with `SEN` or `LND`, respectively, such as
#'   `20150704_LEVEL2_SEN2A_BOA.tif`.
#' @param tiles Optional character vector of FORCE tile names, such as
#'   `"X0069_Y0043"`. Names must correspond to tile directories in
#'   `FORCE_datacube`.
#' @param aoi Optional path to an AOI vector dataset (including an ESRI
#'   Shapefile), or a `SpatVector`. Available tiles intersecting the AOI are
#'   derived from `datacube-definition.prj`. If both `tiles` and `aoi` are
#'   supplied, only explicitly named tiles intersecting the AOI are processed.
#' @param start_date,end_date Bounds of the compositing interval. They default
#'   to
#'   the first and last acquisition dates.
#' @param qai_bits Integer vector of zero-based FORCE QAI bit positions to
#'   reject, between 0 and 15. Defaults to `c(0, 1, 2, 3, 4, 8, 9)`; see
#'   Quality screening below. A supplied vector replaces the defaults;
#'   repeated positions are ignored. Use `integer(0)` to disable flag-based
#'   screening. These are bit positions, not complete decimal QAI values.
#' @param method Composite location, either `"medoid"` or `"geomedian"`.
#' @param use_weights Logical; compute the weighted medoid or weighted
#'   geometric median using cloud-distance and spectral-similarity weights 
#'   derived from the Euclidean distance and Spectral Angle Mapper metrics.
#'   Defaults to `FALSE` for unweighted compositing. When `TRUE`, every BOA
#'   acquisition selected for the requested satellite family and date interval
#'   must have a matching separate FORCE `_DST` file; a missing file causes
#'   an error. When `FALSE`, DST files are not used or required.
#' @param cloud_distance Distance at which the cloud-distance weighting curve is
#'   evaluated, in the units of the FORCE DST product. Used only when
#'   `use_weights = TRUE`.
#' @param bin_method Either `"temporal"` or `"hybrid"`. Temporal bins divide
#'   the requested date interval. Hybrid bins start with balanced numbers of
#'   valid observations and reassign observations that are too far from their
#'   temporal-bin centre.
#' @param n_bins Positive integer number of output bins.
#' @param hybrid_max_days Maximum distance from the initially assigned bin
#'   centre before a hybrid-bin observation is reassigned.
#' @param out_file Optional path to an output FORCE data-cube
#'   root. The input `datacube-definition.prj` is copied to this directory and
#'   products are written beneath `<out_file>/<tile>/`. MED and GEO products
#'   may coexist and share the same medoid-derived `INF` files. If `NULL`,
#'   products are returned without being persisted and `num_cores` must be one.
#' @param overwrite Logical; overwrite existing output products. When `FALSE`,
#'   existing zero-byte MED, GEO, and INF files are removed and regenerated.
#'   Other existing products are reused by filename without comparing input
#'   data, date intervals, binning, weighting, or quality-mask settings.
#'   Use a separate \code{out_file} root for each configuration, or set
#'   \code{overwrite = TRUE} after changing the inputs or settings.
#' @param gdal_options GDAL creation options for an output file.
#' @param num_cores Positive integer number of tiles processed concurrently.
#'   Values
#'   greater than one require multiple selected tiles and a non-`NULL`
#'   `out_file`.
#' @param num_threads Positive integer number of computational threads per tile.
#'   Up to `num_cores * num_threads` computational threads can be active.
#'
#' @return For one tile, a list with `reflectance` and `information`
#'   `SpatRaster` objects. For multiple tiles, a named list of those product
#'   lists. Each bin is a separate FORCE-style file when `out_file` is supplied.
#'   Reflectance filenames use product type `MED` or `GEO`; paired `INF` files
#'   contain QAI, valid-observation count, representative acquisition DOY and
#'   year, difference from target DOY, and representative sensor ID.
#'
#' @section Quality screening:
#' For both Sentinel-2 and Landsat, an observation is excluded from all
#' reflectance bands if any selected QAI bit is set. The defaults reject:
#' \itemize{
#'   \item Bit 0: nodata.
#'   \item Bits 1 and 2: all non-clear cloud states (buffered cloud,
#'     opaque cloud, or cirrus).
#'   \item Bit 3: cloud shadow.
#'   \item Bit 4: snow.
#'   \item Bit 8: subzero reflectance.
#'   \item Bit 9: saturation.
#' }
#' Other flags, including water (bit 5), are not rejected by default.
#' For example, `qai_bits = c(0:5, 8L, 9L)` additionally rejects water.
#' Cloud state is encoded jointly by bits 1 and 2: selecting only one of
#' them tests that bit, not one isolated cloud category. See the
#' [FORCE QAI flags](https://davidfrantz.github.io/tutorials/force-qai/qai/).
#'
#' Non-finite QAI values or non-finite values in any reflectance band always
#' exclude the observation, including when `qai_bits = integer(0)`.
#' Bins with no remaining observations stay `NA`. Use `overwrite = TRUE`
#' to regenerate existing products after changing `qai_bits`. MED and GEO
#' products sharing an INF file should use the same screening settings.
#'
#' @export
hilandyn2_composite <- function(
  FORCE_datacube,
  satellite_family = c("Sentinel-2", "Landsat"),
  tiles = NULL,
  aoi = NULL,
  start_date = NULL,
  end_date = NULL,
  qai_bits = c(0L, 1L, 2L, 3L, 4L, 8L, 9L),
  method = c("medoid", "geomedian"),
  use_weights = FALSE,
  cloud_distance = 1500,
  bin_method = c("temporal", "hybrid"),
  n_bins = 1L,
  hybrid_max_days = 30,
  out_file = NULL,
  overwrite = FALSE,
  gdal_options = c("COMPRESS=ZSTD", "BIGTIFF=YES"),
  num_cores = 1L,
  num_threads = 1L
) {
  method <- match.arg(method)
  bin_method <- match.arg(bin_method)
  satellite_family <- match.arg(satellite_family)
  qai_mask <- .force_qai_mask(qai_bits)
  if (!is.logical(use_weights) || length(use_weights) != 1L ||
    is.na(use_weights)) {
    stop("use_weights must be TRUE or FALSE", call. = FALSE)
  }
  num_cores <- as.integer(num_cores)
  if (length(num_cores) != 1L || is.na(num_cores) ||
    num_cores < 1L) {
    stop(
      "num_cores must be one positive integer",
      call. = FALSE
    )
  }

  tile_paths <- .force_select_tiles(
    FORCE_datacube,
    tiles,
    aoi
  )
  tile_names <- names(tile_paths)
  output_root <- if (is.null(out_file))
    NULL else {
    .force_prepare_output_datacube(FORCE_datacube, out_file)
  }
  if (num_cores > 1L && length(tile_paths) > 1L && is.null(out_file)) {
    stop(
      "out_file must be an output data-cube root for parallel tile processing",
      call. = FALSE
    )
  }

  arguments <- list(
    satellite_family = satellite_family,
    use_weights = use_weights,
    method = method,
    n_bins = n_bins,
    bin_method = bin_method,
    start_date = start_date,
    end_date = end_date,
    cloud_distance = cloud_distance,
    hybrid_max_days = hybrid_max_days,
    overwrite = overwrite,
    num_threads = num_threads,
    gdal_options = gdal_options,
    qai_mask = qai_mask
  )
  tasks <- .force_composite_tasks(tile_paths, output_root)

  if (num_cores > 1L && length(tasks) > 1L) {
    previous_plan <- future::plan()
    on.exit(
      future::plan(previous_plan),
      add = TRUE
    )
    future::plan(future::multisession, workers = min(num_cores, length(tasks)))
    responses <- future.apply::future_lapply(
      tasks,
      function(task, arguments) {
        worker <- get(
          ".force_composite_tile_worker",
          envir = asNamespace("hilandyn2"),
          inherits = FALSE
        )
        do.call(worker, c(task, list(arguments = arguments)))
      },
      arguments = arguments,
      future.packages = "hilandyn2",
      future.globals = FALSE,
      future.seed = NULL
    )
  } else {
    responses <- lapply(
      tasks, function(task) {
        do.call(
          .force_composite_tile_worker,
          c(task, list(arguments = arguments))
        )
      }
    )
  }

  for (i in seq_along(responses)) {
    if (length(responses[[i]]$warnings)) {
      for (warning_message in responses[[i]]$warnings) {
        warning(
          "[",
          tile_names[i],
          "] ",
          warning_message,
          call. = FALSE
        )
      }
    }
  }
  results <- lapply(
    responses, function(response) {
      lapply(
        response$value, function(product) {
          if (!is.character(product))
          return(product)
          rasters <- lapply(product, terra::rast)
          if (length(rasters) == 1L)
          rasters[[1L]] else do.call(c, rasters)
        }
      )
    }
  )
  names(results) <- tile_names
  if (length(results) == 1L)
    return(invisible(results[[1L]]))
  invisible(results)
}

.force_qai_mask <- function(qai_bits) {
  if (!is.numeric(qai_bits) || is.complex(qai_bits) ||
    !is.null(dim(qai_bits)) || any(!is.finite(qai_bits)) ||
    any(qai_bits < 0 | qai_bits > 15 | qai_bits != floor(qai_bits))) {
    stop(
      "qai_bits must be a vector of integer bit positions from 0 to 15",
      call. = FALSE
    )
  }
  as.integer(sum(2^unique(qai_bits)))
}

.hilandyn2_composite_tile <- function(
  FORCE_datacube,
  satellite_family,
  use_weights = FALSE,
  method = c("medoid", "geomedian"),
  n_bins = 1L,
  bin_method = c("temporal", "hybrid"),
  start_date = NULL,
  end_date = NULL,
  cloud_distance = 1500,
  hybrid_max_days = 15,
  out_file = NULL,
  overwrite = FALSE,
  num_threads = 1L,
  gdal_options = c("COMPRESS=ZSTD", "BIGTIFF=YES"),
  qai_mask = 799L
) {

  method <- match.arg(method)
  bin_method <- match.arg(bin_method)
  n_bins <- as.integer(n_bins)
  num_threads <- as.integer(num_threads)
  if (!is.logical(use_weights) ||
    length(use_weights) != 1L || is.na(use_weights)) {
    stop(
      "use_weights must be TRUE or FALSE",
      call. = FALSE
    )
  }
  if (length(n_bins) != 1L || is.na(n_bins) ||
    n_bins < 1L) {
    stop(
      "n_bins must be one positive integer",
      call. = FALSE
    )
  }
  if (length(num_threads) != 1L || is.na(num_threads) ||
    num_threads < 1L) {
    stop(
      "num_threads must be one positive integer",
      call. = FALSE
    )
  }
  if (!is.numeric(cloud_distance) ||
    length(cloud_distance) != 1L || !is.finite(cloud_distance) ||
    cloud_distance <= 0) {
    stop(
      "cloud_distance must be one finite positive number",
      call. = FALSE
    )
  }
  if (!is.numeric(hybrid_max_days) ||
    length(hybrid_max_days) != 1L || !is.finite(hybrid_max_days) ||
    hybrid_max_days < 0) {
    stop(
      "hybrid_max_days must be one finite non-negative number",
      call. = FALSE
    )
  }

  datacube <- .force_datacube_input(
    FORCE_datacube,
    use_weights,
    satellite_family,
    start_date = start_date,
    end_date = end_date
  )
  FORCE_BOA <- datacube$boa
  FORCE_QAI <- datacube$qai
  FORCE_DST_data <- datacube$dst
  dates <- datacube$dates
  acquisition_sensors <- datacube$sensor_ids
  acquisitions <- length(dates)
  sensor <- datacube$sensor
  band_names <- sensor$band_names
  spectral_bands <- length(band_names)

  if (terra::nlyr(FORCE_QAI) != acquisitions) {
    stop(
      "FORCE QAI data must have one layer per acquisition",
      call. = FALSE
    )
  }
  if (use_weights && (is.null(FORCE_DST_data) ||
    terra::nlyr(FORCE_DST_data) != acquisitions)) {
    stop(
      "FORCE DST data must have one layer per acquisition",
      call. = FALSE
    )
  }
  matching_geometry <- terra::compareGeom(
    FORCE_BOA,
    FORCE_QAI,
    stopOnError = FALSE
  )
  if (use_weights) {
    matching_geometry <- matching_geometry &&
      terra::compareGeom(
        FORCE_BOA,
        FORCE_DST_data,
        stopOnError = FALSE
      )
  }
  if (!matching_geometry) {
    stop(
      "all input rasters must have matching geometry",
      call. = FALSE
    )
  }

  order_dates <- order(dates)
  dates <- dates[order_dates]
  acquisition_sensors <- acquisition_sensors[order_dates]
  sr_groups <- matrix(
    seq_len(terra::nlyr(FORCE_BOA)),
    ncol = spectral_bands,
    byrow = TRUE
  )
  input <- do.call(
    c, lapply(
      order_dates, function(i) {
        acquisition <- list(
          FORCE_BOA[[sr_groups[i,
          ]]]
        )
        if (use_weights)
          acquisition <- c(
          acquisition,
          list(FORCE_DST_data[[i]])
        )
        acquisition <- c(
          acquisition,
          list(FORCE_QAI[[i]])
        )
        do.call(c, acquisition)
      }
    )
  )
  input_sources <- terra::sources(input)

  start_date <- if (is.null(start_date))
    min(dates) else as.Date(start_date)
  end_date <- if (is.null(end_date))
    max(dates) else as.Date(end_date)
  if (is.na(start_date) ||
    is.na(end_date) ||
    end_date < start_date) {
    stop(
      "start_date and end_date must define a valid interval",
      call. = FALSE
    )
  }
  total_days <- as.integer(end_date - start_date) +
    1L
  if (n_bins > total_days) {
    stop(
      "n_bins cannot exceed the number of days in the compositing interval",
      call. = FALSE
    )
  }
  edges <- floor(
    seq(
      0,
      total_days,
      length.out = n_bins +
        1L
    )
  )
  bin_starts <- start_date +
    edges[-length(edges)]
  bin_ends <- start_date + edges[-1L] -
    1L
  bin_ends[n_bins] <- end_date
  intervals <- cbind(
    as.integer(bin_starts - as.Date("1970-01-01")),
    as.integer(bin_ends - as.Date("1970-01-01"))
  )
  acquisition_days <- as.integer(dates - as.Date("1970-01-01"))
  acquisition_doys <- as.integer(format(dates, "%j"))
  acquisition_years <- as.integer(format(dates, "%Y"))
  target_dates <- bin_starts +
    floor(as.integer(bin_ends - bin_starts)/2)
  target_days <- as.integer(target_dates - as.Date("1970-01-01"))
  target_doys <- as.integer(format(target_dates, "%j"))
  output_files <- .force_composite_product_files(
    out_file,
    basename(FORCE_datacube),
    target_dates,
    sensor$target_sensor,
    method
  )
  persistent <- !is.null(out_file)
  if (persistent && !overwrite) {
    .force_remove_zero_byte_files(c(
      output_files$reflectance,
      output_files$information
    ))
  }
  write_reflectance <- rep(
    !persistent || overwrite,
    n_bins
  ) |
    !file.exists(output_files$reflectance)
  write_information <- rep(
    !persistent || overwrite,
    n_bins
  ) |
    !file.exists(output_files$information)
  if (!any(
    c(
      write_reflectance,
      write_information
    )
  )) {
    return(invisible(output_files))
  }

  information_names <- c(
    "QAI",
    "NOBS",
    "DOY",
    "YEAR",
    "D_TDOY",
    "SENSOR"
  )
  method_id <- match(method, c("medoid", "geomedian")) -
    1L
  hybrid_binning <- identical(bin_method, "hybrid")
  products <- .force_write_composites(
    input = input,
    input_sources = input_sources,
    reflectance_files = output_files$reflectance,
    information_files = output_files$information,
    write_reflectance = write_reflectance,
    write_information = write_information,
    band_names = band_names,
    information_names = information_names,
    spectral_bands = spectral_bands,
    acquisition_days = acquisition_days,
    acquisition_doys = acquisition_doys,
    acquisition_years = acquisition_years,
    acquisition_sensors = acquisition_sensors,
    intervals = intervals,
    target_days = target_days,
    target_doys = target_doys,
    bins = n_bins,
    method = method_id,
    hybrid_binning = hybrid_binning,
    use_weights = use_weights,
    requested_distance = cloud_distance,
    hybrid_max_days = hybrid_max_days,
    num_threads = num_threads,
    overwrite = overwrite,
    gdal_options = gdal_options,
    qai_mask = qai_mask
  )
  if (persistent)
    return(invisible(output_files))
  stack_products <- function(products) {
    if (length(products) == 1L)
      products[[1L]] else do.call(c, products)
  }
  invisible(
    list(
      reflectance = stack_products(products$reflectance),
      information = stack_products(products$information)
    )
  )
}

.force_write_composites <- function(
  input,
  input_sources,
  reflectance_files,
  information_files,
  write_reflectance,
  write_information,
  band_names,
  information_names,
  spectral_bands,
  acquisition_days,
  acquisition_doys,
  acquisition_years,
  acquisition_sensors,
  intervals,
  target_days,
  target_doys,
  bins,
  method,
  hybrid_binning,
  use_weights,
  requested_distance,
  hybrid_max_days,
  num_threads,
  overwrite,
  gdal_options,
  qai_mask
) {

  terra_tempdir <- terra::terraOptions(print = FALSE)$tempdir
  if (!length(terra_tempdir) ||
    !dir.exists(terra_tempdir)) {
    terra_tempdir <- tempdir()
  }
  write_copies <- 2L
  output_layers <- spectral_bands +
    length(information_names)
  output_template <- terra::rast(
    input, nlyrs = output_layers *
      bins
  )
  sizing_template <- c(input, output_template)
  blocks <- terra::blocks(
    sizing_template,
    n = write_copies
  )
  rm(
    sizing_template,
    output_template
  )

  reflectance <- vector("list", bins)
  information <- vector("list", bins)
  reflectance_open <- rep(FALSE, bins)
  information_open <- rep(FALSE, bins)
  input_open <- FALSE
  on.exit(
    {
      for (bin in which(reflectance_open)) {
        try(
          terra::writeStop(reflectance[[bin]]),
          silent = TRUE
        )
      }
      for (bin in which(information_open)) {
        try(
          terra::writeStop(information[[bin]]),
          silent = TRUE
        )
      }
      if (input_open) try(
        terra::readStop(input),
        silent = TRUE
      )
    }, add = TRUE
  )

  for (bin in seq_len(bins)) {
    if (write_reflectance[bin]) {
      reflectance_file <- reflectance_files[bin]
      if (!nzchar(reflectance_file)) {
        reflectance_file <- tempfile(
          "hilandyn2_composite_",
          tmpdir = terra_tempdir,
          fileext = ".tif"
        )
      }
      reflectance[[bin]] <- terra::rast(
        input, nlyrs = spectral_bands
      )
      names(reflectance[[bin]]) <- band_names
      terra::writeStart(
        reflectance[[bin]],
        filename = reflectance_file,
        overwrite = overwrite,
        n = write_copies,
        sources = input_sources,
        wopt = list(
          filetype = "GTiff",
          datatype = "FLT4S",
          NAflag = -9999,
          gdal = gdal_options
        )
      )
      reflectance_open[bin] <- TRUE
    }
    if (write_information[bin]) {
      information_file <- information_files[bin]
      if (!nzchar(information_file)) {
        information_file <- tempfile(
          "hilandyn2_information_",
          tmpdir = terra_tempdir,
          fileext = ".tif"
        )
      }
      information[[bin]] <- terra::rast(
        input,
        nlyrs = length(information_names)
      )
      names(information[[bin]]) <- information_names
      terra::writeStart(
        information[[bin]],
        filename = information_file,
        overwrite = overwrite,
        n = write_copies,
        sources = input_sources,
        wopt = list(
          filetype = "GTiff",
          datatype = "INT2S",
          NAflag = -9999,
          gdal = gdal_options
        )
      )
      information_open[bin] <- TRUE
    }
  }

  terra::readStart(input)
  input_open <- TRUE

  for (chunk in seq_len(blocks$n)) {
    rows <- blocks$nrows[chunk]
    values <- input@pntr$readValues(
      blocks$row[chunk] -
        1L,
      rows,
      0L,
      terra::ncol(input)
    )
    result <- force_composite_chunk_cpp(
      values,
      rows *
        terra::ncol(input),
      terra::nlyr(input),
      spectral_bands,
      acquisition_days,
      acquisition_doys,
      acquisition_years,
      acquisition_sensors,
      intervals,
      target_days,
      target_doys,
      bins,
      method,
      hybrid_binning,
      use_weights,
      requested_distance,
      hybrid_max_days,
      num_threads,
      qai_mask
    )
    for (bin in seq_len(bins)) {
      offset <- (bin -
        1L) * output_layers
      if (write_reflectance[bin]) {
        terra::writeValues(
          reflectance[[bin]],
          result[,
          offset +
            seq_len(spectral_bands),
          drop = FALSE],
          blocks$row[chunk],
          rows
        )
      }
      if (write_information[bin]) {
        terra::writeValues(
          information[[bin]],
          result[,
          offset +
            spectral_bands +
            seq_along(information_names),
          drop = FALSE],
          blocks$row[chunk],
          rows
        )
      }
    }
    rm(values, result)
    gc(verbose = FALSE)
  }

  terra::readStop(input)
  input_open <- FALSE
  for (bin in seq_len(bins)) {
    if (write_reflectance[bin]) {
      reflectance[[bin]] <- terra::writeStop(reflectance[[bin]])
      reflectance_open[bin] <- FALSE
    }
    if (write_information[bin]) {
      information[[bin]] <- terra::writeStop(information[[bin]])
      information_open[bin] <- FALSE
    }
  }
  list(
    reflectance = reflectance,
    information = information
  )
}

.force_composite_tasks <- function(tile_paths, output_root) {
  lapply(
    seq_along(tile_paths),
    function(i) {
      list(
        tile_path = unname(tile_paths[i]),
        tile_name = names(tile_paths)[i],
        out_file = output_root
      )
    }
  )
}

.force_composite_tile_worker <- function(
  tile_path,
  tile_name,
  out_file,
  arguments
) {
  warning_messages <- character()
  tile_function <- get(
    ".hilandyn2_composite_tile",
    envir = asNamespace("hilandyn2"),
    inherits = FALSE
  )
  result <- withCallingHandlers(
    do.call(
      tile_function, c(
        list(
          FORCE_datacube = tile_path,
          out_file = out_file
        ),
        arguments
      )
    ),
    warning = function(condition) {
      warning_messages <<- c(warning_messages, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(
    value = result,
    warnings = warning_messages,
    tile = tile_name
  )
}
