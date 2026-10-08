# Internal tiled dispatcher used by hilandyn2_map().
.hilandyn2_map_tiled <- function(
  FORCE_datacube,
  satellite_family = c("Sentinel-2", "Landsat"),
  composite_product = NULL,
  out_path = NULL,
  temp_dir = NULL,
  ts_ids,
  tiles = NULL,
  spectral_bands = NULL,
  spectral_indices = NULL,
  include_adjacent = TRUE,
  num_cores = 1L,
  win_side = 3L,
  num_threads = 1L,
  num_copies = 16L,
  overwrite = FALSE,
  map_arguments = list()
) {

  satellite_family <- match.arg(satellite_family)
  product <- .force_map_composite_spec(
    FORCE_datacube,
    satellite_family,
    composite_product
  )
  analysis_spec <- .force_analysis_spec(
    satellite_family,
    spectral_bands,
    spectral_indices
  )
  num_cores <- as.integer(num_cores)
  num_threads <- as.integer(num_threads)
  num_copies <- as.integer(num_copies)
  win_side <- as.integer(win_side)
  if (length(num_cores) != 1L || is.na(num_cores) || num_cores < 1L) {
    stop("num_cores must be one positive integer", call. = FALSE)
  }
  if (length(num_threads) != 1L || is.na(num_threads) || num_threads < 1L) {
    stop("num_threads must be one positive integer", call. = FALSE)
  }
  if (length(num_copies) != 1L || is.na(num_copies) || num_copies < 1L) {
    stop("num_copies must be one positive integer", call. = FALSE)
  }
  if (length(win_side) != 1L || is.na(win_side) || win_side < 1L ||
      win_side %% 2L == 0L) {
    stop("win_side must be one positive odd integer", call. = FALSE)
  }
  if (!is.logical(include_adjacent) || length(include_adjacent) != 1L ||
      is.na(include_adjacent)) {
    stop("include_adjacent must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.null(out_path) &&
      (!is.character(out_path) || length(out_path) != 1L || is.na(out_path) ||
       !nzchar(out_path))) {
    stop("out_path must be NULL or one non-empty path", call. = FALSE)
  }
  if (!is.null(out_path)) {
    out_path <- normalizePath(
      out_path,
      winslash = "/",
      mustWork = FALSE
    )
    dir.create(
      out_path,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  if (!is.null(temp_dir) &&
      (!is.character(temp_dir) || length(temp_dir) != 1L ||
       is.na(temp_dir) || !nzchar(temp_dir))) {
    stop("temp_dir must be NULL or one non-empty path", call. = FALSE)
  }
  if (is.null(temp_dir)) {
    temp_root <- terra::terraOptions(print = FALSE)$tempdir
    if (!length(temp_root) || is.na(temp_root) || !dir.exists(temp_root)) {
      temp_root <- base::tempdir()
    }
  } else {
    temp_root <- normalizePath(
      temp_dir,
      winslash = "/",
      mustWork = FALSE
    )
    dir.create(
      temp_root,
      recursive = TRUE,
      showWarnings = FALSE
    )
    if (!dir.exists(temp_root)) {
      stop("temp_dir could not be created", call. = FALSE)
    }
  }

  force_tiles <- .force_tile_directories(
    FORCE_datacube,
    product$pattern,
    "FORCE_datacube"
  )
  available_names <- sort(names(force_tiles))

  if (is.null(tiles)) {
    target_names <- available_names
  } else {
    if (!is.character(tiles) || !length(tiles) || anyNA(tiles) ||
        any(!nzchar(tiles)) || anyDuplicated(toupper(tiles))) {
      stop("tiles must contain unique non-empty tile names", call. = FALSE)
    }
    matched <- match(toupper(tiles), toupper(available_names))
    if (anyNA(matched)) {
      stop(
        "tiles not found in FORCE_datacube: ",
        paste(tiles[is.na(matched)], collapse = ", "),
        call. = FALSE
      )
    }
    target_names <- available_names[matched]
  }
  if (is.null(out_path) && (length(target_names) > 1L || num_cores > 1L)) {
    stop("out_path is required for multiple or parallel tile processing",
         call. = FALSE)
  }

  force_files <- .map_tile_file_table(
    force_tiles[available_names],
    product$pattern,
    "FORCE_datacube"
  )

  tasks <- lapply(target_names, function(tile) {
    neighbors <- if (include_adjacent) {
      .map_adjacent_tiles(tile, available_names)
    } else {
      tile
    }
    list(
      tile = tile,
      neighbors = neighbors,
      force_files = force_files[neighbors],
      pad_focal_edges = !include_adjacent,
      tile_temp_dir = tempfile(
        pattern = paste0(
          "hilandyn2_",
          tile,
          "_"
        ), tmpdir = temp_root
      ),
      tile_out_path = if (is.null(out_path)) NULL else
        file.path(out_path, tile)
    )
  })
  worker_arguments <- list(
    ts_ids = ts_ids,
    satellite_family = satellite_family,
    composite_product = product$product,
    analysis_spec = analysis_spec,
    win_side = win_side,
    num_threads = num_threads,
    num_copies = num_copies,
    overwrite = overwrite,
    map_arguments = map_arguments
  )

  if (num_cores > 1L && length(tasks) > 1L) {
    previous_plan <- future::plan()
    on.exit(future::plan(previous_plan), add = TRUE)
    future::plan(future::multisession, workers = min(num_cores, length(tasks)))
    responses <- future.apply::future_lapply(
      tasks,
      function(task, worker_arguments) {
        worker <- get(
          ".map_tile_worker",
          envir = asNamespace("hilandyn2"),
          inherits = FALSE
        )
        do.call(worker, c(list(task = task), worker_arguments))
      },
      worker_arguments = worker_arguments,
      future.packages = "hilandyn2",
      future.chunk.size = 1L,
      future.seed = NULL
    )
  } else {
    responses <- lapply(tasks, function(task) {
      do.call(.map_tile_worker, c(list(task = task), worker_arguments))
    })
  }

  for (i in seq_along(responses)) {
    for (warning_message in responses[[i]]$warnings) {
      warning(
        "[",
        target_names[i],
        "] ",
        warning_message,
        call. = FALSE
      )
    }
  }
  results <- lapply(responses, function(response) {
    if (is.character(response$value)) terra::rast(response$value)
    else response$value
  })
  names(results) <- target_names
  if (length(results) == 1L) return(invisible(results[[1L]]))
  invisible(results)
}

.map_tile_file_key <- function(files) {
  sub(
    "_TILE_X-?[0-9]+_Y-?[0-9]+(?=\\.(tif|tiff|vrt)$)",
    "",
    basename(files),
    ignore.case = TRUE,
    perl = TRUE
  )
}

.map_tile_file_table <- function(
  tile_paths,
  pattern,
  argument
) {
  lapply(tile_paths, function(path) {
    files <- sort(list.files(
      path,
      pattern = pattern,
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    ))
    keys <- .map_tile_file_key(files)
    if (anyDuplicated(toupper(keys))) {
      stop(
        argument,
        " contains duplicate dataset names in tile ",
        basename(path),
        call. = FALSE
      )
    }
    stats::setNames(files, toupper(keys))
  })
}

.map_adjacent_tiles <- function(tile, available_tiles) {
  coordinates <- .force_tile_coordinates(available_tiles)
  target <- .force_tile_coordinates(tile)[1L, ]
  adjacent <- abs(coordinates[, "x"] - target["x"]) <= 1L &
    abs(coordinates[, "y"] - target["y"]) <= 1L
  available_tiles[adjacent]
}

.map_tile_vrt_input <- function(
  files_by_tile,
  target_tile,
  half_window,
  argument
) {
  target_files <- files_by_tile[[target_tile]]
  if (is.null(target_files) || !length(target_files)) {
    stop(
      argument,
      " has no files for target tile ",
      target_tile,
      call. = FALSE
    )
  }
  reference <- terra::rast(unname(target_files[1L]))[[1L]]
  reference_extent <- as.vector(terra::ext(reference))
  resolution <- terra::res(reference)
  buffered_extent <- c(
    reference_extent[1L] - half_window * resolution[1L], # xmin
    reference_extent[3L] - half_window * resolution[2L], # ymin
    reference_extent[2L] + half_window * resolution[1L], # xmax
    reference_extent[4L] + half_window * resolution[2L]  # ymax
  )
  if (any(!is.finite(c(buffered_extent, resolution))) ||
      buffered_extent[1L] >= buffered_extent[3L] ||
      buffered_extent[2L] >= buffered_extent[4L] || any(resolution <= 0)) {
    stop(
      argument,
      " has invalid raster geometry for target tile ",
      target_tile,
      ": extent (xmin, ymin, xmax, ymax) = ",
      paste(buffered_extent, collapse = ", "),
      "; resolution = ",
      paste(resolution, collapse = ", "),
      call. = FALSE
    )
  }
  target_keys <- names(target_files)
  tile_order <- c(target_tile, setdiff(names(files_by_tile), target_tile))
  files_by_tile <- files_by_tile[tile_order]

  datasets <- lapply(target_keys, function(key) {
    sources <- unname(vapply(
      files_by_tile,
      function(tile_files) {
      matched <- tile_files[match(key, names(tile_files))]
      if (length(matched) && !is.na(matched)) matched else NA_character_
    },
      character(1)
    ))
    sources <- sources[!is.na(sources)]
    target_raster <- terra::rast(target_files[[key]])
    band_indices <- seq_len(terra::nlyr(target_raster))
    band_options <- as.vector(rbind(
      rep("-b", length(band_indices)), band_indices
    ))
    terra::vrt(
      sources,
      options = c(
        "-te",
        buffered_extent,
        band_options,
        "-vrtnodata",
        "-9999"
      ),
      overwrite = TRUE
    )
  })
  do.call(c, datasets)
}

.map_output_filename <- function(
  out_path,
  bands,
  window_cells,
  ts_ids
) {
  file.path(out_path, paste0(
    "HILANDYN2_",
    bands,
    "_BND_",
    window_cells,
    "_CEL_",
    ts_ids[1L],
    "_",
    ts_ids[length(ts_ids)],
    ".tif"
  ))
}

.map_output_is_complete <- function(path) {
  if (!file.exists(path)) return(FALSE)
  information <- file.info(path)
  if (is.na(information$size) || information$size <= 0) return(FALSE)
  output <- tryCatch(
    suppressWarnings(terra::rast(path)),
    error = function(condition) NULL
  )
  !is.null(output) && terra::nrow(output) > 0L && terra::ncol(output) > 0L &&
    terra::nlyr(output) > 0L && terra::hasValues(output)
}

.map_tile_worker <- function(
  task,
  ts_ids,
  satellite_family,
  composite_product,
  analysis_spec,
  win_side,
  num_threads,
  num_copies,
  overwrite,
  map_arguments
) {
  warning_messages <- character()
  half_window <- (win_side - 1L) %/% 2L
  analysis_bands <- analysis_spec$count
  expected_output <- if (is.null(task$tile_out_path)) {
    NULL
  } else {
    .map_output_filename(
      task$tile_out_path,
      analysis_bands,
      win_side * win_side,
      ts_ids
    )
  }
  if (!overwrite && !is.null(expected_output) && file.exists(expected_output)) {
    if (.map_output_is_complete(expected_output)) {
      return(list(
        value = expected_output,
        warnings = warning_messages,
        tile = task$tile,
        skipped = TRUE
      ))
    }
    stop(
      "existing output is unreadable or empty: ",
      expected_output,
      ". Use overwrite = TRUE to replace it.",
      call. = FALSE
    )
  }

  previous_temp_dir <- terra::terraOptions(print = FALSE)$tempdir
  tile_temp_dir <- task$tile_temp_dir
  dir.create(
    tile_temp_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  if (!dir.exists(tile_temp_dir)) {
    stop(
      "tile temporary directory could not be created: ",
      tile_temp_dir,
      call. = FALSE
    )
  }
  cleanup_tile_temp <- FALSE
  terra::terraOptions(tempdir = tile_temp_dir)
  on.exit({
    if (length(previous_temp_dir) && !is.na(previous_temp_dir)) {
      terra::terraOptions(tempdir = previous_temp_dir)
    }
    if (cleanup_tile_temp && dir.exists(tile_temp_dir)) {
      unlink(
        tile_temp_dir,
        recursive = TRUE,
        force = TRUE
      )
    }
  }, add = TRUE)

  force_input <- if (isTRUE(task$pad_focal_edges)) {
    do.call(
      c,
      lapply(unname(task$force_files[[task$tile]]), terra::rast)
    )
  } else {
    .map_tile_vrt_input(
      task$force_files,
      task$tile,
      half_window,
      "FORCE_datacube"
    )
  }
  target_source <- unname(task$force_files[[task$tile]][1L])
  target <- terra::rast(target_source)[[1L]]
  expected_geometry <- target
  if (!is.null(map_arguments$roi_vec)) {
    roi <- terra::vect(map_arguments$roi_vec)
    if (terra::crs(roi) != terra::crs(target)) {
      roi <- terra::project(roi, terra::crs(target))
    }
    expected_geometry <- terra::crop(target, roi)
  }
  result <- withCallingHandlers(
    do.call(
      .hilandyn2_map_tile,
      c(
        list(
          FORCE_datacube = force_input,
          satellite_family = satellite_family,
          composite_product = composite_product,
          analysis_spec = analysis_spec,
          out_path = task$tile_out_path,
          ts_ids = ts_ids,
          overwrite = overwrite,
          pad_focal_edges = isTRUE(task$pad_focal_edges),
          win_side = win_side,
          num_threads = num_threads,
          num_copies = num_copies
        ),
        map_arguments
      )
    ),
    warning = function(condition) {
      warning_messages <<- c(warning_messages, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  if (!terra::compareGeom(
    result,
    expected_geometry,
    crs = TRUE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE,
    stopOnError = FALSE
  )) {
    stop(
      "mapped output does not match the expected target/ROI geometry",
      call. = FALSE
    )
  }
  if (!is.null(task$tile_out_path)) {
    rm(
      force_input,
      target,
      expected_geometry
    )
    gc(FALSE)
  }
  cleanup_tile_temp <- !is.null(task$tile_out_path)
  list(
    value = if (is.null(task$tile_out_path)) result else
      unname(terra::sources(result)[1L]),
    warnings = warning_messages,
    tile = task$tile,
    skipped = FALSE
  )
}
