.s2_band_names <- c(
  "BLUE",
  "GREEN",
  "RED",
  "REDEDGE1",
  "REDEDGE2",
  "REDEDGE3",
  "NIR1",
  "NIR2",
  "SWIR1",
  "SWIR2"
)

.landsat_band_names <- c(
  "BLUE",
  "GREEN",
  "RED",
  "NIR",
  "SWIR1",
  "SWIR2"
)

.force_index_names <- c(
  "NDVI",
  "NDRE1",
  "NDRE2",
  "NDRE3",
  "NDMI",
  "NBR",
  "MSI",
  "TCB",
  "TCG",
  "TCW",
  "TCA",
  "CRSWIR",
  "CRRE",
  "DRS",
  "IRECI",
  "NMDI"
)
.s2_index_names <- .force_index_names

.landsat_index_names <- c(
  "NDVI",
  "NDMI",
  "NBR",
  "MSI",
  "TCB",
  "TCG",
  "TCW",
  "TCA",
  "DRS",
  "NMDI"
)

.force_sensor_s2 <- 0L
.force_sensor_landsat <- 1L

.spectral_change_directions <- c(
  BLUE = 1,
  GREEN = 1,
  RED = 1,
  REDEDGE1 = 1,
  REDEDGE2 = 1,
  REDEDGE3 = 1,
  NIR = -1,
  NIR1 = -1,
  NIR2 = -1,
  SWIR1 = 1,
  SWIR2 = 1,
  NDVI = -1,
  NDRE1 = -1,
  NDRE2 = -1,
  NDRE3 = -1,
  NDMI = -1,
  NBR = -1,
  MSI = 1,
  TCB = 1,
  TCG = -1,
  TCW = -1,
  TCA = -1,
  CRSWIR = 1,
  CRRE = 1,
  DRS = 1,
  IRECI = -1,
  NMDI = -1
)

.normalise_spectral_names <- function(x) {
  x <- toupper(trimws(x))
  sub(
    "\\.[0-9]+$",
    "",
    x
  )
}

.match_spectral_selection <- function(
  selected,
  available,
  argument
) {
  if (is.null(selected))
    return(NULL)
  if (!is.character(selected)) {
    stop(
      argument,
      " must be a character vector",
      call. = FALSE
    )
  }
  selected <- .normalise_spectral_names(selected)
  if (anyDuplicated(selected)) {
    stop(
      argument,
      " contains duplicate names",
      call. = FALSE
    )
  }
  unknown <- setdiff(selected, available)
  if (length(unknown)) {
    stop(
      argument,
      " contains unknown names: ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  available[available %in% selected]
}

.layer_groups <- function(
  x,
  time_steps,
  argument
) {
  layers_per_step <- terra::nlyr(x)/time_steps
  if (layers_per_step%%1 != 0) {
    stop(
      argument,
      " does not contain the same number of layers at every time step",
      call. = FALSE
    )
  }
  matrix(
    seq_len(terra::nlyr(x)),
    ncol = layers_per_step,
    byrow = TRUE
  )
}

.force_family_spec <- function(satellite_family) {
  satellite_family <- match.arg(satellite_family, c("Sentinel-2", "Landsat"))
  if (satellite_family == "Sentinel-2") {
    return(
      list(
        argument = satellite_family,
        name = "Sentinel-2",
        filename_prefix = "SEN",
        target_sensor = "SEN2L",
        id = .force_sensor_s2,
        band_names = .s2_band_names
      )
    )
  }
  list(
    argument = satellite_family,
    name = "Landsat",
    filename_prefix = "LND",
    target_sensor = "LNDLG",
    id = .force_sensor_landsat,
    band_names = .landsat_band_names
  )
}

.force_product_keys <- function(files, product) {
  basenames <- basename(files)
  suffix <- paste0(
    "_",
    product,
    "\\.(tif|tiff|vrt)$"
  )
  valid <- grepl(
    paste0("^[0-9]{8}_.+", suffix),
    basenames,
    ignore.case = TRUE
  )
  if (any(!valid)) {
    stop(
      "invalid FORCE ",
      product,
      " filename",
      call. = FALSE
    )
  }
  sub(
    suffix,
    "",
    basenames,
    ignore.case = TRUE
  )
}

.force_boa_band_spec <- function(x, family) {
  layer_names <- .normalise_spectral_names(names(x))
  positions <- match(family$band_names, layer_names)
  if (!anyNA(positions)) {
    family$positions <- positions
    return(family)
  }
  if (terra::nlyr(x) != length(family$band_names)) {
    stop(
      family$name,
      " FORCE BOA data must contain ",
      length(family$band_names),
      " layers per acquisition",
      call. = FALSE
    )
  }
  family$positions <- seq_along(family$band_names)
  family
}

.force_datacube_input <- function(
  path,
  use_weights,
  satellite_family,
  start_date = NULL,
  end_date = NULL
) {
  if (!is.character(path) ||
    length(path) != 1L || is.na(path) ||
    !dir.exists(path)) {
    stop("FORCE_datacube must be one existing directory", call. = FALSE)
  }
  family <- .force_family_spec(satellite_family)
  files <- list.files(
    path,
    pattern = "\\.(tif|tiff|vrt)$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  family_pattern <- paste0(
    "(^|_)",
    family$filename_prefix,
    "[^_]*(_|$)"
  )
  files <- files[grepl(
    family_pattern,
    basename(files),
    ignore.case = TRUE
  )]
  boa_files <- sort(files[grepl(
    "_BOA\\.(tif|tiff|vrt)$",
    files,
    ignore.case = TRUE
  )])
  qai_files <- files[grepl(
    "_QAI\\.(tif|tiff|vrt)$",
    files,
    ignore.case = TRUE
  )]
  dst_files <- files[grepl(
    "_DST\\.(tif|tiff|vrt)$",
    files,
    ignore.case = TRUE
  )]
  if (!length(boa_files)) {
    stop(
      "FORCE_datacube does not contain BOA files for satellite_family = \"",
      family$argument,
      "\"",
      call. = FALSE
    )
  }
  if (!length(qai_files)) {
    stop(
      "FORCE_datacube does not contain QAI files for satellite_family = \"",
      family$argument,
      "\"",
      call. = FALSE
    )
  }

  boa_keys <- .force_product_keys(boa_files, "BOA")
  if (anyDuplicated(boa_keys)) {
    stop(
      "FORCE_datacube contains duplicate BOA acquisitions; provide one tile",
      call. = FALSE
    )
  }
  dates <- as.Date(
    substr(
      boa_keys,
      1L,
      8L
    ),
    format = "%Y%m%d"
  )
  if (anyNA(dates)) {
    stop("FORCE BOA filenames contain invalid acquisition dates", call. = FALSE)
  }
  normalize_bound <- function(value, label) {
    if (is.null(value))
      return(NULL)
    value <- as.Date(value)
    if (length(value) != 1L || is.na(value)) {
      stop(
        label,
        " must be NULL or one valid date",
        call. = FALSE
      )
    }
    value
  }
  start_date <- normalize_bound(start_date, "start_date")
  end_date <- normalize_bound(end_date, "end_date")
  if (!is.null(start_date) &&
    !is.null(end_date) &&
    end_date < start_date) {
    stop("start_date and end_date must define a valid interval", call. = FALSE)
  }
  selected <- rep(TRUE, length(dates))
  if (!is.null(start_date))
    selected <- selected & dates >= start_date
  if (!is.null(end_date))
    selected <- selected & dates <= end_date
  if (!any(selected)) {
    stop(
      "FORCE_datacube contains no BOA acquisitions in the requested interval",
      call. = FALSE
    )
  }
  boa_files <- boa_files[selected]
  boa_keys <- boa_keys[selected]
  dates <- dates[selected]

  qai_keys <- .force_product_keys(qai_files, "QAI")
  if (anyDuplicated(qai_keys)) {
    stop(
      "FORCE_datacube contains duplicate QAI acquisitions; provide one tile",
      call. = FALSE
    )
  }
  qai_match <- match(boa_keys, qai_keys)
  if (anyNA(qai_match)) {
    stop("each FORCE BOA acquisition must have a matching QAI file",
      call. = FALSE)
  }
  qai_files <- qai_files[qai_match]

  if (use_weights) {
    dst_keys <- .force_product_keys(dst_files, "DST")
    if (anyDuplicated(dst_keys)) {
      stop(
        "FORCE_datacube contains duplicate DST acquisitions; provide one tile",
        call. = FALSE
      )
    }
    dst_match <- match(boa_keys, dst_keys)
    if (anyNA(dst_match)) {
      stop(
        "use_weights = TRUE requires matching FORCE DST files; ",
        "missing for BOA: ",
        paste(boa_files[is.na(dst_match)], collapse = ", "),
        call. = FALSE
      )
    }
    dst_files <- dst_files[dst_match]
  }

  boa_sources <- lapply(boa_files, terra::rast)
  specs <- lapply(
    boa_sources,
    .force_boa_band_spec,
    family = family
  )
  sensor <- family
  boa <- do.call(
    c, lapply(
      seq_along(boa_sources),
      function(i) {
        boa_sources[[i]][[specs[[i]]$positions]]
      }
    )
  )
  names(boa) <- rep(sensor$band_names, times = length(boa_sources))
  qai <- terra::rast(qai_files)

  dst <- if (use_weights) terra::rast(dst_files) else NULL
  sensor_tokens <- toupper(sub(
    "^[0-9]{8}_[^_]+_([^_]+)$",
    "\\1",
    boa_keys,
    perl = TRUE
  ))
  sensor_levels <- if (identical(family$argument, "Sentinel-2")) {
    c(
      "SEN2A",
      "SEN2B",
      "SEN2C",
      "SEN2D"
    )
  } else {
    c(
      "LND04",
      "LND05",
      "LND07",
      "LND08",
      "LND09"
    )
  }
  sensor_ids <- match(sensor_tokens, sensor_levels)
  if (anyNA(sensor_ids)) {
    stop(
      "unsupported FORCE sensor token in BOA filename: ",
      paste(
        unique(sensor_tokens[is.na(sensor_ids)]),
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  list(
    boa = boa,
    qai = qai,
    dst = dst,
    dates = dates,
    sensor = sensor,
    sensor_tokens = sensor_tokens,
    sensor_ids = as.integer(sensor_ids)
  )
}

.force_time_series_spec <- function(
  x,
  time_steps,
  family
) {
  groups <- .layer_groups(
    x,
    time_steps,
    "FORCE reflectance data"
  )
  first_names <- .normalise_spectral_names(
    names(x)[groups[1,
      ]]
  )
  positions <- match(family$band_names, first_names)
  if (anyNA(positions)) {
    if (ncol(groups) != length(family$band_names)) {
      stop(
        "FORCE data for satellite_family = \"",
        family$argument,
        "\" must contain ",
        length(family$band_names),
        " reflectance bands per time step",
        call. = FALSE
      )
    }
    positions <- seq_along(family$band_names)
  }
  family$groups <- groups
  family$positions <- positions
  family
}

.force_analysis_spec <- function(
  satellite_family,
  spectral_bands = NULL,
  spectral_indices = NULL
) {
  family <- .force_family_spec(satellite_family)
  selected_bands <- if (is.null(spectral_bands)) {
    family$band_names
  } else {
    .match_spectral_selection(
      spectral_bands,
      family$band_names,
      "spectral_bands"
    )
  }
  available_indices <- if (family$id == .force_sensor_s2) {
    .s2_index_names
  } else {
    .landsat_index_names
  }
  selected_indices <- if (is.null(spectral_indices)) {
    character()
  } else {
    .match_spectral_selection(
      spectral_indices,
      available_indices,
      "spectral_indices"
    )
  }
  if (!length(selected_bands) &&
    !length(selected_indices)) {
    stop(
      "no spectral variables were selected",
      call. = FALSE
    )
  }

  list(
    family = family,
    band_names = selected_bands,
    band_positions = match(selected_bands, family$band_names),
    index_names = selected_indices,
    index_ids = match(selected_indices, .force_index_names) -
      1L,
    variable_names = c(selected_bands, selected_indices),
    external_variables = length(selected_bands),
    computed_variables = length(selected_indices),
    count = length(selected_bands) +
      length(selected_indices)
  )
}

.select_grouped_layers <- function(
  x,
  groups,
  positions,
  layer_names
) {
  selected <- do.call(
    c, lapply(
      seq_len(nrow(groups)),
      function(i) {
        x[[groups[i, positions]]]
      }
    )
  )
  names(selected) <- rep(layer_names, times = nrow(groups))
  selected
}

.canonical_force_raster <- function(
  x,
  time_steps,
  family
) {
  spec <- .force_time_series_spec(
    x,
    time_steps,
    family
  )
  .select_grouped_layers(
    x,
    spec$groups,
    spec$positions,
    spec$band_names
  )
}

.force_map_chunk_values <- function(
  values,
  input_rows,
  input_columns,
  time_steps,
  selected_band_positions,
  index_ids,
  sensor,
  rescale_reflectance,
  num_threads,
  mask_values = NULL
) {
  spectral_bands <- if (sensor == .force_sensor_s2) {
    length(.s2_band_names)
  } else {
    length(.landsat_band_names)
  }
  if (input_columns != spectral_bands * time_steps) {
    stop(
      "mapping chunk has an invalid number of FORCE reflectance layers",
      call. = FALSE
    )
  }

  if (is.null(mask_values))
    mask_values <- numeric()
  force_map_chunk_values_cpp(
    values,
    mask_values,
    input_rows,
    input_columns,
    time_steps,
    selected_band_positions,
    index_ids,
    sensor,
    rescale_reflectance,
    num_threads
  )
}

.infer_change_directions <- function(variable_names) {
  normalised <- .normalise_spectral_names(variable_names)
  directions <- unname(.spectral_change_directions[normalised])
  list(names = normalised, directions = directions)
}

.resolve_change_directions <- function(
  variable_names,
  cng_dir,
  external_variables,
  computed_variables
) {
  inferred <- .infer_change_directions(variable_names)
  total_variables <- length(variable_names)
  if (is.null(cng_dir)) {
    if (anyNA(inferred$directions)) {
      unknown <- unique(variable_names[is.na(inferred$directions)])
      stop(
        "cannot infer change direction for: ",
        paste(unknown, collapse = ", "),
        "; provide cng_dir for externally supplied rasters",
        call. = FALSE
      )
    }
    return(inferred$directions)
  }

  if (!is.numeric(cng_dir) ||
    anyNA(cng_dir) ||
    any(!cng_dir %in% c(-1, 1))) {
    stop(
      "cng_dir must contain only -1 and 1",
      call. = FALSE
    )
  }
  if (computed_variables > 0L && length(cng_dir) == external_variables) {
    computed_directions <- inferred$directions[external_variables +
      seq_len(computed_variables)]
    return(c(cng_dir, computed_directions))
  }
  if (length(cng_dir) != total_variables) {
    stop(
      "length of cng_dir must match either all variables or only the ",
      "externally supplied variables",
      call. = FALSE
    )
  }
  cng_dir
}
