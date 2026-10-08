.force_tile_pattern <- "^X-?[0-9]+_Y-?[0-9]+$"

.force_tile_directories <- function(
  path,
  pattern = "_BOA\\.(tif|tiff|vrt)$",
  argument = "FORCE_datacube",
  allow_root = TRUE
) {
  if (!is.character(path) ||
    length(path) != 1L || is.na(path) ||
    !dir.exists(path)) {
    stop(
      argument,
      " must be one existing directory",
      call. = FALSE
    )
  }
  if (!is.character(pattern) ||
    length(pattern) != 1L || is.na(pattern) ||
    !nzchar(pattern)) {
    stop(
      argument,
      " pattern must be one non-empty regular expression",
      call. = FALSE
    )
  }
  path <- normalizePath(
    path,
    winslash = "/",
    mustWork = TRUE
  )
  directories <- list.dirs(
    path,
    full.names = TRUE,
    recursive = FALSE
  )
  tile_names <- basename(directories)
  is_tile <- grepl(
    .force_tile_pattern,
    tile_names,
    ignore.case = TRUE
  )
  directories <- directories[is_tile]
  tile_names <- tile_names[is_tile]
  has_files <- vapply(
    directories,
    function(x) {
      length(
        list.files(
          x,
          pattern = pattern,
          recursive = TRUE,
          ignore.case = TRUE
        )
      ) >
        0L
    },
    logical(1)
  )
  directories <- directories[has_files]
  tile_names <- tile_names[has_files]

  if (!length(directories) &&
    isTRUE(allow_root)) {
    root_has_files <- length(
      list.files(
        path,
        pattern = pattern,
        recursive = FALSE,
        ignore.case = TRUE
      )
    ) >
      0L
    if (root_has_files) {
      directories <- path
      tile_names <- basename(path)
    }
  }
  if (!length(directories)) {
    stop(
      argument,
      " does not contain matching FORCE tile rasters",
      call. = FALSE
    )
  }

  ord <- order(toupper(tile_names))
  stats::setNames(directories[ord], toupper(tile_names[ord]))
}

.force_datacube_definition <- function(path) {
  definition_file <- file.path(path, "datacube-definition.prj")
  if (!file.exists(definition_file)) {
    stop(
      "AOI tile selection requires datacube-definition.prj in FORCE_datacube",
      call. = FALSE
    )
  }
  lines <- trimws(readLines(definition_file, warn = FALSE))
  read_value <- function(field) {
    found <- grepl(
      paste0(
        "^",
        field,
        "[[:space:]]*="
      ),
      lines
    )
    if (sum(found) != 1L) {
      stop(
        "invalid datacube-definition.prj: missing ",
        field,
        call. = FALSE
      )
    }
    trimws(sub(
      "^[^=]+=[[:space:]]*",
      "",
      lines[found]
    ))
  }
  numeric_value <- function(field) {
    value <- suppressWarnings(as.numeric(read_value(field)))
    if (is.na(value) ||
      !is.finite(value)) {
      stop(
        "invalid numeric value for ",
        field,
        " in datacube-definition.prj",
        call. = FALSE
      )
    }
    value
  }
  tile_size_x <- numeric_value("TILE_SIZE_X")
  tile_size_y <- numeric_value("TILE_SIZE_Y")
  if (tile_size_x <= 0 || tile_size_y <= 0) {
    stop("FORCE tile sizes must be positive", call. = FALSE)
  }
  list(
    projection = read_value("PROJECTION"),
    origin_x = numeric_value("ORIGIN_MAP_X"),
    origin_y = numeric_value("ORIGIN_MAP_Y"),
    tile_size_x = tile_size_x,
    tile_size_y = tile_size_y
  )
}

.force_tile_coordinates <- function(tile_names) {
  matches <- regexec(
    "^X(-?[0-9]+)_Y(-?[0-9]+)$",
    toupper(tile_names),
    perl = TRUE
  )
  parts <- regmatches(
    toupper(tile_names),
    matches
  )
  valid <- lengths(parts) ==
    3L
  if (any(!valid)) {
    stop(
      "AOI selection requires FORCE tile directories named X####_Y####",
      call. = FALSE
    )
  }
  cbind(
    x = as.integer(
      vapply(
        parts,
        `[[`,
        character(1),
        2L
      )
    ),
    y = as.integer(
      vapply(
        parts,
        `[[`,
        character(1),
        3L
      )
    )
  )
}

.force_tile_grid <- function(tile_names, definition) {
  coordinates <- .force_tile_coordinates(tile_names)
  tile_x <- coordinates[, "x"]
  tile_y <- coordinates[, "y"]
  xmin <- definition$origin_x + tile_x * definition$tile_size_x
  xmax <- xmin + definition$tile_size_x
  ymax <- definition$origin_y - tile_y * definition$tile_size_y
  ymin <- ymax - definition$tile_size_y
  polygons <- vapply(
    seq_along(tile_names),
    function(i) {
      sprintf(
        "POLYGON ((%.15g %.15g, %.15g %.15g, %.15g %.15g, %.15g %.15g, %.15g %.15g))",
        xmin[i],
        ymin[i],
        xmax[i],
        ymin[i],
        xmax[i],
        ymax[i],
        xmin[i],
        ymax[i],
        xmin[i],
        ymin[i]
      )
    },
    character(1)
  )
  grid <- terra::vect(polygons, crs = definition$projection)
  grid$tile <- tile_names
  grid
}

.force_read_aoi <- function(aoi) {
  if (inherits(aoi, "SpatVector"))
    return(aoi)
  if (!is.character(aoi) ||
    length(aoi) != 1L || is.na(aoi) ||
    !file.exists(aoi)) {
    stop(
      "aoi must be one existing vector filename or a SpatVector",
      call. = FALSE
    )
  }
  terra::vect(aoi)
}

.force_select_tiles <- function(
  path,
  tiles = NULL,
  aoi = NULL
) {
  available <- .force_tile_directories(path)

  if (!is.null(tiles)) {
    if (!is.character(tiles) ||
      !length(tiles) ||
      anyNA(tiles) ||
      any(!nzchar(tiles)) ||
      anyDuplicated(toupper(tiles))) {
      stop(
        "tiles must be a non-empty character vector of unique tile names",
        call. = FALSE
      )
    }
    tile_match <- match(
      toupper(tiles),
      toupper(names(available))
    )
    if (anyNA(tile_match)) {
      stop(
        "tiles not found in FORCE_datacube: ",
        paste(
          tiles[is.na(tile_match)],
          collapse = ", "
        ),
        call. = FALSE
      )
    }
    available <- available[tile_match]
  }

  if (!is.null(aoi)) {
    definition <- .force_datacube_definition(path)
    grid <- .force_tile_grid(
      names(available),
      definition
    )
    aoi <- .force_read_aoi(aoi)
    if (!nzchar(terra::crs(aoi))) {
      stop("aoi must have a coordinate reference system", call. = FALSE)
    }
    aoi <- terra::project(aoi, definition$projection)
    intersects <- terra::relate(
      grid,
      aoi,
      relation = "intersects"
    )
    available <- available[rowSums(intersects) >
      0L]
  }

  if (!length(available)) {
    stop(
      "tile selection does not intersect any available FORCE tile",
      call. = FALSE
    )
  }
  available
}

.force_definition_source <- function(path) {
  path <- normalizePath(
    path,
    winslash = "/",
    mustWork = TRUE
  )
  candidates <- file.path(
    c(path, dirname(path)),
    "datacube-definition.prj"
  )
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) {
    stop(
      "FORCE output requires datacube-definition.prj in FORCE_datacube ",
      "or its parent directory",
      call. = FALSE
    )
  }
  normalizePath(
    candidates[1L],
    winslash = "/",
    mustWork = TRUE
  )
}

.force_prepare_output_datacube <- function(input_root, output_root) {
  if (!is.character(output_root) ||
    length(output_root) != 1L || is.na(output_root) ||
    !nzchar(output_root)) {
    stop(
      "out_file must be NULL or one non-empty output directory",
      call. = FALSE
    )
  }
  if (grepl(
    "\\.(tif|tiff|vrt)$",
    output_root,
    ignore.case = TRUE
  )) {
    stop(
      "out_file must be a FORCE data-cube root, not a raster filename",
      call. = FALSE
    )
  }
  input_root <- normalizePath(
    input_root,
    winslash = "/",
    mustWork = TRUE
  )
  output_root <- normalizePath(
    output_root,
    winslash = "/",
    mustWork = FALSE
  )
  if (identical(
    tolower(input_root),
    tolower(output_root)
  )) {
    stop(
      "out_file must be separate from the Level-2 FORCE_datacube",
      call. = FALSE
    )
  }

  source_definition <- .force_definition_source(input_root)
  dir.create(
    output_root,
    recursive = TRUE,
    showWarnings = FALSE
  )
  target_definition <- file.path(output_root, "datacube-definition.prj")
  if (file.exists(target_definition)) {
    source_lines <- readLines(source_definition, warn = FALSE)
    target_lines <- readLines(target_definition, warn = FALSE)
    if (!identical(source_lines, target_lines)) {
      stop(
        "out_file contains a different datacube-definition.prj",
        call. = FALSE
      )
    }
  } else if (!file.copy(source_definition, target_definition)) {
    stop("could not copy datacube-definition.prj to out_file", call. = FALSE)
  }

  output_root
}

.force_composite_product_files <- function(
  output_root,
  tile_name,
  target_dates,
  target_sensor,
  method
) {
  if (is.null(output_root)) {
    return(
      list(
        reflectance = rep("", length(target_dates)),
        information = rep("", length(target_dates))
      )
    )
  }
  if (!grepl(
    .force_tile_pattern,
    tile_name,
    ignore.case = TRUE
  )) {
    stop(
      "FORCE output requires tile directories named X####_Y####",
      call. = FALSE
    )
  }
  tile_dir <- file.path(output_root, tile_name)
  dir.create(
    tile_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  prefix <- paste0(
    format(target_dates, "%Y%m%d"),
    "_LEVEL3_",
    target_sensor,
    "_"
  )
  product <- if (identical(method, "medoid"))
    "MED" else "GEO"
  list(
    reflectance = file.path(tile_dir, paste0(
      prefix,
      product,
      ".tif"
    )),
    information = file.path(tile_dir, paste0(prefix, "INF.tif"))
  )
}

.force_remove_zero_byte_files <- function(paths) {
  existing <- file.exists(paths)
  if (!any(existing)) return(invisible(character()))
  existing_paths <- paths[existing]
  information <- file.info(existing_paths, extra_cols = FALSE)
  empty <- existing_paths[
    !is.na(information$size) &
      !is.na(information$isdir) &
      !information$isdir &
      information$size == 0
  ]
  if (!length(empty)) return(invisible(character()))
  unlink(empty)
  remaining <- empty[file.exists(empty)]
  if (length(remaining)) {
    stop(
      "could not remove zero-byte output file(s): ",
      paste(remaining, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(empty)
}

.force_map_composite_spec <- function(
  path,
  satellite_family,
  composite_product = NULL
) {
  if (!is.character(path) ||
    length(path) != 1L || is.na(path) ||
    !dir.exists(path)) {
    stop(
      "FORCE_datacube must be one existing directory",
      call. = FALSE
    )
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
  products <- toupper(
    sub(
      ".*_(BAP|MED|GEO)\\.(tif|tiff|vrt)$",
      "\\1",
      basename(files),
      ignore.case = TRUE,
      perl = TRUE
    )
  )
  valid <- products %in% c(
    "BAP",
    "MED",
    "GEO"
  )
  files <- files[valid]
  products <- products[valid]
  available <- unique(products)
  if (!length(available)) {
    stop(
      "FORCE_datacube does not contain BAP, MED, or GEO composite files for ",
      "satellite_family = \"",
      family$argument,
      "\"",
      call. = FALSE
    )
  }

  if (is.null(composite_product)) {
    if (length(available) > 1L) {
      stop(
        "FORCE_datacube contains multiple composite products (",
        paste(
          sort(available),
          collapse = ", "
        ),
        "); select one with composite_product",
        call. = FALSE
      )
    }
    composite_product <- available
  } else {
    if (!is.character(composite_product) ||
      length(composite_product) != 1L || is.na(composite_product) ||
      !nzchar(composite_product)) {
      stop(
        "composite_product must be NULL or one of BAP, MED, and GEO",
        call. = FALSE
      )
    }
    composite_product <- toupper(composite_product)
    if (!composite_product %in%
      c(
        "BAP",
        "MED",
        "GEO"
      )) {
      stop(
        "composite_product must be NULL or one of BAP, MED, and GEO",
        call. = FALSE
      )
    }
    if (!composite_product %in%
      available) {
      stop(
        "FORCE_datacube does not contain composite_product = \"",
        composite_product,
        "\" for satellite_family = \"",
        family$argument,
        "\"",
        call. = FALSE
      )
    }
  }

  selected <- sort(files[products == composite_product])
  list(
    product = composite_product,
    family = family,
    files = selected,
    pattern = paste0(
      family$filename_prefix,
      "[^_]*_",
      composite_product,
      "\\.(tif|tiff|vrt)$"
    )
  )
}

.force_map_raster <- function(
  FORCE_datacube,
  satellite_family,
  composite_product,
  time_steps
) {
  family <- .force_family_spec(satellite_family)
  if (inherits(
    FORCE_datacube,
    "SpatRaster"
  )) {
    return(
      .canonical_force_raster(
        FORCE_datacube,
        time_steps,
        family
      )
    )
  }

  product <- .force_map_composite_spec(
    FORCE_datacube,
    satellite_family,
    composite_product
  )
  if (length(product$files) != time_steps) {
    stop(
      "FORCE_datacube contains ",
      length(product$files),
      " ",
      product$product,
      " files for satellite_family = \"",
      family$argument,
      "\", but ts_ids contains ",
      time_steps,
      " time steps; each selected tile must contain one file per time step",
      call. = FALSE
    )
  }
  sources <- lapply(product$files, terra::rast)
  specs <- lapply(
    sources,
    .force_boa_band_spec,
    family = family
  )
  x <- do.call(
    c, lapply(
      seq_along(sources),
      function(i) {
        sources[[i]][[specs[[i]]$positions]]
      }
    )
  )
  names(x) <- rep(
    family$band_names,
    times = time_steps
  )
  x
}
