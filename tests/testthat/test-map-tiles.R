make_map_tile_data <- function(steps = 8L) {
  root <- tempfile("map_tiles_")
  dir.create(root)
  bands <- c(
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
  tile_specs <- list(
    X0000_Y0000 = list(
      extent = c(
        4771000,
        4771030,
        4549716.2554659806,
        4549746.2554659806
      ),
      offset = 0
    ),
    X0001_Y0000 = list(
      extent = c(
        4771030,
        4771060,
        4549716.2554659806,
        4549746.2554659806
      ),
      offset = 10
    )
  )
  for (tile in names(tile_specs)) {
    directory <- file.path(root, tile)
    dir.create(directory)
    specification <- tile_specs[[tile]]
    for (step in seq_len(steps)) {
      x <- terra::rast(
        nrows = 3,
        ncols = 3,
        nlyrs = length(bands),
        xmin = specification$extent[1],
        xmax = specification$extent[2],
        ymin = specification$extent[3],
        ymax = specification$extent[4],
        crs = "EPSG:3857"
      )
      band_values <- seq(
        1000,
        10000,
        by = 1000
      ) +
        specification$offset + step
      terra::values(x) <- matrix(
        rep(band_values, each = terra::ncell(x)), nrow = terra::ncell(x)
      )
      names(x) <- bands
      date <- as.Date(sprintf("%d-07-01", 2016 + step))
      terra::writeRaster(
        x,
        file.path(directory, paste0(
          format(date, "%Y%m%d"), "_LEVEL3_SEN2L_MED.tif"
        )),
        overwrite = TRUE
      )
    }
  }
  root
}

test_that("adjacent FORCE tiles populate the extended focal window", {
  root <- make_map_tile_data()
  product <- hilandyn2:::.force_map_composite_spec(
    root,
    "Sentinel-2",
    "MED"
  )
  paths <- hilandyn2:::.force_tile_directories(
    root,
    product$pattern,
    "FORCE_datacube"
  )
  files <- hilandyn2:::.map_tile_file_table(
    paths,
    product$pattern,
    "FORCE_datacube"
  )
  neighbors <- hilandyn2:::.map_adjacent_tiles(
    "X0000_Y0000", names(paths)
  )
  input <- hilandyn2:::.map_tile_vrt_input(
    files[neighbors],
    "X0000_Y0000",
    1L,
    "FORCE_datacube"
  )

  first_layer <- terra::as.matrix(input[[1L]], wide = TRUE)
  expect_equal(dim(first_layer), c(5L, 5L))
  expect_true(all(is.na(first_layer[, 1L])))
  expect_equal(first_layer[2:4, 5L], rep(1011, 3))
  expect_equal(terra::nlyr(input), 80L)
})

test_that("vector AOI output retains the cropped target geometry", {
  root <- make_map_tile_data(steps = 16L)
  target_file <- list.files(
    file.path(root, "X0000_Y0000"), full.names = TRUE
  )[1L]
  target <- terra::rast(target_file)[[1L]]
  target_extent <- terra::ext(target)
  aoi <- terra::as.polygons(terra::ext(
    terra::xmin(target_extent),
    terra::xmin(target_extent) + 20,
    terra::ymax(target_extent) - 20,
    terra::ymax(target_extent)
  ))
  terra::crs(aoi) <- terra::crs(target)
  aoi_file <- tempfile(fileext = ".gpkg")
  terra::writeVector(
    aoi,
    aoi_file,
    overwrite = TRUE
  )

  result <- hilandyn2_map(
    FORCE_datacube = root,
    satellite_family = "Sentinel-2",
    composite_product = "MED",
    out_path = tempfile("map_tile_aoi_"),
    tiles = "X0000_Y0000",
    include_adjacent = FALSE,
    roi_vec = aoi_file,
    ts_ids = as.numeric(paste0(2017:2032, "01")),
    spectral_bands = character(0),
    spectral_indices = "NDVI",
    num_cores = 1L,
    num_threads = 1L,
    win_side = 3L,
    cell_weights = FALSE,
    rmse_out = FALSE,
    est_mat_out = FALSE,
    wavelet_level = 0
  )

  expected <- terra::crop(target, aoi)
  expect_true(terra::compareGeom(
    result,
    expected,
    crs = TRUE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE,
    stopOnError = FALSE
  ))
})
