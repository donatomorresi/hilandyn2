test_that(
  "band and index masking preserve valid cells and time-step ordering", {
    source <- matrix(
      seq(
        500,
        by = 13,
        length.out = 120
      ),
      nrow = 6
    )
    expected_indices <- do.call(
      cbind, lapply(
        0:1, function(step) {
          offset <- step * 10L
          cbind(
            source[, offset + c(1L, 3L)],
            (source[, offset + 7L] - source[, offset + 3L])/(source[, offset +
            7L] + source[, offset + 3L]),
            (source[, offset + 8L] - source[,
            offset + 10L])/(source[, offset + 8L] + source[, offset + 10L])
          )
        }
      )
    )
    masks <- list(
      rep(1, 6),
      c(
        0,
        Inf,
        NA,
        NaN,
        -Inf,
        1
      ),
      rep(NA_real_, 6)
    )
    for (mask in masks) {
      results <- list(
        hilandyn2:::mask_chunk_values_cpp(
          as.numeric(source),
          mask,
          6L,
          20L
        ),
        hilandyn2:::.force_map_chunk_values(
          as.numeric(source),
          6L,
          20L,
          2L,
          c(1L, 3L),
          c(0L, 5L),
          0L,
          FALSE,
          1L,
          mask
        )
      )
      for (i in seq_along(results)) {
        expected <- list(source, expected_indices)[[i]]
        expected[is.na(mask),
          ] <- NA_real_
        actual <- matrix(results[[i]], nrow = 6)
        expect_equal(
          is.na(actual),
          is.na(expected)
        )
        expect_equal(
          actual[!is.na(expected)],
          expected[!is.na(expected)]
        )
      }
    }
  }
)

test_that(
  "mask coverage warns and retains uncovered cells but rejects a wrong grid",
  {
    target <- terra::rast(
      nrows = 4,
      ncols = 5,
      xmin = 0,
      xmax = 5,
      ymin = 0,
      ymax = 4,
      crs = "EPSG:3857"
    )
    mask <- terra::rast(
      nrows = 2,
      ncols = 2,
      xmin = 1,
      xmax = 3,
      ymin = 1,
      ymax = 3,
      crs = terra::crs(target)
    )
    terra::values(mask) <- c(
      NA,
      0,
      7,
      NA
    )
    path <- tempfile(fileext = ".tif")
    mask <- terra::writeRaster(
      mask,
      path,
      overwrite = TRUE
    )
    original <- terra::values(mask)
    expect_warning(
      specification <- hilandyn2:::.map_raster_mask_spec(mask, target),
      "does not fully cover.*value 1"
    )
    terra::readStart(mask)
    on.exit(
      terra::readStop(mask),
      add = TRUE
    )
    expected <- matrix(
      1,
      nrow = 6,
      ncol = 7
    )
    expected[3:4, 3:4] <- matrix(
      c(
        NA,
        0,
        7,
        NA
      ),
      nrow = 2,
      byrow = TRUE
    )
    actual <- matrix(
      hilandyn2:::.map_read_mask_chunk(
        specification,
        1L,
        6L,
        7L,
        1L
      ),
      nrow = 6,
      byrow = TRUE
    )
    expect_equal(
      is.na(actual),
      is.na(expected)
    )
    expect_equal(
      actual[!is.na(expected)],
      expected[!is.na(expected)]
    )
    expect_equal(
      terra::values(terra::rast(path)),
      original
    )

    disjoint <- terra::rast(
      nrows = 2,
      ncols = 2,
      xmin = 10,
      xmax = 12,
      ymin = 10,
      ymax = 12,
      crs = terra::crs(target)
    )
    terra::values(disjoint) <- NA_real_
    expect_warning(
      specification <- hilandyn2:::.map_raster_mask_spec(disjoint, target),
      "does not fully cover"
    )
    terra::readStart(disjoint)
    on.exit(
      terra::readStop(disjoint),
      add = TRUE
    )
    expect_equal(
      hilandyn2:::.map_read_mask_chunk(
        specification,
        1L,
        4L,
        5L
      ),
      rep(1, 20)
    )
    wrong_grid <- terra::rast(
      nrows = 2,
      ncols = 2,
      xmin = 0.5,
      xmax = 2.5,
      ymin = 0,
      ymax = 2,
      crs = terra::crs(target)
    )
    expect_error(
      hilandyn2:::.map_raster_mask_spec(wrong_grid, target),
      "same CRS"
    )
  }
)
