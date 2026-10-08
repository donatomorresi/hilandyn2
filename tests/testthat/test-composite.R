make_datacube <- function(
  acquisition_values,
  dates,
  bands = 10L,
  sensor = "SEN2A",
  qai_values = rep(0, length(dates)),
  output_dir = tempfile("force_datacube_"),
  dst_values = NULL
) {
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  band_names <- if (bands == 10L) {
    c(
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
  } else {
    c(
      "BLUE",
      "GREEN",
      "RED",
      "NIR",
      "SWIR1",
      "SWIR2"
    )
  }
  for (i in seq_along(dates)) {
    prefix <- paste0(
      format(dates[i], "%Y%m%d"),
      "_LEVEL2_",
      sensor
    )
    boa <- terra::rast(
      nrows = 1,
      ncols = 1,
      nlyrs = bands
    )
    terra::values(boa) <- matrix(
      acquisition_values[i],
      nrow = 1,
      ncol = bands
    )
    names(boa) <- band_names
    terra::writeRaster(
      boa,
      file.path(output_dir, paste0(prefix, "_BOA.tif")),
      overwrite = TRUE
    )
    qai <- terra::rast(boa[[1L]])
    terra::values(qai) <- qai_values[i]
    names(qai) <- "QAI"
    terra::writeRaster(
      qai,
      file.path(output_dir, paste0(prefix, "_QAI.tif")),
      overwrite = TRUE
    )
    if (!is.null(dst_values) && !is.na(dst_values[i])) {
      dst <- terra::rast(boa[[1L]])
      terra::values(dst) <- dst_values[i]
      names(dst) <- "DST"
      terra::writeRaster(
        dst,
        file.path(output_dir, paste0(prefix, "_DST.tif")),
        overwrite = TRUE
      )
    }
  }
  output_dir
}

make_force_definition <- function(path) {
  writeLines(
    c(
      "PROJECTION = EPSG:3857",
      "ORIGIN_GEO_X = 0",
      "ORIGIN_GEO_Y = 0",
      "ORIGIN_MAP_X = 0",
      "ORIGIN_MAP_Y = 200",
      "TILE_SIZE_X = 100",
      "TILE_SIZE_Y = 100"
    ),
    file.path(path, "datacube-definition.prj")
  )
}

test_that("zero-byte composite outputs are removed", {
  directory <- tempfile("zero_byte_outputs_")
  dir.create(directory)
  empty <- file.path(directory, "empty.tif")
  populated <- file.path(directory, "populated.tif")
  file.create(empty)
  writeBin(as.raw(1:4), populated)

  removed <- hilandyn2:::.force_remove_zero_byte_files(
    c(
      empty,
      populated,
      file.path(directory, "missing.tif")
    )
  )

  expect_identical(removed, empty)
  expect_false(file.exists(empty))
  expect_true(file.exists(populated))
})

test_that(
  "MED and GEO exclude invalid observations and share medoid-derived INF", {
    dates <- as.Date(c(
      "2020-06-01",
      "2020-06-11",
      "2020-06-21",
      "2020-06-25"
    ))
    datacube <- make_datacube(
      c(
        1,
        3,
        9,
        NA
      ),
      dates,
      qai_values = c(
        0,
        1,
        32,
        0
      )
    )
    for (method in c("medoid", "geomedian")) {
      result <- hilandyn2_composite(
        datacube,
        method = method,
        start_date = as.Date("2020-06-01"),
        end_date = as.Date("2020-06-30"),
        num_threads = 1
      )
      expect_equal(
        as.numeric(terra::values(result$reflectance)),
        rep(if (method == "medoid") 1 else 5, 10),
        tolerance = 1e-07
      )
      expect_equal(
        as.numeric(terra::values(result$information)),
        c(
          0,
          2,
          153,
          2020,
          -14,
          1
        )
      )
      for (qai_bits in list(c(
        0:5,
        8L,
        9L
      ), integer(0))) {
        filtered <- hilandyn2_composite(
          datacube,
          method = method,
          start_date = as.Date("2020-06-01"),
          end_date = as.Date("2020-06-30"),
          qai_bits = qai_bits,
          num_threads = 1
        )
        expect_equal(
          as.numeric(terra::values(filtered$reflectance)),
          rep(if (length(qai_bits)) 1 else 3, 10),
          tolerance = 1e-07
        )
        expect_equal(
          as.numeric(terra::values(filtered$information))[1:2],
          if (length(qai_bits)) c(0, 1) else c(1, 3)
        )
      }
    }
  }
)

test_that("weighted composites require matching DST files", {
  expect_identical(formals(hilandyn2_composite)$use_weights, FALSE)
  for (invalid in list(
    NULL,
    NA,
    logical(0),
    c(TRUE, FALSE),
    1,
    "TRUE"
  )) {
    expect_error(hilandyn2_composite(tempfile(), use_weights = invalid),
      "use_weights must be TRUE or FALSE")
  }
  dates <- as.Date(c(
    "2020-06-01",
    "2020-06-11",
    "2020-06-21",
    "2020-06-25"
  ))
  for (family in c("Sentinel-2", "Landsat")) {
    bands <- if (family == "Sentinel-2") 10L else 6L
    sensor <- if (family == "Sentinel-2") "SEN2A" else "LND08"
    missing <- make_datacube(
      c(
        1,
        3,
        5,
        9
      ),
      dates,
      bands,
      sensor
    )
    partial <- make_datacube(
      c(
        1,
        3,
        5,
        9
      ),
      dates,
      bands,
      sensor,
      dst_values = c(
        NA,
        NA,
        3000,
        3000
      )
    )
    complete <- make_datacube(
      c(
        1,
        3,
        5,
        9
      ),
      dates,
      bands,
      sensor,
      dst_values = c(
        0,
        0,
        3000,
        3000
      )
    )
    for (method in c("medoid", "geomedian")) {
      for (datacube in c(missing, partial)) {
        expect_error(hilandyn2_composite(
          datacube,
          satellite_family = family,
          method = method,
          use_weights = TRUE
        ), "use_weights = TRUE requires matching FORCE DST files")
      }
      unweighted <- hilandyn2_composite(
        partial,
        satellite_family = family,
        method = method,
        use_weights = FALSE
      )
      expect_equal(
        as.numeric(terra::values(unweighted$reflectance)),
        rep(if (method == "medoid") 3 else 4.5, bands),
        tolerance = 1e-07
      )
      weighted <- hilandyn2_composite(
        complete,
        satellite_family = family,
        method = method,
        use_weights = TRUE
      )
      expect_equal(
        as.numeric(terra::values(weighted$reflectance)),
        rep(5, bands),
        tolerance = 1e-07
      )
      expect_equal(as.numeric(terra::values(weighted$information))[1:3],
        c(
          0,
          4,
          173
        ))
      selected <- hilandyn2_composite(
        partial,
        satellite_family = family,
        method = method,
        use_weights = TRUE,
        start_date = dates[4L]
      )
      expect_equal(
        as.numeric(terra::values(selected$reflectance)),
        rep(9, bands),
        tolerance = 1e-07
      )
    }
  }
})

test_that("QAI screening validates bit positions and rejects selected flags", {
  defaults <- c(
    0:4,
    8L,
    9L
  )
  expect_identical(hilandyn2:::.force_qai_mask(defaults), 799L)
  expect_identical(hilandyn2:::.force_qai_mask(c(5L, 5L)), 32L)
  expect_identical(hilandyn2:::.force_qai_mask(integer(0)), 0L)
  expect_identical(hilandyn2:::.force_qai_mask(0:15), 65535L)
  for (invalid in list(
    -1,
    16,
    0.5,
    NA_real_,
    Inf,
    "4",
    TRUE,
    NULL,
    list(4),
    matrix(4),
    1i
  )) {
    expect_error(hilandyn2_composite(tempfile(), qai_bits = invalid),
      "qai_bits must be")
  }

  qai <- c(
    0,
    2^(0:15),
    6,
    34,
    NA_real_,
    Inf,
    0,
    0
  )
  for (bands in c(6L, 10L)) {
    input <- cbind(matrix(
      100,
      nrow = length(qai),
      ncol = bands
    ), qai)
    input[length(qai) - 1L, 1L] <- NA_real_
    input[length(qai), bands] <- Inf
    for (bits in list(
      defaults,
      c(defaults, 5L),
      integer(0)
    )) {
      mask <- hilandyn2:::.force_qai_mask(bits)
      valid <- is.finite(qai)
      valid[valid] <- bitwAnd(as.integer(qai[valid]), mask) == 0L
      valid[(length(qai) - 1L):length(qai)] <- FALSE
      result <- hilandyn2:::force_composite_chunk_cpp(
        as.numeric(input),
        nrow(input),
        ncol(input),
        bands,
        1L,
        153L,
        2020L,
        1L,
        matrix(c(1L, 1L), nrow = 1),
        1L,
        153L,
        1L,
        0L,
        FALSE,
        FALSE,
        1500,
        15,
        1L,
        qai_mask = mask
      )
      expect_true(all(is.na(result[!valid, , drop = FALSE])))
      expect_true(all(result[valid, seq_len(bands)] == 100))
      expect_equal(as.numeric(result[valid, bands + 1L]), qai[valid])
      expect_true(all(result[valid, bands + 2L] == 1))
    }
  }
})

test_that(
  "Landsat temporal bins preserve empty bins and FORCE output names", {
    root <- tempfile("force_root_")
    dir.create(root)
    make_force_definition(root)
    tile <- file.path(root, "X0000_Y0000")
    dates <- as.Date(c("2020-06-01", "2020-06-21"))
    make_datacube(
      c(1, 9),
      dates,
      bands = 6L,
      sensor = "LND08",
      output_dir = tile
    )
    make_datacube(
      c(100, 900),
      dates,
      output_dir = tile
    )
    output <- tempfile("force_output_")
    result <- hilandyn2_composite(
      root,
      satellite_family = "Landsat",
      method = "medoid",
      n_bins = 3,
      start_date = as.Date("2020-06-01"),
      end_date = as.Date("2020-06-30"),
      out_file = output,
      num_threads = 1
    )
    values <- as.numeric(terra::values(result$reflectance))
    expect_equal(values[1:6], rep(1, 6))
    expect_true(all(is.na(values[7:12])))
    expect_equal(values[13:18], rep(9, 6))
    expect_true(
      all(
        is.na(
          terra::values(result$information)[1,
          7:12]
        )
      )
    )
    prefix <- paste0(
      c(
        "20200605",
        "20200615",
        "20200625"
      ),
      "_LEVEL3_LNDLG_"
    )
    for (product in c("reflectance", "information")) {
      paths <- terra::sources(result[[product]])
      expect_true(all(file.exists(paths)))
      expect_equal(
        unname(basename(paths)),
        paste0(
          prefix, if (product == "reflectance")
          "MED.tif" else "INF.tif"
        )
      )
    }
    expect_true(file.exists(file.path(output, "datacube-definition.prj")))
  }
)

test_that(
  "hybrid bins reassign observations using temporal centres", {
    dates <- as.Date(c(
      "2020-06-01",
      "2020-06-02",
      "2020-06-03",
      "2020-06-30"
    ))
    datacube <- make_datacube(
      c(
        1,
        3,
        5,
        9
      ),
      dates
    )

    result <- hilandyn2_composite(
      datacube,
      method = "medoid",
      n_bins = 2,
      bin_method = "hybrid",
      start_date = as.Date("2020-06-01"),
      end_date = as.Date("2020-06-30"),
      hybrid_max_days = 1,
      num_threads = 1
    )
    result_values <- terra::values(result$reflectance)[1,
      ]

    expect_equal(
      as.numeric(result_values[1:10]),
      rep(3, 10)
    )
    expect_equal(
      as.numeric(result_values[11:20]),
      rep(9, 10)
    )
  }
)
