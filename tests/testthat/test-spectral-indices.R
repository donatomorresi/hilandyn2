test_that(
  "Sentinel-2 index formulas and ordering are correct",
  {
    source <- matrix(
      seq(
        1000,
        10000,
        by = 1000
      ),
      nrow = 1
    )
    specification <- hilandyn2:::.force_analysis_spec(
      "Sentinel-2",
      character(),
      c(
        "NBR",
        "NDVI",
        "MSI"
      )
    )
    result <- hilandyn2:::.force_map_chunk_values(
      as.numeric(source),
      1L,
      10L,
      1L,
      integer(),
      specification$index_ids,
      specification$family$id,
      TRUE,
      1L
    )

    expect_equal(
      specification$index_names,
      c(
        "NDVI",
        "NBR",
        "MSI"
      )
    )
    expect_equal(
      as.numeric(result),
      c(
        0.4,
        -1/9,
        1.125
      ),
      tolerance = 1e-07
    )
  }
)

test_that(
  "FORCE Landsat indices use the six legacy bands",
  {
    source <- matrix(
      seq(
        1000,
        6000,
        by = 1000
      ),
      nrow = 1
    )
    specification <- hilandyn2:::.force_analysis_spec(
      "Landsat",
      character(),
      c(
        "NBR",
        "NDVI",
        "MSI"
      )
    )
    result <- hilandyn2:::.force_map_chunk_values(
      as.numeric(source),
      1L,
      6L,
      1L,
      integer(),
      specification$index_ids,
      specification$family$id,
      TRUE,
      1L
    )

    expect_equal(
      specification$index_names,
      c(
        "NDVI",
        "NBR",
        "MSI"
      )
    )
    expect_equal(
      as.numeric(result),
      c(
        1/7,
        -0.2,
        1.25
      ),
      tolerance = 1e-07
    )

    tc_specification <- hilandyn2:::.force_analysis_spec(
      "Landsat",
      character(),
      c(
        "TCA",
        "TCW",
        "TCG",
        "TCB"
      )
    )
    tasseled_cap <- hilandyn2:::.force_map_chunk_values(
      as.numeric(source),
      1L,
      6L,
      1L,
      integer(),
      tc_specification$index_ids,
      tc_specification$family$id,
      TRUE,
      1L
    )
    brightness <- sum(
      (1:6/10) * c(
        0.3037,
        0.2793,
        0.4743,
        0.5585,
        0.5082,
        0.1863
      )
    )
    greenness <- sum(
      (1:6/10) * c(
        -0.2848,
        -0.2435,
        -0.5436,
        0.7243,
        0.084,
        -0.18
      )
    )
    wetness <- sum(
      (1:6/10) * c(
        0.1509,
        0.1973,
        0.3279,
        0.3406,
        -0.7112,
        -0.4572
      )
    )
    expect_equal(
      tc_specification$index_names,
      c(
        "TCB",
        "TCG",
        "TCW",
        "TCA"
      )
    )
    expect_equal(
      as.numeric(tasseled_cap),
      c(
        brightness,
        greenness,
        wetness,
        atan(greenness/brightness)
      ),
      tolerance = 1e-07
    )
    expect_error(
      hilandyn2:::.force_analysis_spec(
        "Landsat",
        character(),
        "NDRE1"
      ),
      "unknown names: NDRE1"
    )
  }
)

test_that(
  "change directions are inferred and explicit overrides are preserved", {
    variables <- c(
      "BLUE",
      "NIR",
      "NDMI",
      "MSI"
    )
    expect_equal(
      hilandyn2:::.resolve_change_directions(
        variables,
        NULL,
        2L,
        2L
      ),
      c(
        1,
        -1,
        -1,
        1
      )
    )
    expect_equal(
      hilandyn2:::.resolve_change_directions(
        variables,
        c(-1, 1),
        2L,
        2L
      ),
      c(
        -1,
        1,
        -1,
        1
      )
    )
    expect_equal(
      hilandyn2:::.resolve_change_directions(
        variables,
        rep(1, 4),
        2L,
        2L
      ),
      rep(1, 4)
    )
  }
)
