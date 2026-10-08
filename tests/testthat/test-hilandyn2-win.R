monthly_ids <- function(n, start = as.Date("2020-01-01")) {
  as.integer(format(seq(
    start,
    by = "month",
    length.out = n
  ), "%Y%m"))
}

test_that("hilandyn2_win detects a known downward change", {
  nt <- 72L
  ts_ids <- monthly_ids(nt)
  time <- seq_len(nt)
  signal <- matrix(
    0.7 + 0.05 * sin(2 * pi * time / 12) - 0.3 * (time > 36L),
    nrow = 1L
  )

  result <- hilandyn2_win(
    signal,
    nb = 1L,
    nc = 1L,
    nt = nt,
    nr = 1L,
    ts_ids = ts_ids,
    cng_dir = -1,
    cell_weights = FALSE,
    noise_zscore = 0,
    wavelet_level = 0L,
    est_mat_out = TRUE,
    verbose = FALSE
  )

  expect_equal(unname(result$D_NUM), 1)
  expect_equal(unname(result$D_MAX_TS), ts_ids[37L])
  expect_equal(which(unname(result$CPT_ID) == 101), 37L)
  expect_equal(dim(result$EST), c(1L, nt))
})

test_that("MODWPT reduces seasonal structure and respects the ACF threshold", {
  nt <- 72L
  ts_ids <- monthly_ids(nt)
  time <- seq_len(nt)
  seasonal_signal <- matrix(
    0.7 + 0.15 * sin(2 * pi * time / 12) +
      0.02 * cos(2 * pi * time / 6) + 0.003 * time,
    nrow = 1L
  )
  seasonal_reference <- seasonal_signal * 1

  adjusted <- hilandyn2_win(
    seasonal_signal,
    nb = 1L,
    nc = 1L,
    nt = nt,
    nr = 1L,
    ts_ids = ts_ids,
    cng_dir = -1,
    cell_weights = FALSE,
    noise_zscore = 0,
    wavelet_filters = c("LA4", "LA8"),
    wavelet_level = -1L,
    wavelet_filter_penalty = 0.1,
    acf_min = 0.2,
    sgn_mat_out = TRUE,
    verbose = FALSE
  )

  expect_true(is.finite(adjusted$S_ORIG))
  expect_true(is.finite(adjusted$S_ADJ))
  expect_gt(unname(adjusted$S_ORIG - adjusted$S_ADJ), 0.05)
  expect_gt(max(abs(adjusted$SGN - seasonal_reference)), 0.05)

  weak_signal <- matrix(
    0.7 + 0.003 * time + 0.001 * sin(1.7 * time),
    nrow = 1L
  )
  weak_reference <- weak_signal * 1
  skipped <- hilandyn2_win(
    weak_signal,
    nb = 1L,
    nc = 1L,
    nt = nt,
    nr = 1L,
    ts_ids = ts_ids,
    cng_dir = -1,
    cell_weights = FALSE,
    noise_zscore = 0,
    wavelet_filters = c("LA4", "LA8"),
    wavelet_level = -1L,
    wavelet_filter_penalty = 0.1,
    acf_min = 0.2,
    sgn_mat_out = TRUE,
    verbose = FALSE
  )

  expect_lt(unname(skipped$S_ORIG), 0.2)
  expect_equal(
    unname(skipped$S_ADJ),
    unname(skipped$S_ORIG),
    tolerance = 1e-12
  )
  expect_equal(
    skipped$SGN,
    weak_reference,
    tolerance = 1e-12
  )
})

test_that("threshold search reports when no candidate satisfies cpt_max", {
  nt <- 72L
  time <- seq_len(nt)
  signal <- matrix(
    0.7 + 0.05 * sin(2 * pi * time / 12) - 0.3 * (time > 36L),
    nrow = 1L
  )

  expect_error(
    hilandyn2_win(
      x = signal,
      nb = 1L,
      nc = 1L,
      nt = nt,
      nr = 1L,
      ts_ids = monthly_ids(nt),
      cng_dir = -1,
      cell_weights = FALSE,
      noise_zscore = 0,
      wavelet_level = 0L,
      th_min = 0,
      th_max = 1e-6,
      cpt_max = 0L,
      verbose = FALSE
    ),
    "No candidate satisfies cpt_max; increase th_max or cpt_max.",
    fixed = TRUE
  )
})

