#' High-dimensional detection of Landscape Dynamics 2 engine
#'
#' \code{hilandyn2_win()} segments a high-dimensional \emph{N x T} time-series
#' matrix obtained from a single spatial kernel, where \emph{N} is the
#' number of input bands multiplied by the number of cells and \emph{T}
#' is the number of time points. Before segmentation, it can filter
#' impulsive noise using a robust baseline and adjust seasonality using
#' the maximal overlap discrete wavelet packet transform (MODWPT).
#'
#' Typically, \code{hilandyn2_win()} is used to test the algorithm using
#' satellite time-series data (Sentinel-2 and Landsat) extracted from a single 
#' spatial kernel. It processes intra-annual high-dimensional time series to 
#' detect changes in spectral trends. Changes in the intercept, slope or both of
#' linear trends are detected using the High-dimensional Trend Segmentation
#' (HiTS) procedure proposed by
#' \insertCite{maeng2019adaptive;textual}{hilandyn2}.
#' The HiTS procedure aims at detecting changepoints in a piecewise linear
#' signal where their number and location are unknown.
#' High-dimensional time series can include single or multiple reflectance
#' bands and spectral indices, hereafter referred to as bands.
#' Input data can be a \code{numeric} \code{vector} extracted from a
#' \code{SpatRaster} or a \code{matrix}. Each row corresponds to a band
#' and a cell in the spatial kernel.
#'
#' Impulsive noise is identified from deviations relative to a robust
#' temporal baseline and filtered before segmentation. The
#' \code{noise_zscore} parameter controls the filtering threshold; use zero
#' to disable noise filtering. Missing values are filled from spatial
#' neighbors where possible, then from temporally nearby observations.
#'
#' Optional MODWPT seasonal adjustment reduces recurring within-year
#' variation before HiTS segmentation. Candidate wavelet filters are
#' specified by \code{wavelet_filters}; \code{wavelet_level = NULL} selects
#' the decomposition level automatically, while \code{wavelet_level = 0}
#' disables seasonal adjustment. The \code{acf_min} parameter sets the
#' minimum seasonal autocorrelation required to attempt adjustment.
#'
#' @param x numeric vector or matrix. Input time series where each row of the
#'   matrix contains time-ordered data relative to a band and a cell in the
#'   spatial kernel.
#' @param nb integer. Number of input bands in the time series.
#' @param nc integer. The number of cells within the spatial kernel.
#' @param nt integer. Number of time steps in the time series.
#' @param nr integer. Number of rows (variables) in the high-dimensional time
#'   series.
#' @param ts_ids numeric vector of length \code{nt}. Increasing, unique
#'   six-digit identifiers containing the four-digit year followed by the
#'   two-digit position within the year. See the example for more details.
#' @param cng_dir numeric vector of length \code{nr}. Direction of
#'   disturbance-related change for each row of \code{x}: 1 for an increase
#'   and -1 for a decrease from pre-change to post-change observations.
#'   Repeat each band's direction \code{nc} times for a multi-cell kernel.
#' @param cell_weights logical. Whether to compute cell-based weights within the
#'   spatial kernel.
#' @param ts_len_min integer. Minimum number of observations in the input time
#'   series. Time series shorter than this value are ignored by the algorithm.
#' @param gap_len_max integer. Maximum number of missing consecutive
#'   observations allowed. Time series with longer data gaps are ignored by the
#'   algorithm.
#' @param noise_zscore Z-score threshold used by the impulsive-noise filter.
#' @param wavelet_filters Character vector of candidate wavelet filters.
#' @param wavelet_level Optional wavelet decomposition level. Use `0` to disable
#'   seasonal adjustment.
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
#' @param rmse_out,len_seg_out,slo_seg_out,bln_mat_out
#'   Logical switches controlling optional output elements.
#' @param zsc_mat_out,mag_mat_out,sgn_mat_out,est_mat_out
#'   Logical switches controlling optional output elements.
#' @param nout Retained for compatibility. The output length is computed
#'   internally from the requested outputs; this argument is ignored.
#' @param verbose Logical; print model-selection progress.
#' @param debug Logical; save the input data associated with a processing error
#'   to an RDS file in the working directory before rethrowing the error.
#'
#'
#' @return A list containing the elements below. Disabled optional outputs are
#'   \code{NULL}; per-band time-series outputs are \code{nb} by \code{nt}
#'   matrices. Unprocessable windows return missing values for enabled outputs.
#'   Invalid input and unexpected processing failures raise errors.
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
#' @references
#' \insertRef{Morresi2024}{hilandyn2}
#'
#' @examples
#' data(s2_windows)
#' bands <- c("NBR", "MSI", "CRSWIR", "CRRE", "NMDI")
#' nc <- 9L # 3 by 3 kernel; nine consecutive rows per index
#' ts_ids <- as.numeric(paste0(rep(2016:2025, each = 4),
#'                             sprintf("%02d", rep(1:4, times = 10))))
#' cng_dir <- rep(c(-1, 1, 1, 1, -1), each = nc)
#' 
#' # Analyse each window with and without automatic seasonal adjustment.
#' results <- vector("list", length(s2_windows))
#' names(results) <- names(s2_windows)
#' for (i in seq_along(s2_windows)) {
#'   results[[i]] <- lapply(c(0L, -1L), function(level) {
#'     hilandyn2_win(
#'       x = s2_windows[[i]] * 1, # fresh copy for each analysis
#'       nb = length(bands), nc = nc, nt = length(ts_ids),
#'       nr = length(bands) * nc, ts_ids = ts_ids, cng_dir = cng_dir,
#'       cell_weights = TRUE, sd_method = "mean_abs_dev",
#'       noise_zscore = 99, th_min = 1.4, th_max = 1.4,
#'       wavelet_filters = c("LA4", "LA6", "LA8"), wavelet_level = level,
#'       sgn_mat_out = TRUE, est_mat_out = TRUE, verbose = FALSE
#'     )
#'   })
#'   names(results[[i]]) <- c("Unadjusted", "Seasonal adjustment enabled")
#' }
#' 
#' # Each window produces a page with paired panels and shared index scales.
#' old_par <- par(mfrow = c(length(bands), 2), mar = c(2, 3, 2, 1),
#'                oma = c(2, 0, 2, 0))
#' ticks <- seq(1, length(ts_ids), by = 8)
#' change_cols <- c("101" = "firebrick", "102" = "darkgreen", "9" = "grey50")
#' for (i in seq_along(results)) {
#'   pair <- results[[i]]
#'   for (b in seq_along(bands)) {
#'     limits <- range(pair[[1]]$SGN[b, ], pair[[1]]$EST[b, ],
#'                     pair[[2]]$SGN[b, ], pair[[2]]$EST[b, ], na.rm = TRUE)
#'     for (j in seq_along(pair)) {
#'       out <- pair[[j]]
#'       plot(out$SGN[b, ], type = "p", pch = 1, cex = 0.6,
#'            col = "red", ylim = limits, xaxt = "n",
#'            xlab = "", ylab = bands[b])
#'       lines(out$EST[b, ], col = "blue", lwd = 2)
#'       changes <- which(!is.na(out$CPT_ID))
#'       abline(v = changes, lty = 2,
#'              col = change_cols[as.character(out$CPT_ID[changes])])
#'       axis(1, at = ticks, labels = ts_ids[ticks], cex.axis = 0.7)
#'       if (b == 1L) title(main = names(pair)[j], cex.main = 0.8)
#'       if (b == length(bands)) {
#'         legend("bottomleft", c("Input", "Fitted"),
#'                col = c("red", "blue"), pch = c(1, NA),
#'                lty = c(NA, 1), bty = "n", cex = 0.7)
#'       }
#'     }
#'   }
#'   mtext(names(results)[i], side = 3, outer = TRUE)
#'   mtext("Time ID; dashed lines: disturbance (red), greening (green), other (grey)",
#'         side = 1, outer = TRUE, cex = 0.7)
#' }
#' par(old_par)
#'
#' @export

hilandyn2_win <- function(
  x,
  nb,
  nc,
  nt,
  nr,
  ts_ids,
  cng_dir,
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
  rmse_out = FALSE,
  len_seg_out = FALSE,
  slo_seg_out = FALSE,
  bln_mat_out = FALSE,
  zsc_mat_out = FALSE,
  mag_mat_out = FALSE,
  sgn_mat_out = FALSE,
  est_mat_out = FALSE,
  nout = 1L,
  verbose = TRUE,
  debug = FALSE
)
{
  
  dimensions <- list(
    nb = nb,
    nc = nc,
    nt = nt,
    nr = nr
  )
  for (argument in names(dimensions)) {
    value <- dimensions[[argument]]
    if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
        !is.finite(value) || value < 1 || value != floor(value)) {
      stop(
        argument,
        " must be one positive integer",
        call. = FALSE
      )
    }
  }
  if (nr != nb * nc) {
    stop("nr must equal nb * nc", call. = FALSE)
  }
  if (!is.numeric(x) || is.complex(x) || length(x) != nr * nt) {
    stop("x must contain nr * nt numeric values", call. = FALSE)
  }
  if (!is.numeric(ts_ids) || is.complex(ts_ids) || length(ts_ids) != nt ||
      any(!is.finite(ts_ids)) || any(ts_ids <= 0 | ts_ids != floor(ts_ids)) ||
      anyDuplicated(ts_ids) || is.unsorted(ts_ids)) {
    stop("ts_ids must contain nt increasing, unique positive integers",
         call. = FALSE)
  }
  if (!is.numeric(cng_dir) || is.complex(cng_dir) || length(cng_dir) != nr ||
      anyNA(cng_dir) || any(!cng_dir %in% c(-1, 1))) {
    stop("cng_dir must contain nr values, each either -1 or 1", call. = FALSE)
  }

  dim(x) <- c(nr, nt)
  
  
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
  
  out_nm <- valn
  if (length(vec1n) > 0L)
    out_nm <- c(
      out_nm, paste0(
        rep(vec1n, each = nb),
        "_",
        seq(1, nb)
      )
    )
  out_nm <- c(
    out_nm, paste0(
      rep(vec2n, each = nt),
      "_",
      ts_ids
    )
  )
  if (length(matn) > 0L)
    out_nm <- c(
      out_nm, paste0(
        rep(
          matn, each = nb *
            nt
        ),
        "_",
        seq(1, nb),
        "_",
        rep(ts_ids, each = nb)
      )
    )
  
  # Assign values to NULL parameters
  if (is.null(th_max))
    th_max <- -1
  if (is.null(seg_len_min))
    seg_len_min <- -1
  if (is.null(wavelet_level))
    wavelet_level <- -1
  
  nout <- length(out_nm)
  
  tryCatch(
    {
      v <- hilandyn2_int_cpp(
        x,
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
        verbose
      )
      
      names(v) <- out_nm
      
      out_list <- list(
        D_MAX_TS = v[1],
        D_MAX_MG = v[2],
        D_FRS_TS = v[3],
        D_FRS_MG = v[4],
        D_LST_TS = v[5],
        D_LST_MG = v[6],
        D_NUM = v[7],
        G_MAX_TS = v[8],
        G_MAX_MG = v[9],
        G_FRS_TS = v[10],
        G_FRS_MG = v[11],
        G_LST_TS = v[12],
        G_LST_MG = v[13],
        G_NUM = v[14],
        S_ORIG = v[15],
        S_ADJ = v[16],
        TH_CONST = v[17],
        RMSE = if (any(grepl("RMSE", names(v)))) v[grep("RMSE", names(v))],
        CPT_ID = v[grep("CPT_ID", names(v))],
        REL_MAG = v[grep("REL_MAG", names(v))],
        LEN = if (any(grepl("LEN", names(v)))) v[grep("LEN", names(v))],
        SLO = if (any(grepl("SLO", names(v)))) 
          matrix(
            v[grep("SLO", names(v))],
            nb,
            nt
          ),
        BLN = if (any(grepl("BLN", names(v)))) 
          matrix(
            v[grep("BLN", names(v))],
            nb,
            nt
          ),
        ZSC = if (any(grepl("ZSC", names(v)))) 
          matrix(
            v[grep("ZSC", names(v))],
            nb,
            nt
          ),
        MAG = if (any(grepl("^MAG", names(v)))) 
          matrix(
            v[grep("^MAG", names(v))],
            nb,
            nt
          ),
        SGN = if (any(grepl("SGN", names(v)))) 
          matrix(
            v[grep("SGN", names(v))],
            nb,
            nt
          ),
        EST = if (any(grepl("EST", names(v)))) 
          matrix(
            v[grep("EST", names(v))],
            nb,
            nt
          )
      )
      
      return(out_list)
      
    }, error = function(e) {
      
      if (debug) {
        timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
        debug_filename <- paste0(
          "input_data_",
          timestamp,
          ".rds"
        )
        
        saveRDS(x, file = debug_filename)
        message(
          "Input data has been saved to: ",
          debug_filename
        )
      }
      stop(e)
    }
  )
}
