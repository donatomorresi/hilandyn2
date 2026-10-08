#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <vector>
#include "hilandyn2_misc.h"
#include "hits.h"
#include "hits_misc.h"
#include "hilandyn2_deseason.h"


arma::umat rle(const arma::uvec& x) {
  if (x.is_empty()) return arma::umat(0, 4);

  const arma::uvec change_points = arma::find(
    x.tail(x.n_elem - 1) != x.head(x.n_elem - 1));
  const arma::uvec starts = arma::join_vert(
    arma::uvec{0}, change_points + 1);
  const arma::uvec ends = arma::join_vert(
    change_points, arma::uvec{x.n_elem - 1});

  arma::umat out(starts.n_elem, 4);
  out.col(0) = x.elem(starts);
  out.col(1) = ends - starts + 1;
  out.col(2) = starts;
  out.col(3) = ends;
  return out;
}


arma::uvec ordered_temporally_closest_candidates(
  const arma::uword t,
  const arma::uword n,
  const arma::uword period
) {
  
  // Build candidate column indices ordered by temporal closeness to t:
  // 1) same-year columns first (excluding t), sorted by |c - t|
  // 2) other-year columns next, sorted by |c - t|
  
  if (n <= 1 || period == 0) return arma::uvec();
  
  const arma::uword year_t   = t / period;
  const arma::uword y_start  = year_t * period;
  const arma::uword y_end_ex = std::min(y_start + period, n);
  
  // all columns
  arma::uvec all_cols = arma::regspace<arma::uvec>(0, n - 1);
  
  // same year mask
  arma::uvec same_year =
      arma::find((all_cols >= y_start) % (all_cols < y_end_ex));

  // remove t from same_year
  arma::uvec same_year_no_t = same_year.elem(arma::find(same_year != t));
  
  // other years = all except same-year range
  arma::uvec other_years =
      arma::find((all_cols < y_start) + (all_cols >= y_end_ex));

  // sort same_year_no_t by |col - t|
  arma::vec dist_same = arma::abs(
      arma::conv_to<arma::vec>::from(same_year_no_t) - static_cast<double>(t));
  arma::uvec ord_same = arma::sort_index(dist_same, "ascend");
  arma::uvec same_sorted = same_year_no_t.elem(ord_same);
  
  // sort other_years by |col - t|
  arma::vec dist_other = arma::abs(arma::conv_to<arma::vec>::from(other_years) -
                                   static_cast<double>(t));
  arma::uvec ord_other = arma::sort_index(dist_other, "ascend");
  arma::uvec other_sorted = other_years.elem(ord_other);
  
  // concatenate [same_sorted ; other_sorted]
  arma::uvec out(same_sorted.n_elem + other_sorted.n_elem);
  if (same_sorted.n_elem > 0) {
    out.subvec(0, same_sorted.n_elem - 1) = same_sorted;
  }
  if (other_sorted.n_elem > 0) {
    out.subvec(same_sorted.n_elem, out.n_elem - 1) = other_sorted;
  }
  return out;
}

void fill_seasonal_nan(
  arma::mat &Y,
  const arma::uword period,
  const arma::uword min_obs
) {

  if (period == 0)   throw std::invalid_argument("period must be positive");
  if (min_obs == 0) throw std::invalid_argument("min_obs must be positive");

  const arma::uword p = Y.n_rows;
  const arma::uword n = Y.n_cols;
  if (n == 0 || p == 0) return;

  for (arma::uword r = 0; r < p; ++r) {
    // Row fallback median from original row finite values
    arma::rowvec row0 = Y.row(r);
    arma::uvec finite0 = arma::find_finite(row0);

    const bool has_row_fallback = !finite0.is_empty();
    const double row_fallback =
        has_row_fallback ? arma::as_scalar(arma::median(row0.elem(finite0)))
                         : 0.0;

    // Missing indices in this row
    arma::uvec miss_idx = arma::find_nonfinite(Y.row(r));

    for (arma::uword k = 0; k < miss_idx.n_elem; ++k) {
      const arma::uword t = miss_idx(k);

      // ordered candidate columns by temporal closeness rule
      arma::uvec cand_cols =
          ordered_temporally_closest_candidates(t, n, period);
      if (cand_cols.is_empty()) {
        Y(r, t) = row_fallback;
        continue;
      }

      // candidate values
      arma::rowvec cand_vals = Y(arma::uvec{r}, cand_cols);
      arma::uvec finite_cand_pos = arma::find_finite(cand_vals);

      if (finite_cand_pos.n_elem >= min_obs) {
        // take only the first min_obs finite candidates
        arma::uvec take_pos = finite_cand_pos.head(min_obs);
        arma::vec take_vals = cand_vals.elem(take_pos);

        Y(r, t) = arma::as_scalar(arma::mean(take_vals));
      } else {
        // not enough finite neighbors
        Y(r, t) = row_fallback;
      }
    }
  }
}

double mad(arma::vec x) {
  const double constant = 1.4826;
  arma::vec med(x.n_elem, arma::fill::value(median(x)));
  return arma::median(arma::abs(x - med)) * constant;
}


arma::vec mad_rows(
  const arma::mat& X,
  double consistency_c,
  double eps
) {
  arma::vec med = arma::median(X, 1);                    // d x 1
  arma::mat abs_dev = arma::abs(X.each_col() - med);    // d x n
  arma::vec mad = arma::median(abs_dev, 1);             // d x 1
  mad = arma::clamp(mad, eps, arma::datum::inf);
  return consistency_c * mad;
}


double bivar(const arma::vec& x) {

  const double constant = 1.4826;
  double med = arma::median(x);
  arma::vec med_vec(x.n_elem, arma::fill::value(med));
  arma::vec med_dif = x - med_vec;
  double mad = arma::median(arma::abs(med_dif)) * constant;
  arma::vec u = arma::abs(
      med_dif / arma::vec(x.n_elem, arma::fill::value(9 * 0.6744898 * mad)));
  arma::vec av(x.n_elem, arma::fill::zeros);
  av.elem(arma::find(u < 1)).ones();
  arma::vec u2 = arma::square(u);
  arma::vec five_vec = arma::vec(x.n_elem, arma::fill::value(5));

  double top = x.n_elem * arma::sum(av % arma::square(med_dif) %
                                    arma::pow((arma::ones(x.n_elem) - u2), 4));
  double bot = arma::sum(av % (arma::ones(x.n_elem) - u2) %
                         (arma::ones(x.n_elem) - five_vec % u2));
  double bi = top / std::pow(bot, 2);
  return bi;
}


arma::vec sd_est(const arma::mat& X,
                 const std::string& method) {
  
  arma::mat dif_m = arma::diff(X, 2, 1);
  arma::vec sd(X.n_rows);
  
  if (method == "median_abs_dev") {
    sd = mad_rows(dif_m);
  }
  else if (method == "mean_abs_dev") {
    dif_m.each_col() -= arma::median(dif_m, 1);
    sd = arma::mean(arma::abs(dif_m), 1);
  }
  else if (method == "biweight_mid") {
    for (arma::uword i = 0; i < sd.n_elem; i++) {
      sd[i] = std::sqrt(bivar(dif_m.row(i).t()));
    }
  }
  else if (method == "std_dev") {
    sd = arma::stddev(dif_m, 0, 1); // normalizing by N-1
  }
  else {
    throw std::invalid_argument(
        "Unrecognised method for estimating standard deviation. Please use one "
        "of the following: median_abs_dev; mean_abs_dev; biweight_mid; "
        "std_dev");
  }
  return sd;
}


void wgt_sam(
  arma::vec& weights,
  const arma::mat& X,
  const arma::uword& nb,
  const arma::uword& nc,
  const arma::uword& hnc
) {

  arma::vec sam_vec(nc, arma::fill::zeros);
  arma::uvec ind_vec = arma::regspace<arma::uvec>(0, X.n_rows - 1);
  arma::umat ind_mat = arma::reshape(ind_vec, nc, nb);

  for (arma::uword j = 0; j < nc; j++) {
    if (j == hnc) {
      continue;
    }
    else {
      arma::vec ngb_sam_vec(nb);  // for storing SAM values of different bands

      for (arma::uword n = 0; n < nb; n++) {
        double s = arma::dot(X.row(ind_mat(hnc, n)), X.row(ind_mat(j, n))) /
                   (arma::norm(X.row(ind_mat(hnc, n))) *
                    arma::norm(X.row(ind_mat(j, n))));

        if (s < -1) {
          s = -1;
        }
        
        if (s > 1) {
          s = 1;
        }
        ngb_sam_vec[n] = std::acos(s);
      }
      sam_vec[j] = arma::sum(ngb_sam_vec);
    }
  }

  arma::vec cell_wgt =
      arma::ones<arma::vec>(nc) - (sam_vec / arma::sum(sam_vec));
  weights = arma::repmat(cell_wgt, nb, 1);
}


arma::vec get_rmse(
  const arma::mat& X,
  const arma::mat& Y,
  const arma::uvec& foc_ind
) {
  return arma::sqrt(
      arma::sum(arma::square(X.rows(foc_ind) - Y.rows(foc_ind)), 1) / X.n_cols);
}


arma::uvec get_ts_period(const arma::uvec& ts_ids) {
  arma::uvec years = ts_ids / 100;
  arma::uvec result(years.n_elem, arma::fill::none);
  arma::uvec unique_years = arma::unique(years);

  // Calculate the period of each unique year and store it in a vector by
  // repeating the values
  for (arma::uword i = 0; i < unique_years.n_elem; ++i) {
    arma::uvec indices = arma::find(years == unique_years(i));
    arma::uword period = indices.n_elem;
    result.elem(indices).fill(period);
  }
  return result;
}


inline arma::uword reflect_index(long long i, arma::uword T) {
  if (T <= 1) return 0;
  
  const long long period = 2LL * static_cast<long long>(T) - 2LL;
  
  i %= period;
  if (i < 0) i += period;
  
  if (i >= static_cast<long long>(T)) {
    i = period - i;
  }
  
  return static_cast<arma::uword>(i);
}


inline arma::uword normalise_filter_width(arma::uword w) {
  if (w == 0) w = 1;
  if (w % 2 == 0) ++w;
  if (w < 3) w = 3;
  return w;
}


arma::mat median_fw_bw_filter(const arma::mat& X, arma::uword w) {
  const arma::uword B = X.n_rows;
  const arma::uword T = X.n_cols;

  arma::mat Y(B, T, arma::fill::zeros);
  if (T == 0) return Y;
  if (w == 0) w = 1;

  arma::mat Yf(B, T, arma::fill::zeros);
  arma::mat Yb(B, T, arma::fill::zeros);

  arma::uvec idx(w);

  // Forward-looking reflected median:
  // nominal window [t, t + w - 1]
  for (arma::uword t = 0; t < T; ++t) {
    for (arma::uword k = 0; k < w; ++k) {
      idx(k) = reflect_index(
        static_cast<long long>(t) + static_cast<long long>(k),
        T
      );
    }

    Yf.col(t) = arma::median(X.cols(idx), 1);
  }

  // Backward-looking reflected median:
  // nominal window [t - w + 1, t]
  for (arma::uword t = 0; t < T; ++t) {
    const long long start =
      static_cast<long long>(t) - static_cast<long long>(w) + 1LL;

    for (arma::uword k = 0; k < w; ++k) {
      idx(k) = reflect_index(start + static_cast<long long>(k), T);
    }

    Yb.col(t) = arma::median(X.cols(idx), 1);
  }

  // Median-of-3 fusion: forward median, backward median, original value
  for (arma::uword b = 0; b < B; ++b) {
    for (arma::uword t = 0; t < T; ++t) {
      double a = Yf(b, t);
      double c = Yb(b, t);
      double x = X(b, t);

      if (a > c) std::swap(a, c);
      if (c > x) std::swap(c, x);
      if (a > c) std::swap(a, c);
      Y(b, t) = c;
    }
  }

  return Y;
}


arma::mat median_fw_bw_filter_adaptive(const arma::mat& X,
                                       const arma::uvec& w_by_col) {
  const arma::uword B = X.n_rows;
  const arma::uword T = X.n_cols;
  
  arma::mat Y(B, T, arma::fill::zeros);
  if (T == 0) return Y;
  if (w_by_col.n_elem != T) {
    throw std::invalid_argument(
        "median_fw_bw_filter_adaptive: w_by_col length must match X.n_cols"
    );
  }
  
  arma::mat Yf(B, T, arma::fill::zeros);
  arma::mat Yb(B, T, arma::fill::zeros);
  
  // Forward-looking reflected median: nominal window [t, t + wt - 1]
  for (arma::uword t = 0; t < T; ++t) {
    arma::uword wt = normalise_filter_width(w_by_col(t));
    arma::uvec idx(wt);
    
    for (arma::uword k = 0; k < wt; ++k) {
      idx(k) = reflect_index(
        static_cast<long long>(t) + static_cast<long long>(k),
        T
      );
    }
    
    Yf.col(t) = arma::median(X.cols(idx), 1);
  }
  
  // Backward-looking reflected median: nominal window [t - wt + 1, t]
  for (arma::uword t = 0; t < T; ++t) {
    arma::uword wt = normalise_filter_width(w_by_col(t));
    arma::uvec idx(wt);
    const long long start =
      static_cast<long long>(t) - static_cast<long long>(wt) + 1LL;
    
    for (arma::uword k = 0; k < wt; ++k) {
      idx(k) = reflect_index(start + static_cast<long long>(k), T);
    }
    
    Yb.col(t) = arma::median(X.cols(idx), 1);
  }
  
  // Median-of-3 fusion: forward median, backward median, original value
  for (arma::uword b = 0; b < B; ++b) {
    for (arma::uword t = 0; t < T; ++t) {
      double a = Yf(b, t);
      double c = Yb(b, t);
      double x = X(b, t);
      
      if (a > c) std::swap(a, c);
      if (c > x) std::swap(c, x);
      if (a > c) std::swap(a, c);
      
      Y(b, t) = c;
    }
  }
  
  return Y;
}


denoise_res denoise_input(
  arma::mat X,
  double noise_zscore,
  const arma::uvec& window_widths
) {

  denoise_res result;

  if (window_widths.n_elem != X.n_cols) {
    throw std::invalid_argument(
        "denoise_input: window_widths length must match X.n_cols");
  }

  arma::mat Y = median_fw_bw_filter_adaptive(X, window_widths);

  arma::mat resid = X - Y;
  arma::vec mad_scale = mad_rows(arma::diff(X, 1, 1)) / std::sqrt(2.0);
  resid.each_col() /= mad_scale;      // row-wise standardization

  arma::mat zscores = arma::abs(resid);
  if (noise_zscore > 0) {
    const arma::uvec outlier_indices = arma::find(zscores >= noise_zscore);
    X.elem(outlier_indices) = Y.elem(outlier_indices);
  }
  result.zscores = zscores;
  result.Y_bln = Y;
  result.X_flt = X;
  return result;
}


arma::mat pairwise_signs_window(
  const arma::mat& Y,
  arma::uword t,
  arma::uword period,
  double prev_cp = arma::datum::nan,
  double next_cp = arma::datum::nan
) {
  arma::uword T = Y.n_cols; // time
  arma::uword B = Y.n_rows; // bands
  if (t >= T) return arma::mat(B, 0); // empty
  
  // pre window bounds (inclusive)
  arma::uword pre_a = (t >= period) ? (t - period) : 0;
  if (std::isfinite(prev_cp))
    pre_a = std::max(pre_a, static_cast<arma::uword>(prev_cp) + 1);
  arma::uword pre_b = (t == 0) ? 0 : (t - 1);
  
  // post window bounds (inclusive)
  arma::uword post_a = t;
  arma::uword post_b = std::min<arma::uword>(T - 1, t + period - 1);
  if (std::isfinite(next_cp))
    post_b = std::min(post_b, static_cast<arma::uword>(next_cp) - 1);

  if (pre_a > pre_b || post_a > post_b) return arma::mat(B, 0);
  
  arma::uword pre_len  = pre_b  - pre_a  + 1;
  arma::uword post_len = post_b - post_a + 1;
  arma::uword L = std::min(pre_len, post_len);
  if (L == 0) return arma::mat(B, 0);
  
  // Align the pair matching adjacent to the changepoint:
  // use the LAST L points of pre and the FIRST L points of post.
  arma::uword pre_start  = pre_b - L + 1;
  arma::uword pre_end    = pre_b;
  arma::uword post_start = post_a;
  arma::uword post_end   = post_a + L - 1;
  
  arma::mat pre  = Y.cols(pre_start,  pre_end);   // B x L
  arma::mat post = Y.cols(post_start, post_end);  // B x L
  
  return arma::sign(post - pre);
}


arma::vec majority_vote_sign(const arma::mat& S) {
  if (S.n_cols == 0) return arma::vec(S.n_rows, arma::fill::zeros);
  arma::vec tally = arma::sum(S, 1);     // B x 1 (sum across columns)
  return arma::sign(tally);              // majority sign per band
}


std::string get_current_timestamp() {
  auto now = std::chrono::system_clock::now();
  auto in_time_t = std::chrono::system_clock::to_time_t(now);
  std::stringstream ss;
  ss << std::put_time(std::localtime(&in_time_t), "%Y%m%d_%H%M%S");
  return ss.str();
}


void save_to_csv(arma::vec data, std::string filename) {
  std::ofstream file;
  file.open(filename);
  
  for (arma::uword i = 0; i < data.n_elem; ++i) {
    file << data[i];
    if (i < data.n_elem - 1) {
      file << ",";
    }
  }
  file << "\n";
  file.close();
}
