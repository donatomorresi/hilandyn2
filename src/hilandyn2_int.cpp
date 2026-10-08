#include <RcppArmadillo.h>
#include <chrono>
#include "hits.h"
#include "hits_misc.h"
#include "hilandyn2_misc.h"
#include "hilandyn2_deseason.h"
#ifdef _OPENMP
#include <omp.h>
#endif

				   
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]


// [[Rcpp::export]]
arma::vec hilandyn2_int_cpp(
  arma::mat& X,
  arma::uword nb,
  arma::uword nc,
  const arma::uword nt,
  const arma::uword nr,
  const arma::uvec ts_ids,
  const bool cell_weights,
  const bool rmse_out,
  const bool len_seg_out,
  const bool slo_seg_out,
  const bool bln_mat_out,
  const bool zsc_mat_out,
  const bool mag_mat_out,
  const bool sgn_mat_out,
  const bool est_mat_out,
  const arma::uword nout,
  const arma::uword ts_len_min,
  const arma::uword gap_len_max,
  const std::string sd_method,
  const arma::vec cng_dir,
  const double noise_zscore,
  const std::vector<std::string>& wavelet_filters,
  const arma::sword wavelet_level,
  const double wavelet_filter_penalty,
  const double acf_min,
  double th_min,
  double th_max,
  double dBIC_min,
  double seg_len_min,
  const arma::uword cpt_max,
  bool verbose
) {

  // min length check
  if (X.n_cols < ts_len_min) {
    return arma::vec(nout, arma::fill::value(arma::datum::nan));
  }

  // zeroes check
  X.replace(0.0, arma::datum::nan);

  // nodata check
  arma::uvec nf_ind = arma::find_nonfinite(X);
  
  if (nf_ind.n_elem == X.n_elem) {
    return arma::vec(nout, arma::fill::value(arma::datum::nan));
  }
  
  // positional indices
  const arma::uword hnb = nb == 1 ? 1 : nb / 2;
  arma::uword hnc = (nc - 1) / 2;
  arma::uvec row_ind = arma::regspace<arma::uvec>(0, nr - 1);
  arma::umat ind_mat = arma::reshape(row_ind, nc, nb);
  arma::uvec foc_ind = ind_mat.row(hnc).t();

  // nodata processing
  arma::umat X_na(nr, nt, arma::fill::zeros);
  X_na(nf_ind).ones();
  arma::uvec na_row_cnt =
      arma::sum(X_na, 0).t(); // count NA rows per column (length equal to nt)

  arma::uword cc_num = nt;   // number of complete columns
  arma::uvec cc_ind = arma::regspace<arma::uvec>(0, nt - 1);
  arma::uvec ec_ind = cc_ind;
  arma::uvec na_col_cnd(nt);
  arma::umat gap_mat;
  
  if (any(na_row_cnt > 0)) {
    
    // Stage 1: spatial gap-fill (process gaps in columns)
    bool filled = false;
    arma::uvec na_row_ind = arma::find(na_row_cnt > 0);
    arma::vec na_row_perc =
        arma::conv_to<arma::vec>::from(na_row_cnt(na_row_ind)) / nr;

    // iterate over columns (time points)
    arma::uvec submat_ind;
    arma::vec Xi;
    for (arma::uword i = 0; i < na_row_perc.n_elem; ++i) {
      arma::uword col_i = na_row_ind(i);
      Xi = X.col(col_i);

      // fill values using summary value across cells
      if (na_row_perc[i] < 0.5) {
        
        // iterate over bands
        arma::vec Xi_b;
        arma::uvec finite_ind;
        for (arma::uword b = 0; b < nb; ++b) {
          arma::uvec band_ind = ind_mat.col(b);
          Xi_b = Xi.rows(band_ind);
          finite_ind = arma::find_finite(Xi_b);
          
          if (finite_ind.is_empty()) continue;

          submat_ind = band_ind + nr * col_i;
          submat_ind.shed_rows(finite_ind);
          double fill_value = arma::mean(Xi_b(finite_ind));
          X(submat_ind) =
              arma::vec(submat_ind.n_elem, arma::fill::value(fill_value));
          filled = true;
        }
      }
    }

    // Stage 2: check columns with missing values
    if (filled) {
      nf_ind = arma::find_nonfinite(X);
      X_na.zeros();
      X_na(nf_ind).ones();
      na_row_cnt = arma::sum(X_na, 0).t();
    }
    
    na_col_cnd = na_row_cnt > 0;
    cc_num = arma::sum(na_col_cnd == 0);
    cc_ind = arma::find(na_col_cnd == 0);
    ec_ind = arma::find(na_col_cnd == 1);
    
    gap_mat = rle(na_col_cnd);
    gap_mat.shed_rows(arma::find(gap_mat.col(0) == 0));
    
    if (any(gap_mat.col(1) > gap_len_max)) {
      return arma::vec(nout, arma::fill::value(arma::datum::nan));
    }
  }

  // Get time series period
  arma::uvec ts_period = get_ts_period(ts_ids);
  arma::uword period = arma::max(ts_period);
  
  // Create empty outputs
  arma::mat Y_bln;
  arma::mat Y_bln_foc;
  arma::mat zscores;
  arma::mat zscores_foc;
  
  if (bln_mat_out) {
    Y_bln_foc.set_size(nb, nt);
    Y_bln_foc.fill(arma::datum::nan);
  }
  
  if (zsc_mat_out) {
    zscores.set_size(arma::size(X));
    zscores.fill(arma::datum::nan);
    zscores_foc.set_size(nb, nt);
    zscores_foc.fill(arma::datum::nan);
  }

  // Temporal noise filter
  if (noise_zscore > 0 || cc_num < nt) {
    Y_bln.set_size(arma::size(X));
    Y_bln.fill(arma::datum::nan);
    
    arma::uvec bln_win_width(cc_num, arma::fill::value(period));
    arma::uvec years = ts_ids / 100;
    arma::uvec cc_years = years.elem(cc_ind);
    arma::uvec unique_years = arma::unique(years);
    arma::uvec year_complete_cnt(unique_years.n_elem, arma::fill::zeros);
    
    for (arma::uword i = 0; i < unique_years.n_elem; ++i) {
      arma::uvec year_ind = arma::find(years == unique_years(i));
      year_complete_cnt(i) = arma::sum(na_col_cnd.elem(year_ind) == 0);
    }

    const arma::uword wide_win_width = std::min(period + period / 2, cc_num);

    for (arma::uword i = 0; i < unique_years.n_elem; ++i) {
      arma::uvec year_cc_ind = arma::find(cc_years == unique_years(i));
      
      if (year_complete_cnt(i) < period && !year_cc_ind.is_empty()) {
        bln_win_width.elem(year_cc_ind).fill(wide_win_width);
      }
    }
    
    denoise_res denoise_out = denoise_input(X.cols(cc_ind), 
                                            noise_zscore, bln_win_width);
    // Denoising
    if (noise_zscore > 0) {
      X.cols(cc_ind) = denoise_out.X_flt;
    }
    
    Y_bln.cols(cc_ind) = denoise_out.Y_bln;

    // Gap-fill
    if (cc_num < nt) {
      fill_seasonal_nan(X, period, 2);
    }
  
    // Output baseline
    if (bln_mat_out) {
      Y_bln_foc = Y_bln.rows(foc_ind);
    }
    
    // Output zscores
    if (zsc_mat_out) {
      zscores.cols(cc_ind) = denoise_out.zscores;
      zscores_foc = zscores.rows(foc_ind);
    }
  }

  // Seasonal adjustment
  double s_orig = arma::datum::nan;
  double s_adj = arma::datum::nan;
  
  if (wavelet_level != 0) {
    MODWPTDeseasonResult deseason_out =
        deseasonalise_modwpt(X, period, wavelet_filters, wavelet_level,
                             wavelet_filter_penalty, acf_min);

    s_orig = deseason_out.acf_original;
    s_adj = deseason_out.acf_adjusted;

    if (verbose) {
      std::cout << std::boolalpha
                << "Used adjusted signal: " << deseason_out.used_adjusted
                << "\nOriginal ACF score=" << s_orig
                << "\nAdjusted ACF score=" << s_adj 
                << "\n"
                << std::endl;
    }

    X = deseason_out.X;
  }
  
  // Compute cell-wise weights
  arma::vec weights(nr, arma::fill::ones);
  if (cell_weights) {
    wgt_sam(weights, X, nb, nc, hnc);
  }

  // Estimate standard deviation
  arma::vec sd = sd_est(X, sd_method);

  if (any(sd == 0) || !sd.is_finite()) {
    return arma::vec(nout, arma::fill::value(arma::datum::nan));
  }
  
  // Detect changepoints and estimate signal
  arma::uvec cpt;
  double th_selected = arma::datum::nan;
  if (seg_len_min == -1) {
    seg_len_min = static_cast<double>(period);
  }
  arma::mat Y = hd_bts_cpt(X, sd, cpt, th_selected, weights, 
                           th_max, th_min, dBIC_min,
                           seg_len_min, cpt_max, verbose);
  
  arma::mat X_foc = X.rows(foc_ind);
  arma::mat Y_foc = Y.rows(foc_ind);
  // Default outputs
  arma::vec cpt_id(nt, arma::fill::value(arma::datum::nan));
  arma::vec rel_mag_med(nt, arma::fill::value(arma::datum::nan));

  // Optional outputs
  arma::vec rmse_vec;
  arma::vec len_seg;
  arma::mat slo_seg;
  arma::mat mag_foc;

  // Compute RMSE
  if (rmse_out) {
    rmse_vec = get_rmse(X, Y, foc_ind);
  }
  
  // Segment stats
  if (len_seg_out || slo_seg_out) {
    arma::uvec seg_beg = arma::join_vert(arma::zeros<arma::uvec>(1), cpt);
    
    if (len_seg_out) {
      len_seg.set_size(nt);
      len_seg.fill(arma::datum::nan);
      for (arma::uword i = 0; i < seg_beg.n_elem; ++i) {
        if (seg_beg.n_elem == 1) {
          len_seg[seg_beg[i]] = nt;
        } else if (i == seg_beg.n_elem - 1) {
          len_seg[seg_beg[i]] = nt - seg_beg[i];
        } else {
          len_seg[seg_beg[i]] = seg_beg[i + 1] - seg_beg[i];
        }
      }
    }
    
    if (slo_seg_out) {
      slo_seg.set_size(nb, nt);
      slo_seg.fill(arma::datum::nan);
      for (arma::uword i = 0; i < seg_beg.n_elem; ++i) {
        const arma::uword s0 = seg_beg[i];
        if (s0 + 1 < nt) {
          slo_seg.col(s0) = Y_foc.col(s0 + 1) - Y_foc.col(s0);
        } else {
          slo_seg.col(s0).zeros();
        }
      }
    }
  }
  
  // Output magnitude matrix
  if (mag_mat_out) {
    mag_foc.set_size(nb, nt);
    mag_foc.fill(arma::datum::nan);
  }
  
  // Change magnitude and type
  double d_max_ti = arma::datum::nan, d_max_mg = arma::datum::nan,
         d_frs_ti = arma::datum::nan, d_frs_mg = arma::datum::nan,
         d_lst_ti = arma::datum::nan, d_lst_mg = arma::datum::nan,
         d_num = arma::datum::nan, g_max_ti = arma::datum::nan,
         g_max_mg = arma::datum::nan, g_frs_ti = arma::datum::nan,
         g_frs_mg = arma::datum::nan, g_lst_ti = arma::datum::nan,
         g_lst_mg = arma::datum::nan, g_num = arma::datum::nan;

  if (!cpt.is_empty()) {
    arma::mat mag1_est(nc, nb);
    arma::umat cnd_d_mat(nc, nb), cnd_g_mat(nc, nb);
    arma::mat cng_sgn_mat;
    arma::vec cng_dir_cpti(nr);
    arma::vec rel_mag_ci(nr);
    arma::vec mag_foc_cpti(nb);
    arma::vec std_mag_foc_cpti(nb);
    arma::vec glb_mag_vec(nt);
    double prev_cpt, next_cpt;

    // Compute per-band std dev on fitted values
    arma::vec sd_foc = arma::stddev(Y_foc, 0, 1);

    for (arma::uword i = 0; i < cpt.n_elem; i++) {
      const arma::uword cpti = cpt[i];
      
      const arma::uword ts_id_cpti = ts_ids(cpti);
      arma::uvec cpti_pre_ind = arma::find(ts_ids == ts_id_cpti - 100);
      if (cpti_pre_ind.is_empty()) cpti_pre_ind = {cpti > 0 ? cpti - 1 : cpti};
      const arma::uword cpti_pre = cpti_pre_ind[0];
      mag1_est(row_ind) = Y.col(cpti_pre) - Y.col(cpti);
      
      mag_foc_cpti = mag1_est.row(hnc).t();
      std_mag_foc_cpti = mag_foc_cpti / sd_foc;
      glb_mag_vec[cpti] = arma::norm(std_mag_foc_cpti, 2);
      
      if (mag_mat_out) {
        mag_foc.col(cpti) = mag_foc_cpti; 
      }
      
      // Compute direction of change using pre/post pairs
      prev_cpt = cpt.in_range(i - 1) ? static_cast<double>(cpt(i - 1))
                                     : arma::datum::nan;
      next_cpt = cpt.in_range(i + 1) ? static_cast<double>(cpt(i + 1))
                                     : arma::datum::nan;
      cng_sgn_mat = pairwise_signs_window(X, cpti, period, prev_cpt, next_cpt);
      cng_dir_cpti = majority_vote_sign(cng_sgn_mat);
      
      cnd_d_mat(row_ind) = cng_dir_cpti == cng_dir;
      cnd_g_mat(row_ind) = cng_dir_cpti == -cng_dir;

      arma::uvec cnd_d = arma::sum(cnd_d_mat, 1);
      arma::uvec cnd_g = arma::sum(cnd_g_mat, 1);
      
      rel_mag_ci = (mag1_est.as_col() / Y.col(cpti_pre)) * 100.0;
      arma::uvec fin = arma::find_finite(rel_mag_ci);
      
      if (!fin.is_empty()) {
        rel_mag_med[cpti] = arma::median(arma::abs(rel_mag_ci(fin)));
      }
      
      if (cnd_d[hnc] >= hnb) {
        cpt_id[cpti] = 101;    // disturbance
      }
      else if (cnd_g[hnc] >= hnb) {
        cpt_id[cpti] = 102;    // greening
      }
      else {
        cpt_id[cpti] = 9;
      }
    }

    // summaries
    arma::uvec d_ind = arma::find(cpt_id == 101);
    arma::uvec g_ind = arma::find(cpt_id == 102);

    if (!d_ind.is_empty()) {
      arma::uword d_max_ind = d_ind[arma::index_max(glb_mag_vec(d_ind))];
      d_max_ti = ts_ids[d_max_ind];
      d_max_mg = rel_mag_med[d_max_ind];
      d_frs_ti = ts_ids[d_ind[0]];
      d_frs_mg = rel_mag_med[d_ind[0]];
      d_lst_ti = ts_ids[d_ind[d_ind.n_elem - 1]];
      d_lst_mg = rel_mag_med[d_ind[d_ind.n_elem - 1]];
      d_num = d_ind.n_elem;
    }

    if (!g_ind.is_empty()) {
      arma::uword g_max_ind = g_ind[arma::index_max(glb_mag_vec(g_ind))];
      g_max_ti = ts_ids[g_max_ind];
      g_max_mg = rel_mag_med[g_max_ind];
      g_frs_ti = ts_ids[g_ind[0]];
      g_frs_mg = rel_mag_med[g_ind[0]];
      g_lst_ti = ts_ids[g_ind[g_ind.n_elem - 1]];
      g_lst_mg = rel_mag_med[g_ind[g_ind.n_elem - 1]];
      g_num = g_ind.n_elem;
    }
  }

  arma::vec sng_val = {
    d_max_ti, d_max_mg, d_frs_ti, d_frs_mg, d_lst_ti, d_lst_mg, d_num,
    g_max_ti, g_max_mg, g_frs_ti, g_frs_mg, g_lst_ti, g_lst_mg, g_num,
    s_orig, s_adj, th_selected
  };
  
  // Allocate output
  arma::vec out(nout, arma::fill::value(arma::datum::nan));
  arma::uword pos = 0;
  
  // Append outputs
  out.subvec(pos, pos + sng_val.n_elem - 1) = sng_val; 
  pos += sng_val.n_elem;
  
  // nc vector
  if (rmse_out) {
    out.subvec(pos, pos + rmse_vec.n_elem - 1) = rmse_vec;
    pos += rmse_vec.n_elem;
  }
  
  // nt vectors
  out.subvec(pos, pos + cpt_id.n_elem - 1) = cpt_id;
  pos += cpt_id.n_elem;
  
  out.subvec(pos, pos + rel_mag_med.n_elem - 1) = rel_mag_med;
  pos += rel_mag_med.n_elem;
  
  if (len_seg_out) {
    out.subvec(pos, pos + len_seg.n_elem - 1) = len_seg;
    pos += len_seg.n_elem;
  }
  
  // nc * nt matrices
  if (slo_seg_out) {
    out.subvec(pos, pos + slo_seg.n_elem - 1) = slo_seg.as_col();
    pos += slo_seg.n_elem;
  }
  
  if (bln_mat_out) {
    out.subvec(pos, pos + Y_bln_foc.n_elem - 1) = Y_bln_foc.as_col();
    pos += Y_bln_foc.n_elem;
  }
  
  if (zsc_mat_out) {
    out.subvec(pos, pos + zscores_foc.n_elem - 1) = zscores_foc.as_col();
    pos += zscores_foc.n_elem;
  }
  
  if (mag_mat_out) {
    out.subvec(pos, pos + mag_foc.n_elem - 1) = mag_foc.as_col();
    pos += mag_foc.n_elem;
  }
  
  if (sgn_mat_out) {
    out.subvec(pos, pos + X_foc.n_elem - 1) = X_foc.as_col();
    pos += X_foc.n_elem;
  }
  
  if (est_mat_out) {
    out.subvec(pos, pos + Y_foc.n_elem - 1) = Y_foc.as_col();
    pos += Y_foc.n_elem;
  }
  
  return out;
}
 

// [[Rcpp::export]]
arma::mat hilandyn2_map_cpp(
  arma::vec& v,
  const arma::uvec& v_dim,
  const arma::uword& win_side,
  const arma::uword& nb,
  arma::uword& nc,
  const arma::uword& nt,
  const arma::uword& nr,
  const arma::uvec& ts_ids,
  const bool cell_weights,
  const bool rmse_out,
  const bool len_seg_out,
  const bool slo_seg_out,
  const bool bln_mat_out,
  const bool zsc_mat_out,
  const bool mag_mat_out,
  const bool sgn_mat_out,
  const bool est_mat_out,
  const arma::uword& nout,
  const arma::uword ts_len_min,
  const arma::uword gap_len_max,
  const std::string& sd_method,
  const arma::vec& cng_dir,
  const double noise_zscore,
  const std::vector<std::string>& wavelet_filters,
  const arma::sword wavelet_level,
  const double wavelet_filter_penalty,
  const double acf_min,
  double th_min,
  double th_max,
  double dBIC_min,
  double seg_len_min,
  const arma::uword cpt_max,
  const int num_threads,
  bool debug
) {

  if (v_dim.n_elem != 3) {
    Rcpp::stop("mapping input dimensions must contain rows, columns, layers");
  }
  const arma::uword input_rows = v_dim(0);
  const arma::uword input_columns = v_dim(1);
  const arma::uword input_layers = v_dim(2);
  if (win_side < 1 || win_side % 2 == 0 ||
      input_rows < win_side || input_columns < win_side) {
    Rcpp::stop("invalid focal-window dimensions");
  }
  const arma::uword input_cells = input_rows * input_columns;
  if (v.n_elem != input_cells * input_layers) {
    Rcpp::stop("mapping input values do not match their dimensions");
  }

  const arma::uword half_window = (win_side - 1) / 2;
  const arma::uword window_cells = win_side * win_side;
  const arma::uword output_rows = input_rows - 2 * half_window;
  const arma::uword output_columns = input_columns - 2 * half_window;
  const arma::uword output_cells = output_rows * output_columns;
  const arma::uword focal_elements = window_cells * input_layers;
  if (nc != window_cells || nr != nb * nc ||
      input_layers != nb * nt || focal_elements != nr * nt) {
    Rcpp::stop("mapping dimensions are inconsistent");
  }

  arma::umat window_grid(win_side, win_side, arma::fill::none);
  window_grid.each_col() =
      arma::regspace<arma::uvec>(0, win_side - 1);
  window_grid.each_row() +=
      arma::regspace<arma::urowvec>(0, win_side - 1) * input_columns;

  arma::umat layer_grid(window_cells, input_layers, arma::fill::none);
  layer_grid.each_col() = arma::vectorise(window_grid);
  layer_grid.each_row() +=
      arma::regspace<arma::urowvec>(0, input_layers - 1) * input_cells;
  const arma::uvec sample_offsets = arma::vectorise(layer_grid);

  arma::mat Y(nout, output_cells, arma::fill::value(arma::datum::nan));
  
  if (nt < ts_len_min) return Y.t();

  arma::uword error_count = 0;
  std::string first_error;
  const std::string debug_prefix = debug
      ? "input_data_" + get_current_timestamp() + "_" +
        std::to_string(
          std::chrono::steady_clock::now().time_since_epoch().count())
      : "";

#ifndef _OPENMP
  (void)num_threads;
#endif
  
#ifdef _OPENMP
  #pragma omp parallel num_threads(num_threads)
#endif

  {
    arma::mat Xi(nr, nt);
    arma::uvec sample_indices(sample_offsets.n_elem, arma::fill::none);

    arma::uword local_error_count = 0;
    std::string local_first_error;

#ifdef _OPENMP
    #pragma omp for schedule(auto)
#endif

    for (arma::uword i = 0; i < output_cells; i++) {
      const arma::uword output_cell = i;
      const arma::uword output_row = output_cell / output_columns;
      const arma::uword output_column = output_cell % output_columns;
      const arma::uword anchor =
          output_row * input_columns + output_column;
      sample_indices = sample_offsets;
      sample_indices += anchor;
      arma::vec Xi_values(Xi.memptr(), Xi.n_elem, false, true);
      Xi_values = v.elem(sample_indices);

      try {
        Y.col(i) = hilandyn2_int_cpp(
            Xi, nb, nc, nt, nr, ts_ids, cell_weights, rmse_out, len_seg_out,
            slo_seg_out, bln_mat_out, zsc_mat_out, mag_mat_out, sgn_mat_out,
            est_mat_out, nout, ts_len_min, gap_len_max, sd_method, cng_dir,
            noise_zscore, wavelet_filters, wavelet_level,
            wavelet_filter_penalty, acf_min, th_min, th_max, dBIC_min,
            seg_len_min, cpt_max, /*verbose*/ false);
      } 
      catch (const std::exception& e) {
        ++local_error_count;
        if (local_first_error.empty()) {
          local_first_error = std::string(e.what()).substr(0, 512);
        }

        if (debug) {
          const std::string filename = debug_prefix + "_cell_" +
              std::to_string(output_cell + 1) + ".csv";
          save_to_csv(Xi.as_col(), filename);
        }
      }
    }

#ifdef _OPENMP
    #pragma omp critical(hilandyn2_map_errors)
#endif
    {
      error_count += local_error_count;
      if (first_error.empty() && !local_first_error.empty()) {
        first_error = local_first_error;
      }
    }
  }

  if (error_count > 0) {
    Rcpp::warning(
      std::to_string(static_cast<unsigned long long>(error_count)) +
      " of " +
      std::to_string(static_cast<unsigned long long>(output_cells)) +
      " mapping cells failed; their outputs are missing. First error: " +
      first_error
    );
  }

  return Y.t();
}
