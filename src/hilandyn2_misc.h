#ifndef hilandyn_misc_H
#define hilandyn_misc_H

struct denoise_res {
  arma::mat Y_bln;
  arma::mat X_flt;
  arma::mat zscores;
};

arma::umat rle(const arma::uvec& x);

arma::mat interpolate(arma::cube& XY, arma::mat& XI);

arma::mat extrapolate(arma::cube& XY, arma::mat& XI);

void fill_seasonal_nan(
  arma::mat &Y,
  const arma::uword period,
  const arma::uword min_obs
);

double mad(const arma::vec& x);

arma::vec mad_rows(
  const arma::mat &X,
  double consistency_c = 1.4826,
  double eps = 1e-12
);

double bivar(const arma::vec& x);

arma::vec sd_est(const arma::mat& X, const std::string& method);

void wgt_sam(
  arma::vec &weights,
  const arma::mat &X,
  const arma::uword &nb,
  const arma::uword &nc,
  const arma::uword &hnc
);

arma::vec get_rmse(
  const arma::mat &X,
  const arma::mat &Y,
  const arma::uvec &foc_ind
);

arma::uvec get_ts_period(const arma::uvec& ts_ids);

denoise_res denoise_input(
  arma::mat X,
  double noise_zscore,
  const arma::uvec& window_widths
);

arma::mat pairwise_signs_window(
  const arma::mat& Y,
  arma::uword t,
  arma::uword period,
  double prev_cp,
  double next_cp
);

arma::vec majority_vote_sign(const arma::mat& S);

std::string get_current_timestamp();

void save_to_csv(arma::vec data, std::string filename);

#endif
