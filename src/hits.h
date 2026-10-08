#ifndef hits_H
#define hits_H

struct detail_res {
  arma::vec details;
  arma::mat detail_all;
};

struct CandidateCache {
  double th;
  arma::mat Y;
  arma::uvec cpt;  
  double rss;
  double bic;
};

arma::uvec hd_bts_dcmp(
  arma::cube& decomp_hist,
  arma::mat& bts_coeffs,
  const arma::vec& w,
  const double p,
  const arma::uword d,
  const arma::uword n
);

detail_res compute_det(arma::cube& decomp_hist);

arma::cube hd_bts_dns(
  arma::cube decomp_hist,
  const arma::uvec& sameboat,
  const detail_res& details_out,
  const arma::uword& n,
  const double& lambda,
  const double& bal,
  const double& seg_len_min
);

arma::mat hd_bts_inv(
  arma::cube& decomp_hist,
  arma::mat bts_coeffs,
  const arma::uword& n
);

arma::mat hd_bts_cpt(
  arma::mat X,
  const arma::vec& sd,
  arma::uvec& cpt,
  double& th_selected,
  arma::vec& weights,
  const double& th_max,
  const double& th_min,
  const double& dBIC_min,
  double& seg_len_min,
  const arma::uword cpt_max,
  bool verbose
);
#endif
