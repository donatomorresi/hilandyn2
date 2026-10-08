#ifndef hits_misc_H
#define hits_misc_H

struct upd_res {
  arma::vec weights_const;
  arma::vec weights_lin;
  arma::rowvec bts_coeffs;
  arma::vec h;
  arma::vec tc1;
  arma::uvec idx;
};

arma::mat wgt_det(
  arma::mat x,
  arma::umat ind_mat,
  arma::uword nc,
  arma::uword nb
);

arma::vec match(const arma::uvec& a, const arma::uvec& b);

arma::vec filter_bts(const arma::vec& a);

arma::mat filter_bts_mat(const arma::mat& A);

arma::mat orth_matrix(const arma::vec& d);

upd_res updating(
  const arma::umat &ee,
  arma::vec wgt_const,
  arma::vec wgt_lin,
  arma::rowvec bts_c,
  const arma::uvec &idx
);

arma::mat balance_np(
  const arma::uvec &paired,
  const arma::umat &ee,
  const arma::uvec &idx,
  const arma::uword &no_of_current_steps,
  const arma::uword &n
);

arma::mat balance_p(
  const arma::umat &pr,
  const arma::umat &ee_p1,
  const arma::umat &ee_p2,
  const arma::uvec &idx,
  const arma::uword &n
);

arma::uvec finding_cp(
  const arma::cube &decomp_hist,
  const arma::uvec &sameboat,
  const arma::uword &n
);

#endif
