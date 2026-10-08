#ifndef hilandyn2_deseason_H
#define hilandyn2_deseason_H

#include <RcppArmadillo.h>
#include <string>
#include <vector>

struct MODWPTree {
  std::vector<std::vector<arma::mat>> levels;
};

struct MODWPTDeseasonResult {
  arma::mat X;
  double acf_original;
  double acf_adjusted;
  double acf_gain;
  bool used_adjusted;
  arma::imat seasonal_nodes;
  std::string selected_filter;
  double filter_width_penalty;
  double objective_score;
  bool skipped_low_acf;
};

double seasonal_structure_score(const arma::mat& X, arma::uword period);

MODWPTDeseasonResult
deseasonalise_modwpt(
  arma::mat &X,
  arma::uword period,
  const std::vector<std::string> &filter_names,
  int J,
  double penalty_lambda = 0.05,
  double acf_min = 0.1,
  double dt = 1
);

#endif
