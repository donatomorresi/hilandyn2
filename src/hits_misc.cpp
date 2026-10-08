#include <RcppArmadillo.h>
#include "hits_misc.h"


arma::vec match(const arma::uvec& a, const arma::uvec& b) {

  arma::vec res(a.n_elem, arma::fill::value(arma::datum::nan));

  for (arma::uword i = 0; i < a.n_elem; ++i) {
    arma::uvec logi_vec = (b == a[i]);

    if (any(logi_vec)) {
      arma::uvec ind = arma::find(logi_vec);
      res[i] = arma::as_scalar(ind[0]);
    }
  }
  return res;
}


arma::vec filter_bts(const arma::vec& a) {

  double d1 = a[4] * a[0] - a[3] * a[1];  // A zero value causes NaNs
  double d2 = a[1] * a[5] - a[2] * a[4];
  double d3 = a[3] * a[2] - a[0] * a[5];

  if (d1 != 0.0) {

    double d1_sq = d1 * d1;
    double denom = d1_sq + (d2 * d2) + (d3 * d3);
    
    if (denom == 0.0) {
      return arma::zeros<arma::vec>(3);
    }
    
    double w = -std::sqrt(d1_sq / denom);
    double u = w * d2 / d1;
    double v = w * d3 / d1;
    
    return {u, v, w};
    
  } else {
    
    if (d2 == 0.0) {
      return arma::zeros<arma::vec>(3);
    }
    
    double d2_sq = d2 * d2;
    double denom = d2_sq + (d3 * d3);
    
    if (denom == 0.0) {
      return arma::zeros<arma::vec>(3);
    }
    
    double w_f = -std::sqrt(d2_sq / denom);
    double u_f = 0.0;
    double v_f = w_f * d3 / d2;

    return {w_f, v_f, u_f};
  }
}


arma::mat filter_bts_mat(const arma::mat& A) {

  arma::vec a5a1_a4a2_dif =
      A.col(4) % A.col(0) -
      A.col(3) % A.col(1); // NaNs are caused when elements here are zero
  arma::uvec zero_indices = arma::find(a5a1_a4a2_dif == 0);
  arma::mat df(A.n_rows, 3, arma::fill::none);
  
  if (!zero_indices.is_empty()) {
    
    arma::mat A2 = arma::flipud(A);
    
    arma::vec d1_f = A2.col(4) % A2.col(0) - A2.col(3) % A2.col(1);
    arma::vec d2_f = A2.col(1) % A2.col(5) - A2.col(2) % A2.col(4);
    arma::vec d3_f = A2.col(3) % A2.col(2) - A2.col(0) % A2.col(5);
    
    arma::vec d1_f_sq = d1_f % d1_f;
    arma::vec denom = d1_f_sq + (d2_f % d2_f) + (d3_f % d3_f);
    
    arma::vec w_f = -arma::sqrt(d1_f_sq / denom);

    df.col(0) = w_f;
    df.col(1) = w_f % d3_f / d1_f;  //v
    df.col(2) = w_f % d2_f / d1_f;  //u
    
  } else {
    arma::vec a2a6_a3a5_dif = A.col(1) % A.col(5) - A.col(2) % A.col(4);
    arma::vec a4a3_a1a6_dif = A.col(3) % A.col(2) - A.col(0) % A.col(5);

    arma::vec d1_sq = a5a1_a4a2_dif % a5a1_a4a2_dif;
    arma::vec denom = d1_sq + (a2a6_a3a5_dif % a2a6_a3a5_dif) +
                      (a4a3_a1a6_dif % a4a3_a1a6_dif);

    arma::vec w = -arma::sqrt(d1_sq / denom);
    
    df.col(0) = w % a2a6_a3a5_dif / a5a1_a4a2_dif;  //u
    df.col(1) = w % a4a3_a1a6_dif / a5a1_a4a2_dif;  //v
    df.col(2) = w;
  }
  
  return df;
}


arma::mat orth_matrix(const arma::vec& d) {

  arma::mat33 M;
  M.col(0) = d / arma::norm(d);

  double u = M[0], v = M[1], w = M[2];

  arma::vec v1 = {1 - u*u, -u*v, -u*w};
  arma::vec v2 = {0, -w, v};

  M.col(1) = v1 / arma::norm(v1);
  M.col(2) = v2 / arma::norm(v2);

  return M;
}

upd_res updating(
  const arma::umat &ee,
  arma::vec wgt_const,
  arma::vec wgt_lin,
  arma::rowvec bts_c,
  const arma::uvec &idx
) {

  arma::vec wc0 = wgt_const(ee);
  arma::vec wl0 = wgt_lin(ee);
  arma::vec tc0 = bts_c(ee);

  arma::vec m = {wc0[0], wc0[1], wc0[2], wl0[0], wl0[1], wl0[2]};

  arma::vec h = filter_bts(m);
  arma::mat33 M0 = orth_matrix(h).t();
  arma::vec tc1 = M0 * tc0;

  const arma::uvec ee_ind = {ee[2], ee[0], ee[1]};

  // Update the original vectors
  bts_c(ee_ind) = tc1.t();
  wgt_const(ee_ind) = M0 * wc0;
  wgt_lin(ee_ind) = M0 * wl0;
  
  // Create a copy of idx and modify it
  arma::uvec idx1 = idx;
  arma::uvec found_indices = arma::find(idx1 == ee_ind[0], 1, "first");
  if (!found_indices.is_empty()) {
    idx1.shed_row(found_indices(0));
  }

  // Prepare the struct to be returned
  upd_res result;
  result.weights_const = wgt_const;
  result.weights_lin = wgt_lin;
  result.bts_coeffs = bts_c;
  result.h = h;
  result.tc1 = tc1;
  result.idx = idx1;

  return result;
}

arma::mat balance_np(
  const arma::uvec &paired,
  const arma::umat &ee,
  const arma::uvec &idx,
  const arma::uword &no_of_current_steps,
  const arma::uword &n
) {

  arma::umat prd(ee.n_rows, 3);
  
  arma::vec vm = match(vectorise(ee), paired);
  prd(find_nonfinite(vm)).fill(0);
  prd(find_finite(vm)).fill(1);
  
  arma::mat blnc(3, no_of_current_steps);
  arma::uvec f2 = arma::find(arma::sum(prd.cols(0,1), 1) == 2);
  arma::uvec l2 = arma::find(arma::sum(prd.cols(1,2), 1) == 2);
  arma::uvec f2_l2 = arma::join_vert(f2, l2);
  
  arma::uvec nopair = arma::regspace<arma::uvec>(0, ee.n_rows - 1);
  
  static const arma::uvec col1 = {0};
  static const arma::uvec col2 = {1};
  static const arma::uvec col3 = {2};

  arma::vec ee_l2_col2 = arma::conv_to<arma::vec>::from(ee(l2, col2));
  
  if (f2_l2.n_elem > 0) {
    nopair.shed_rows(f2_l2);
  } 
  
  if (f2.n_elem > 0) {

    arma::vec prtn =
        arma::conv_to<arma::vec>::from(ee(f2, col3) - ee(f2, col1));

    arma::vec prtn3 = prtn + 1;
    arma::vec prtn1 = prtn / prtn3;
    arma::vec prtn2 = 1 / prtn3;

    blnc.cols(f2) = arma::join_vert(prtn1, prtn2, prtn3);
  }
  
  if (l2.n_elem > 0) {
    
    arma::vec vm = match(ee(l2, col3), idx);
    arma::vec prtn(vm.n_elem);
    
    for (arma::uword i = 0; i < vm.n_elem; ++i) {
      
      if (std::isnan(vm(i))) {
        double n_ee_na = n - ee(i, 2) + 1;
        prtn(i) = n_ee_na;
      }
      else {
        arma::uword vmi = vm(i);
        
        if (idx.in_range(vmi + 1)) {    
          double idx_vmi = idx(vmi + 1);
          prtn(i) = idx_vmi - ee_l2_col2(i);
        }
        else {
          double n_ee_na = n - ee(i, 2) + 1;
          prtn(i) = n_ee_na;
        }
      }
    }
    
    arma::vec prtn3 = prtn + 1;
    arma::vec prtn1 = 1 / prtn3;
    arma::vec prtn2 = prtn / prtn3;
    blnc.cols(l2) = arma::join_vert(prtn1, prtn2, prtn3);
  }
  
  if (!nopair.is_empty()) {
    double onethird = 1;
    onethird /= 3;
    blnc.cols(nopair).fill(onethird);
  }
  return blnc;
}

arma::mat balance_p(
  const arma::umat &pr,
  const arma::umat &ee_p1,
  const arma::umat &ee_p2,
  const arma::uvec &idx,
  const arma::uword &n
) {

  arma::mat blnc(3, pr.n_cols, arma::fill::zeros);

  arma::vec c1 = arma::conv_to<arma::vec>::from(ee_p1.col(2) - ee_p1.col(0));
  arma::vec c2(ee_p2.n_rows);
  
  arma::vec vm = match(ee_p2.col(2), idx);
  
  for (arma::uword i = 0; i < ee_p2.n_rows; ++i) {
    
    if (std::isnan(vm(i))) {
      
      double n_ee_p1 = n - ee_p1(i, 2) + 1.0;
      c2(i) = n_ee_p1;
    }
    else {
      arma::uword vmi = arma::as_scalar(vm(i));
      
      if (idx.in_range(vmi + 1)) {        
        double idx_vmi = idx(vmi + 1);     
        c2(i) = idx_vmi - ee_p1(i, 2);
      }
      else {
        double n_ee_p1 = n - ee_p1(i, 2) + 1.0;
        c2(i) = n_ee_p1;
      }
    }
  }
  
  arma::vec c1_c2 = c1 + c2;
  
  blnc.row(0) = (c1 / c1_c2).t();
  blnc.row(1) = (c2 / c1_c2).t();
  blnc.row(2) = c1_c2;
  
  return blnc;
}

arma::uvec finding_cp(
  const arma::cube &decomp_hist,
  const arma::uvec &sameboat,
  const arma::uword &n
) {

  arma::umat all_edges(1, 1);
  arma::umat edges =
      arma::conv_to<arma::umat>::from(arma::vectorise(decomp_hist.row(0)));

  if (sameboat.n_elem == 1) {
    all_edges = arma::join_vert(edges, sameboat);
  }
  else {
    all_edges = arma::join_vert(arma::reshape(edges, 3, n - 2), sameboat.t());
  }
  
  arma::umat survived_edges(0, 1);
  
  if (all_edges.n_cols > 1) {
    survived_edges = all_edges.cols(arma::find(
        arma::abs(arma::vectorise(decomp_hist.tube(2, 0))) > arma::datum::eps));
  }
  
  static const arma::uvec row_13 = {0, 2};
  static const arma::uvec row_12 = {0, 1};
  static const arma::uvec row_123 = {0, 1, 2};
  static const arma::uvec iv = {0};
  
  arma::umat cp(0, 1);
  
  if ((survived_edges.n_elem > 0) && (survived_edges.n_cols > 1)) {
    arma::uword i = 0;
    
    while (i < survived_edges.n_cols) {
      arma::uvec iv = {i};
      arma::uvec survived_edges_sub1 =
          arma::diff(arma::vectorise(survived_edges(0, i, arma::size(3, 1))));

      if (i + 1 < survived_edges.n_cols) {
        arma::uvec survived_edges_sub2 = arma::diff(
            arma::vectorise(survived_edges(0, i + 1, arma::size(3, 1))));

        if ((survived_edges(3, i) != 0) && (survived_edges_sub1(0) == 1) &&
            (survived_edges_sub2(0) == 1)) {
          cp = arma::join_vert(cp, survived_edges(row_13, iv));
          i += 2;
          continue; // Continue to the next iteration of the while loop
        } else if ((survived_edges(3, i) != 0) &&
                   (survived_edges_sub1(1) == 1) &&
                   (survived_edges_sub2(0) == 1)) {
          cp = arma::join_vert(cp, survived_edges(row_13, iv + 1));
          i += 2;
          continue; // Continue to the next iteration of the while loop
        }
      }

      bool matched = false;
      if ((survived_edges(3, i) == 0) && (survived_edges_sub1(0) == 1) &&
          (survived_edges_sub1(1) == 1)) {
        arma::mat part_obj0 = decomp_hist.tube(0, 0, arma::size(1, 2));
        arma::umat part_obj1 = arma::conv_to<arma::umat>::from(part_obj0);
        arma::uvec target_pair =
            arma::vectorise(survived_edges(1, i, arma::size(2, 1)));

        for (arma::uword j = 0; j < part_obj1.n_cols; ++j) {
          if (arma::all(part_obj1.col(j) == target_pair)) {
            matched = true;
            break;
          }
        }
      }

      if ((survived_edges(3, i) == 0) && (survived_edges_sub1(0) == 1) &&
          (survived_edges_sub1(1) != 1)) {
        cp = arma::join_vert(cp, survived_edges(row_13, iv));
      } else if ((survived_edges(3, i) == 0) && (survived_edges_sub1(0) == 1) &&
                 (survived_edges_sub1(1) == 1) && matched) {
        cp = arma::join_vert(cp, survived_edges(row_12, iv));
      } else {
        cp = arma::join_vert(cp, survived_edges(row_123, iv));
      }

      i = i + 1;
    }
    
    cp = arma::unique(cp);
    cp.resize(cp.n_elem + 1);
    cp(cp.n_elem - 1) = n + 1;
  }
  else if ((survived_edges.n_elem > 0)) {
    arma::uvec survived_edges_sub0 =
        arma::diff(survived_edges(0, 0, arma::size(3, 1)));

    if ((survived_edges.n_cols == 1) && (survived_edges.row(3)(0) == 0) &&
        (survived_edges_sub0(0) == 1) && (survived_edges_sub0(1) != 1)) {
      cp = arma::join_vert(cp, survived_edges(row_13, iv));
    } else if ((survived_edges.n_cols == 1) &&
               (survived_edges.row(3)(0) == 0) &&
               (survived_edges_sub0(0) == 1) && (survived_edges_sub0(1) == 1)) {
      cp = arma::join_vert(cp, survived_edges(row_12, iv));
    }
  }
  else {
    arma::umat cp(0, 1);
  }
  
  if ((n == 3) && (cp.n_elem > 0) && (survived_edges.n_cols == 1)) {
    cp.shed_row(0);
    cp.resize(cp.n_elem + 1);
    cp(cp.n_elem - 1) = n;
  }
  else {
    cp = cp(arma::find((cp <= n) && (cp > 1)));
  }
  
  if (cp.n_elem > 0) {
    cp -= 1;
  }
  return cp;
}
