#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <stdexcept>
#include "hits.h"
#include "hits_misc.h"
#include "hilandyn2_misc.h"

arma::uvec hd_bts_dcmp(
  arma::cube& decomp_hist,
  arma::mat& bts_coeffs,
  const arma::vec& w,
  const double p,
  const arma::uword d,
  const arma::uword n
) {
  
  arma::uword noe = n - 2;
  arma::vec weights_const = arma::ones<arma::vec>(n);
  arma::vec weights_lin = arma::regspace<arma::vec>(1, n);
  arma::uvec idx = arma::regspace<arma::uvec>(0, n-1);              
  arma::uvec paired;

  arma::umat edges(noe, 4, arma::fill::zeros);
  edges.col(0) = arma::regspace<arma::uvec>(0, n - 3);
  edges.col(1) = arma::regspace<arma::uvec>(1, n - 2);
  edges.col(2) = arma::regspace<arma::uvec>(2, n - 1);

  arma::uword steps_left = n - 2;
  arma::uword current_step = 0;
  arma::uvec sameboat;
  
  static const arma::uvec ind0 = {0};
  static const arma::uvec ind1 = {1};
  static const arma::uvec ind2 = {2};
  static const arma::uvec ind3 = {3};
  static const arma::uvec ind01 = {0,1};
  static const arma::uvec ind012 = {0,1,2};
  static const arma::uvec ind345 = {3,4,5};
  static const arma::uvec ind147 = {1,4,7};
  static const arma::uvec ind258 = {2,5,8};
  
  arma::mat Dmat_pre(edges.n_rows, d, arma::fill::none);
  
  while (edges.n_rows > 0) {
    
    double steps_left_d = steps_left;
    arma::uword max_current_steps = static_cast<arma::uword>(std::ceil(p * steps_left_d));
    arma::uvec removable_nodes(idx.max() + 1, arma::fill::ones);
    arma::uword current_rows = edges.n_rows;

    auto Dmat = Dmat_pre.head_rows(current_rows);
    arma::mat sub_wc_wl, sub_tc, detcoef, M0;
    
    if (any(edges.col(3) > 0)) {
      arma::uvec prv = arma::find(edges.col(3) != 0);
      arma::uword ncol_pr = prv.n_elem / 2;
      arma::umat pr = arma::reshape(prv, 2, ncol_pr);
      arma::uvec pr_row1 = pr.row(0).t();
      arma::uvec pr_row2 = pr.row(1).t();

      sub_wc_wl.set_size(ncol_pr, 6);
      sub_tc.set_size(ncol_pr, 3);
      M0.set_size(ncol_pr, 9);
      arma::mat upd_bts_coeffs(ncol_pr, 3, arma::fill::none);
      arma::mat wc_wl_mat(ncol_pr, 6, arma::fill::none);
      arma::mat details_mat(ncol_pr, 2, arma::fill::none);
      
      arma::uword sub_num = ncol_pr * 3;
      arma::uvec sub_ind = arma::regspace<arma::uvec>(0, sub_num - 1);
      
      for (arma::uword j = 0; j < d; ++j) {
        arma::uvec row_j = {j};

        // start computation of detail coefficients
        arma::uvec ind = edges.submat(pr_row1, ind012).as_col();
        sub_wc_wl(sub_ind) = weights_const(ind);
        sub_wc_wl(sub_ind + sub_num) = weights_lin(ind);
        sub_tc(sub_ind) = bts_coeffs(row_j, ind).t();
        
        detcoef = filter_bts_mat(sub_wc_wl);
        arma::vec details = arma::sum(detcoef % sub_tc, 1);
        // end computation of detail coefficients
        
        for (arma::uword k = 0; k < detcoef.n_rows; ++k) {
          M0.row(k) = orth_matrix(detcoef.row(k).t()).as_row();
        }
        
        upd_bts_coeffs.col(0) = sum(sub_tc % M0.cols(ind147), 1);
        upd_bts_coeffs.col(1) = sum(sub_tc % M0.cols(ind258), 1);
        upd_bts_coeffs.col(2) = bts_coeffs(row_j, edges(pr_row2, ind2)).t();

        wc_wl_mat.col(0) = sum(sub_wc_wl.cols(ind012) % M0.cols(ind147), 1);
        wc_wl_mat.col(1) = sum(sub_wc_wl.cols(ind012) % M0.cols(ind258), 1);
        wc_wl_mat.col(2) = weights_const(edges(pr_row2, ind2));
        wc_wl_mat.col(3) = sum(sub_wc_wl.cols(ind345) % M0.cols(ind147), 1);
        wc_wl_mat.col(4) = sum(sub_wc_wl.cols(ind345) % M0.cols(ind258), 1);
        wc_wl_mat.col(5) = weights_lin(edges(pr_row2, ind2));
        
        details_mat.col(0) = arma::abs(details);
        details_mat.col(1) =
            arma::abs(sum(filter_bts_mat(wc_wl_mat) % upd_bts_coeffs, 1));
        arma::vec p_detail = arma::repelem(arma::max(details_mat, 1), 2, 1);

        if (pr.n_elem != edges.n_rows) {
          arma::uvec edgerow_cd3 =
              arma::regspace<arma::uvec>(0, edges.n_rows - 1);
          edgerow_cd3.shed_rows(arma::vectorise(pr));
          arma::uword edgerow_cd3_num = edgerow_cd3.n_elem;

          // start computation of detail coefficients
          arma::mat sub_wc_wl(edgerow_cd3_num, 6, arma::fill::none); 
          arma::mat sub_tc(edgerow_cd3_num, 3, arma::fill::none);
          arma::mat detcoef(edgerow_cd3_num, 3, arma::fill::none);

          arma::uvec ind = edges.submat(edgerow_cd3, ind012).as_col();
          arma::uword sub_num = edgerow_cd3_num * 3;
          arma::uvec sub_ind = arma::regspace<arma::uvec>(0, sub_num - 1);
          sub_wc_wl(sub_ind) = weights_const(ind);
          sub_wc_wl(sub_ind + sub_num) = weights_lin(ind);
          sub_tc(sub_ind) = bts_coeffs(row_j, ind).t();
          
          detcoef = filter_bts_mat(sub_wc_wl);
          arma::vec details = arma::sum(detcoef % sub_tc, 1);
          // end computation of detail coefficients

          Dmat.col(j) = arma::join_vert(p_detail, details);
        }
        else {
          Dmat.col(j) = p_detail;
        }
      }
    }
    else {
      arma::uvec edgerow = arma::regspace<arma::uvec>(0, edges.n_rows - 1);
      arma::mat sub_m(edges.n_rows, 6, arma::fill::none);

      sub_m.col(0) = weights_const(edges(edgerow, ind0));
      sub_m.col(1) = weights_const(edges(edgerow, ind1));
      sub_m.col(2) = weights_const(edges(edgerow, ind2));
      sub_m.col(3) = weights_lin(edges(edgerow, ind0));
      sub_m.col(4) = weights_lin(edges(edgerow, ind1));
      sub_m.col(5) = weights_lin(edges(edgerow, ind2));
      
      arma::mat detcoef = filter_bts_mat(sub_m);
      
      arma::mat bts_coeffs_sub_0 = bts_coeffs.cols(edges.col(0)).t();
      arma::mat bts_coeffs_sub_1 = bts_coeffs.cols(edges.col(1)).t();
      arma::mat bts_coeffs_sub_2 = bts_coeffs.cols(edges.col(2)).t();

      bts_coeffs_sub_0.each_col() %= detcoef.col(0);
      bts_coeffs_sub_1.each_col() %= detcoef.col(1);
      bts_coeffs_sub_2.each_col() %= detcoef.col(2);
      
      // Sum the results
      Dmat = bts_coeffs_sub_0 + bts_coeffs_sub_1 + bts_coeffs_sub_2;
    }
    
    // weighted L2 norm
    arma::mat Dmat_sq = arma::square(Dmat);
    Dmat_sq.each_row() %= w.t();
    arma::uvec ord_det =
        arma::stable_sort_index(arma::sqrt(arma::sum(Dmat_sq, 1)));

    arma::umat cand(ord_det.n_elem, 2, arma::fill::none);
    cand.col(0) = ord_det;
    cand.col(1) = edges(ord_det, ind3);
    
    arma::uvec eitr;
    arma::uword tei = 0;
    
    if (cand(0, 1) > 0) {
      removable_nodes(edges(ord_det.subvec(0,1), ind0)).zeros();
      removable_nodes(edges(ord_det.subvec(0,1), ind1)).zeros();
      removable_nodes(edges(ord_det.subvec(0,1), ind2)).zeros();
      tei += 1;
      eitr = {0, tei}; 
    }
    else {
      removable_nodes(edges(ord_det(0), 0)) =
          removable_nodes(edges(ord_det(0), 1)) =
              removable_nodes(edges(ord_det(0), 2)) = 0;
      eitr = {0};
    }
    
    while ((eitr.n_elem < max_current_steps) & (tei < noe)) {
      tei += 1;
      
      if (cand(tei, 1) > 0) {
        arma::uvec edges_ind1 = edges(ord_det.subvec(tei, tei + 1), ind0);
        arma::uvec edges_ind2 = edges(ord_det.subvec(tei, tei + 1), ind1);
        arma::uvec edges_ind3 = edges(ord_det.subvec(tei, tei + 1), ind2);
        
        if ((arma::accu(removable_nodes(edges_ind1)) == 2) && 
            (arma::accu(removable_nodes(edges_ind2)) == 2) && 
            (arma::accu(removable_nodes(edges_ind3)) == 2)) {
          
          removable_nodes(edges_ind1).zeros();
          removable_nodes(edges_ind2).zeros();
          removable_nodes(edges_ind3).zeros();
          eitr = arma::join_vert(eitr, arma::uvec({tei, tei + 1}));
          tei += 1;
        }
      }
      else {
        arma::uvec edges_ind1(1), edges_ind2(1), edges_ind3(1);
        
        edges_ind1 = edges(ord_det(tei), 0);
        edges_ind2 = edges(ord_det(tei), 1);
        edges_ind3 = edges(ord_det(tei), 2);
        
        if (removable_nodes(edges_ind1(0)) == 1 && 
            removable_nodes(edges_ind2(0)) == 1 && 
            removable_nodes(edges_ind3(0)) == 1) {
            
          removable_nodes(edges_ind1).zeros();
          removable_nodes(edges_ind2).zeros();
          removable_nodes(edges_ind3).zeros();
          eitr = arma::join_vert(eitr, arma::uvec({tei}));
        }
      }
    }
    
    arma::uvec details_min_ind = ord_det(eitr);
    arma::uword no_of_current_steps = eitr.n_elem;
    arma::umat ee =
        arma::reshape(edges.rows(details_min_ind), no_of_current_steps, 4);
    arma::uvec idx0(idx.n_elem);
    idx0 = idx;
    
    arma::uvec ee_col4_ind = arma::find(ee.col(3) > 0);
    
    if (ee_col4_ind.is_empty()) {
      sameboat = arma::join_vert(sameboat, ee.col(3));
      ee.shed_col(3);
      upd_res udt;
      
      arma::uword slice0 = current_step;
      arma::uword slice1 = current_step + no_of_current_steps - 1;

      for (arma::uword j = 0; j < d; ++j) {
        
        udt = updating(ee, weights_const, weights_lin, bts_coeffs.row(j), idx0);
        bts_coeffs.row(j) = udt.bts_coeffs;

        decomp_hist(arma::span(j * 4), arma::span::all,
                    arma::span(slice0, slice1)) =
            arma::conv_to<arma::mat>::from(ee.t() + 1);
        decomp_hist(arma::span(j * 4 + 1), arma::span::all,
                    arma::span(slice0, slice1)) = udt.h;
        decomp_hist(arma::span(j * 4 + 2), arma::span::all,
                    arma::span(slice0, slice1)) = udt.tc1;
        decomp_hist(arma::span(j * 4 + 3), arma::span::all,
                    arma::span(slice0, slice1)) =
            balance_np(paired, ee, idx0, no_of_current_steps, n);
      }
      
      weights_const = udt.weights_const;
      weights_lin = udt.weights_lin;
      idx = udt.idx;
    }
    else {
      sameboat = arma::join_vert(
          sameboat, arma::vectorise(ee(arma::find(ee.col(3) != 0), ind3)));

      arma::uvec pr0 = arma::find(ee.col(3) != 0);
      arma::umat pr = arma::reshape(pr0, 2, pr0.n_elem / 2);

      ee(pr.row(1), ind01) = ee(pr.row(0), ind01);
      
      arma::umat ee_p1 = ee.rows(pr.row(0));
      arma::umat ee_p2 = ee.rows(pr.row(1));
      ee_p1.shed_col(3); ee_p2.shed_col(3);
      upd_res udt;
      
      arma::uword slice0 = current_step;
      arma::uword slice1 = current_step + ee_p1.n_rows - 1;
      arma::uword slice2 = current_step + ee_p1.n_rows;
      arma::uword slice3 = current_step + pr.n_elem - 1;
      arma::uword slice4 = current_step + pr.n_elem;
      arma::uword slice5 = current_step + no_of_current_steps - 1;
      
      for (arma::uword j = 0; j < d; ++j) {

        udt = updating(ee_p1, weights_const, weights_lin, bts_coeffs.row(j),
                       idx0);

        decomp_hist(arma::span(j * 4), arma::span::all,
                    arma::span(slice0, slice1)) =
            arma::conv_to<arma::mat>::from(ee_p1.t() + 1);
        decomp_hist(arma::span(j * 4 + 1), arma::span::all,
                    arma::span(slice0, slice1)) = udt.h;
        decomp_hist(arma::span(j * 4 + 2), arma::span::all,
                    arma::span(slice0, slice1)) = udt.tc1;
        decomp_hist(arma::span(j * 4 + 3), arma::span::all,
                    arma::span(slice0, slice1)) =
            balance_p(pr, ee_p1, ee_p2, idx0, n);

        udt = updating(ee_p2, udt.weights_const, udt.weights_lin,
                       udt.bts_coeffs, udt.idx);
        idx = udt.idx;

        decomp_hist(arma::span(j * 4), arma::span::all,
                    arma::span(slice2, slice3)) =
            arma::conv_to<arma::mat>::from(ee_p2.t() + 1);
        decomp_hist(arma::span(j * 4 + 1), arma::span::all,
                    arma::span(slice2, slice3)) = udt.h;
        decomp_hist(arma::span(j * 4 + 2), arma::span::all,
                    arma::span(slice2, slice3)) = udt.tc1;
        decomp_hist(arma::span(j * 4 + 3), arma::span::all,
                    arma::span(slice2, slice3)) =
            balance_p(pr, ee_p1, ee_p2, idx, n);

        if (pr.n_elem != ee.n_rows) {
          sameboat = arma::join_vert(
              sameboat, arma::vectorise(ee(arma::find(ee.col(3) == 0), ind3)));
          arma::umat ee_np = ee; ee_np.shed_rows(arma::vectorise(pr));

          udt = updating(ee_np, udt.weights_const, udt.weights_lin,
                         udt.bts_coeffs, idx);
          bts_coeffs.row(j) = udt.bts_coeffs;
          
          arma::uword no_of_current_steps_pr = no_of_current_steps - pr.n_elem;

          decomp_hist(arma::span(j * 4), arma::span::all,
                      arma::span(slice4, slice5)) =
              arma::conv_to<arma::mat>::from(ee_np.t() + 1);
          decomp_hist(arma::span(j * 4 + 1), arma::span::all,
                      arma::span(slice4, slice5)) = udt.h;
          decomp_hist(arma::span(j * 4 + 2), arma::span::all,
                      arma::span(slice4, slice5)) = udt.tc1;
          decomp_hist(arma::span(j * 4 + 3), arma::span::all,
                      arma::span(slice4, slice5)) =
              balance_np(paired, ee_np, idx, no_of_current_steps_pr, n);
        }
        else {
          bts_coeffs.row(j) = udt.bts_coeffs;
        }
      }
      weights_const = udt.weights_const;
      weights_lin = udt.weights_lin;
      idx = udt.idx;
      ee.shed_col(3);
    }
    
    paired =
        arma::unique(arma::join_vert(paired, arma::vectorise(ee.cols(0, 1))));
    arma::vec paired_ee = match(paired, vectorise(ee.col(2)));
    
    if (!arma::find_finite(paired_ee).is_empty()) {
      paired = arma::sort(paired(arma::find_nonfinite(paired_ee)));
    }
    
    arma::umat edges2(idx.n_elem - 2, 3, arma::fill::none);
    
    if (edges2.n_rows == 0) {
      break;
    }
    
    edges = edges2;
    edges.col(0) = idx.subvec(0, idx.n_elem - 3);
    edges.col(1) = idx.subvec(1, idx.n_elem - 2);
    edges.col(2) = idx.subvec(2, idx.n_elem - 1);
    
    arma::vec edges_paired = match(arma::vectorise(edges), paired);
    edges_paired(arma::find_finite(edges_paired)).ones();
    edges_paired(arma::find_nonfinite(edges_paired)).zeros();

    arma::mat matchpair =
        arma::reshape(edges_paired, edges_paired.n_elem / 3, 3).t();
    arma::uvec rs = arma::find(arma::sum(matchpair, 0) == 3);

    if (rs.n_elem > 0) {
      arma::umat edges_rs = edges; 
      edges_rs.shed_rows(rs);
      
      arma::mat matchpair_rs = matchpair; 
      matchpair_rs.shed_cols(rs);

      edges = arma::join_vert(edges.rows(rs), edges_rs);
      matchpair = arma::join_horiz(matchpair.cols(rs), matchpair_rs);

      arma::uvec edges_seq = arma::repelem(
          arma::regspace<arma::uvec>(1, rs.n_elem / 2), 2, 1);
      arma::uvec edges_zero(edges.n_rows - rs.n_elem);
      edges = arma::join_horiz(edges, arma::join_vert(edges_seq, edges_zero));
    }
    else {
      arma::uvec edges_zero(edges.n_rows);
      edges = arma::join_horiz(edges, edges_zero);
    }

    arma::uvec removed = arma::join_vert(
        arma::find(arma::sum(matchpair, 0) == 1),
        arma::find((matchpair.row(0) == 1) && (matchpair.row(1) == 0) &&
                   (matchpair.row(2) == 1)));

    if (removed.n_elem > 0) {
      edges.shed_rows(removed);
    }
    
    noe = edges.n_rows;
    steps_left = steps_left - no_of_current_steps;
    current_step = current_step + no_of_current_steps;

    if (noe == 1) {
      edges.col(3) = 0;
    }
  }
  return sameboat;
}

detail_res compute_det(arma::cube& decomp_hist) {
  const arma::uword d = decomp_hist.n_rows / 4;
  const arma::uword s = decomp_hist.n_slices;
  arma::uvec slice_ind = arma::regspace<arma::uvec>(1, d) * 4 - 2;
  arma::mat detail_all(d, s);
  arma::vec details(s);
  
  if (d > 1) {
    for (arma::uword i = 0; i < s; i++) {
      // L2 norm
      detail_all.col(i) = decomp_hist.slice(i)(slice_ind);
      details[i] = arma::norm(detail_all.col(i));
    }
  }
  else {
    detail_all = arma::vectorise(decomp_hist.tube(2, 0)).t();
    details = arma::abs(detail_all.as_col());
  }
  
  detail_res result;
  result.details = details;
  result.detail_all = detail_all;
  return result;
}

arma::cube hd_bts_dns(
  arma::cube decomp_hist,
  const arma::uvec& sameboat,
  const detail_res& details_out,
  const arma::uword& n,
  const double& lambda,
  const double& bal,
  const double& seg_len_min
) {

  const arma::uword d = decomp_hist.n_rows / 4;
  arma::uvec slice_ind = arma::regspace<arma::uvec>(1, d) * 4 - 2;

  // 1) connected rule
  arma::uvec protect(n, arma::fill::zeros);

  for (arma::uword i = 0; i < (n - 2); i++) {
    
    arma::uword ind_c1 = static_cast<arma::uword>(decomp_hist(0, 0, i)) - 1;
    arma::uword ind_c2 = static_cast<arma::uword>(decomp_hist(0, 1, i)) - 1;
    arma::uword ind_c3 = static_cast<arma::uword>(decomp_hist(0, 2, i)) - 1;

    if ((protect(ind_c1) == 0) && (protect(ind_c2) == 0) &&
        (protect(ind_c3) == 0)) {

      bool seg_len_min_cond;
      double h43 = decomp_hist(3, 2, i);
      if (h43 < 1) {
        seg_len_min_cond = (1.0 >= seg_len_min);
      } else {
        double val1 = decomp_hist(3, 0, i) * h43;
        double val2 = decomp_hist(3, 1, i) * h43;
        seg_len_min_cond = (std::min(val1, val2) >= seg_len_min);
      }

      arma::uword logi_val = (details_out.details(i) > lambda) &&
                             (decomp_hist(3, 0, i) > bal) &&
                             (decomp_hist(3, 1, i) > bal && seg_len_min_cond);
      arma::uvec logi_vec(slice_ind.n_elem, arma::fill::value(logi_val));
      decomp_hist.slice(i)(slice_ind) =
          decomp_hist.slice(i)(slice_ind) % logi_vec;
    }

    if (std::abs(decomp_hist(2, 0, i)) > 0) {     
      protect(ind_c1) = 1;
      protect(ind_c2) = 1;
    }
  }
  
  // 2) two-together rule
  arma::uvec paired0 = arma::find(sameboat != 0);
  
  if (paired0.n_elem > 0)  {
    arma::umat paired = arma::reshape(paired0, 2, paired0.n_elem / 2);
    
    for (arma::uword i = 0; i < paired.n_cols; i++) {

      arma::vec overzero = match(
          paired.col(i), arma::find(arma::abs(decomp_hist.tube(2, 0)) > 0));
      arma::vec zero = match(
          paired.col(i), arma::find(arma::abs(decomp_hist.tube(2, 0)) == 0));
      arma::uvec overzero1 = arma::find_finite(overzero);
      arma::uvec zero1 = arma::find_finite(zero);
      arma::uvec sl = paired(arma::find_nonfinite(overzero), arma::uvec{i});
      
      if ((overzero1.n_elem == 1) & (zero1.n_elem == 1)) {
        decomp_hist.slice(arma::as_scalar(sl))(slice_ind) =
            details_out.detail_all.col(arma::as_scalar(sl));
      }
    }
  }
  
  return decomp_hist;
}

arma::mat hd_bts_inv(
  arma::cube& decomp_hist,
  arma::mat bts_coeffs,
  const arma::uword& n
) {
  
  const arma::uword d = decomp_hist.n_rows / 4;
  const arma::uvec slice_ind_1 = arma::regspace<arma::uvec>(1, d) * 4 - 2;
  const arma::uvec add_row(slice_ind_1.n_elem,
                           arma::fill::value(decomp_hist.n_rows));
  const arma::uvec slice_ind_2 = slice_ind_1 + add_row;
  const arma::uvec slice_ind_3 = slice_ind_2 + add_row;
  const arma::uvec slice_ind_123 =
      arma::join_vert(slice_ind_1, slice_ind_2, slice_ind_3);

  arma::mat33 inv_mat;

  for (arma::uword i = (n - 2); i-- > 0; ) {

    inv_mat =
        orth_matrix(arma::vectorise(decomp_hist(1, 0, i, arma::size(1, 3, 1))));

    arma::uvec ind = arma::conv_to<arma::uvec>::from(arma::vectorise(
                         decomp_hist(0, 0, i, arma::size(1, 3, 1)))) -
                     1;

    decomp_hist.slice(i)(slice_ind_2) = bts_coeffs.col(ind(0));
    decomp_hist.slice(i)(slice_ind_3) = bts_coeffs.col(ind(1));

    arma::mat tmp = arma::reshape(decomp_hist.slice(i)(slice_ind_123),
                                  slice_ind_1.n_elem, 3);
    arma::mat rcstr_tmp(1, 1);
    
    if (d == 1) {
      rcstr_tmp = inv_mat * tmp.as_col();
    }
    else {
      rcstr_tmp = inv_mat * tmp.t();
    }

    bts_coeffs.col(ind(0)) = rcstr_tmp.row(0).t();
    bts_coeffs.col(ind(1)) = rcstr_tmp.row(1).t();
    bts_coeffs.col(ind(2)) = rcstr_tmp.row(2).t();
  }
  return bts_coeffs;
}

// Effective dimension via participation ratio of row covariance
static double participation_ratio_deff(const arma::mat& X) {
  const arma::uword d = X.n_rows, n = X.n_cols;
  arma::mat Xc = X;
  Xc.each_col() -= arma::mean(Xc, 1);
  
  arma::mat C = (Xc * Xc.t()) / std::max(1.0, static_cast<double>(n - 1));
  C.diag() += 1e-8;
  const double s1 = arma::trace(C);
  const double s2 = arma::accu(arma::square(C));
  double deff = (s2 > 0.0) ? (s1 * s1 / s2) : 1.0;
  
  return std::max(1.0, std::min(deff, static_cast<double>(d)));
}

static double BIC_segmented_piecewise_linear(
  double rss,
  arma::uword n,
  arma::uword K,
  double d_eff,
  bool count_break_locations = true
) {
  double nn = std::max(2.0, static_cast<double>(n));
  double de = std::max(1.0, d_eff);
  
  // Discontinuous piecewise linear: 2*(K+1)
  double p_per_dim = static_cast<double>(2 * (K + 1));
  
  // Total regression parameters
  double p_model = de * p_per_dim;
  
  // Include K break-location parameters in count
  if (count_break_locations) p_model += static_cast<double>(K);
  
  // BIC
  double fit = nn * std::log(std::max(rss, 1e-12) / nn);
  double pen = p_model * std::log(nn);
  
  return fit + pen;
}

static bool better_candidate(const CandidateCache& a, const CandidateCache& b) {
  return a.bic < b.bic;
}

static arma::vec build_lambda_path(
  const arma::vec& details,
  double lambda_min,
  double lambda_max
) {
  arma::uvec finite_ind = arma::find_finite(details);
  arma::vec events = details.elem(finite_ind);
  
  if (!events.is_empty()) {
    arma::uvec range_ind = arma::find(
      (events > lambda_min) % (events <= lambda_max)
    );
    events = events.elem(range_ind);
  }
  
  events = arma::sort(arma::unique(events), "descend");
  events.transform([lambda_min](double event) {
    return std::max(
      lambda_min,
      std::nextafter(event, -std::numeric_limits<double>::infinity())
    );
  });
  
  arma::vec lambda_path(events.n_elem + 1);
  lambda_path[0] = lambda_max;
  if (!events.is_empty()) {
    lambda_path.tail(events.n_elem) = events;
  }
  
  return lambda_path;
}

static CandidateCache select_with_early_stop(
  const arma::cube &decomp_hist,
  const arma::uvec &sameboat,
  const detail_res &details_out,
  const arma::mat &X,
  const arma::mat &Xn,
  const arma::vec &sd,
  arma::uword n,
  arma::uword d,
  double d_eff,
  double bal,
  double &seg_len_min,
  double th_min,
  double th_max,
  arma::uword cpt_max,
  double dBIC_min,
  bool verbose
) {
  th_min  = std::max(0.0, th_min);
  
  if (th_max == -1) {
    arma::uvec finite_details = arma::find_finite(details_out.details);
    double detail_max = finite_details.is_empty()
      ? 0.0
    : arma::max(details_out.details(finite_details));
    
    double sqrt_d = std::sqrt(static_cast<double>(d));
    th_max = std::max(th_min, detail_max / sqrt_d);
    while (th_max * sqrt_d < detail_max) {
      th_max = std::nextafter(th_max, std::numeric_limits<double>::infinity());
    }
    
    if (verbose) {
      std::cout << "automatic th_max = " << th_max << "\n";
    }
  }
  else {
    th_max = std::max(th_min, th_max);
  }

  arma::vec mad_scale = mad_rows(arma::diff(X, 1, 1));

  // Score a reconstructed model using the original-scale standardized residuals.
  auto evaluate_candidate = [&](
    double threshold,
    arma::mat fitted,
    arma::uvec changepoints
  ) {
    arma::mat fitted_orig = fitted;
    fitted_orig.each_col() %= sd;
    arma::mat resid = X - fitted_orig;
    resid.each_col() /= mad_scale;
    double rss = arma::accu(arma::square(resid));
    double bic = BIC_segmented_piecewise_linear(
      rss,
      n,
      changepoints.n_elem,
      d_eff
    );

    return CandidateCache{threshold, std::move(fitted),
                          std::move(changepoints), rss, bic};
  };

  if (std::abs(th_max - th_min) <= std::numeric_limits<double>::epsilon()) {
    double lambda = th_max * std::sqrt(static_cast<double>(d));
    arma::cube decomp_it = hd_bts_dns(decomp_hist, sameboat, details_out, n,
                                      lambda, bal, seg_len_min);
    arma::mat Y_it = hd_bts_inv(decomp_it, Xn, n);
    arma::uvec cpt_it = finding_cp(decomp_it, sameboat, n);

    return evaluate_candidate(th_max, std::move(Y_it), std::move(cpt_it));
  }
  
  // best candidate per K found so far
  std::unordered_map<arma::uword, CandidateCache> best_by_K;
  arma::uword K_min_seen = std::numeric_limits<arma::uword>::max();
  
  // current selected model (starts at smallest K observed)
  bool has_cur = false;
  CandidateCache cur{0.0, arma::mat(), arma::uvec(), 0.0, 0.0};
  
  const double sqrt_d = std::sqrt(static_cast<double>(d));
  const double lambda_min = th_min * sqrt_d;
  const double lambda_max = th_max * sqrt_d;
  const arma::vec lambda_path = build_lambda_path(
    details_out.details, lambda_min, lambda_max
  );
  
  if (verbose) {
    std::cout << "exact threshold states = " << lambda_path.n_elem << "\n";
  }
  
  for (arma::uword path_ind = 0; path_ind < lambda_path.n_elem; ++path_ind) {
    const double lambda = lambda_path[path_ind];
    const double th = lambda / sqrt_d;

    arma::cube decomp_it = hd_bts_dns(decomp_hist, sameboat, details_out, n,
                                      lambda, bal, seg_len_min);
    arma::mat Y_it = hd_bts_inv(decomp_it, Xn, n);
    arma::uvec cpt_it = finding_cp(decomp_it, sameboat, n);
    arma::uword K = cpt_it.n_elem;
    
    if (K > cpt_max) {
      if (verbose) {
        std::cout << "K = " << K
                  << "; K > cpt_max. Stop changepoint candidate search."
                  << "\n";
        }
      break;
    }
    
    CandidateCache cand = evaluate_candidate(
      th,
      std::move(Y_it),
      std::move(cpt_it)
    );
    
    // update best_by_K
    auto it = best_by_K.find(K);
    if (it == best_by_K.end() || better_candidate(cand, it->second)) {
      best_by_K[K] = cand;
    }
    
    // track smallest K observed
    if (K < K_min_seen) {
      K_min_seen = K;
      cur = best_by_K.at(K);
      has_cur = true;
    }
    
    if (!has_cur) continue;
    
    // Find smallest available K > K_cur
    arma::uword K_cur = cur.cpt.n_elem;
    arma::uword K_target = std::numeric_limits<arma::uword>::max();
    for (const auto& kv : best_by_K) {
      arma::uword k2 = kv.first;
      if (k2 > K_cur && k2 <= cpt_max) K_target = std::min(K_target, k2);
    }

    if (K_target == std::numeric_limits<arma::uword>::max())
      continue;

    const CandidateCache& best_try = best_by_K[K_target];
    double dBIC = cur.bic - best_try.bic;
    double rel_gain = (cur.rss - best_try.rss) / std::max(cur.rss, 1e-12);
    
    if (verbose) {
      std::cout << "th = " << th
                << "  K current=" << K_cur
                << "  K try = " << K_target
                << "  dBIC = " << dBIC
                << "  relative rss gain = " << rel_gain << "\n";
    }
    
    // Move forward if pass
    if (dBIC >= dBIC_min) {
      cur = best_try;
      continue;
    }

    // Early stop
    break;
  }
  
  if (!has_cur) {
    throw std::runtime_error(
      "No candidate satisfies cpt_max; increase th_max or cpt_max."
    );
  }
  
  return cur;
}

// ============================================================
// Main function
// ============================================================
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
) {
  const double p = 0.01;
  const double bal = 0.0;
  const arma::uword d = X.n_rows;
  const arma::uword n = X.n_cols;
  
  // Normalize
  arma::mat Xn = X;
  Xn.each_col() /= sd;
  
  // Row-wise weights validation
  if (weights.n_elem != d || weights.has_nan()) {
    weights = arma::ones(d);
  }
  
  // Decomposition history
  arma::cube decomp_hist(4 * d, 3, n - 2, arma::fill::zeros);
  arma::uvec sameboat = hd_bts_dcmp(decomp_hist, Xn, weights, p, d, n);
  detail_res details_out = compute_det(decomp_hist);
  
  // Effective dimension
  double d_eff = participation_ratio_deff(X);
  
  CandidateCache best = select_with_early_stop(
    decomp_hist, sameboat, details_out,
    X, Xn, sd, n, d, d_eff, bal, seg_len_min,
    th_min, th_max,
    cpt_max, dBIC_min, verbose
  );
  
  cpt = best.cpt;
  th_selected = best.th;
  arma::mat Y_cur = best.Y;
  Y_cur.each_col() %= sd;
  return Y_cur;
}
