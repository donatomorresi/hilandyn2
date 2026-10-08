#include "hilandyn2_deseason.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>


double row_abs_acf_lag(const arma::rowvec& x, arma::uword lag) {
  const arma::uword n = x.n_elem;
  if (lag == 0 || lag >= n) return arma::datum::nan;
  
  arma::vec a = x.cols(0, n - lag - 1).t();
  arma::vec b = x.cols(lag, n - 1).t();
  
  arma::uvec ok_a = arma::find_finite(a);
  arma::uvec ok_b = arma::find_finite(b);
  arma::uvec ok = arma::intersect(ok_a, ok_b);
  if (ok.n_elem < 3) return arma::datum::nan;
  
  a = a.elem(ok);
  b = b.elem(ok);
  a -= arma::mean(a);
  b -= arma::mean(b);
  
  const double aa = arma::dot(a, a);
  const double bb = arma::dot(b, b);
  const double denom = std::sqrt(aa * bb);
  if (denom <= std::numeric_limits<double>::epsilon()) return arma::datum::nan;
  
  return std::abs(arma::dot(a, b) / denom);
}

double seasonal_structure_score(const arma::mat& X, arma::uword period) {
  if (period == 0 || X.n_cols <= period) return arma::datum::inf;

  const arma::mat differences = arma::diff(X, 1, 1);
  arma::vec scores(X.n_rows, arma::fill::value(arma::datum::nan));
  for (arma::uword r = 0; r < X.n_rows; ++r) {
    scores(r) = row_abs_acf_lag(differences.row(r), period);
  }
  
  arma::uvec ok = arma::find_finite(scores);
  if (ok.n_elem == 0) return arma::datum::inf;
  return arma::median(scores.elem(ok));
}

arma::field<arma::vec> make_filters_from_scaling(arma::vec g) {
  if (g.is_empty()) {
    throw std::invalid_argument(
        "make_filters_from_scaling: empty coefficients");
  }

  g /= std::sqrt(2.0);
  arma::vec h = arma::reverse(g);
  if (h.n_elem > 1) {
    h.elem(arma::regspace<arma::uvec>(1, 2, h.n_elem - 1)) *= -1.0;
  }

  arma::field<arma::vec> filters(2);
  filters(0) = h;
  filters(1) = g;
  return filters;
}

std::string trim_copy(const std::string& x);
std::string upper_copy(std::string x);

arma::field<arma::vec> wavelet_filters_modwt(const std::string& filter_name) {
  const std::string key = upper_copy(trim_copy(filter_name));

  /* Coefficients from the PyWavelets package
   * Repository: https://github.com/PyWavelets/pywt (main branch)
   * File: pywt/_extensions/c/wavelets_coeffs.template.h
   */
  if (key == "DB2") {
    return make_filters_from_scaling(arma::vec{
      4.829629131445341433748715998644486838169524195042022752011715e-01,
      8.365163037378079055752937809168732034593703883484392934953414e-01,
      2.241438680420133810259727622404003554678835181842717613871683e-01,
      -1.294095225512603811744494188120241641745344506599652569070016e-01
    });
  }
  if (key == "LA4" || key == "SYM2") {
    return make_filters_from_scaling(arma::vec{
      -1.294095225512603811744494188120241641745344506599652569070016e-01,
      2.241438680420133810259727622404003554678835181842717613871683e-01,
      8.365163037378079055752937809168732034593703883484392934953414e-01,
      4.829629131445341433748715998644486838169524195042022752011715e-01,
    });
  }
  if (key == "DB3") {
    return make_filters_from_scaling(arma::vec{
      3.326705529500826159985115891390056300129233992450683597084705e-01,
      8.068915093110925764944936040887134905192973949948236181650920e-01,
      4.598775021184915700951519421476167208081101774314923066433867e-01,
      -1.350110200102545886963899066993744805622198452237811919756862e-01,
      -8.544127388202666169281916918177331153619763898808662976351748e-02,
      3.522629188570953660274066471551002932775838791743161039893406e-02
    });
  }
  if (key == "LA6" || key == "SYM3") {
    return make_filters_from_scaling(arma::vec{
      3.522629188570953660274066471551002932775838791743161039893406e-02,
      -8.544127388202666169281916918177331153619763898808662976351748e-02,
      -1.350110200102545886963899066993744805622198452237811919756862e-01,
      4.598775021184915700951519421476167208081101774314923066433867e-01,
      8.068915093110925764944936040887134905192973949948236181650920e-01,
      3.326705529500826159985115891390056300129233992450683597084705e-01,
    });
  }
  /* References:
   * D. B. Percival and A. T. Walden, {Wavelet Methods for
   * Time Series Analysis}, Cambridge University Press, 2000.
   */
  if (key == "LA8" || key == "SYM4") {
    return make_filters_from_scaling(arma::vec{
      -0.0757657147893407,
      -0.0296355276459541,
      0.4976186676324578,
      0.8037387518052163,
      0.2978577956055422,
      -0.0992195435769354,
      -0.0126039672622612,
      0.0322231006040713
    });
  }

  throw std::invalid_argument("Unsupported filter name: " + key);
}

std::string trim_copy(const std::string& x) {
  std::size_t first = 0;
  while (first < x.size() &&
         std::isspace(static_cast<unsigned char>(x[first]))) {
    ++first;
  }

  std::size_t last = x.size();
  while (last > first &&
         std::isspace(static_cast<unsigned char>(x[last - 1]))) {
    --last;
  }

  return x.substr(first, last - first);
}

std::string upper_copy(std::string x) {
  for (char& c : x) {
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  }
  return x;
}

void append_filter_candidate(std::vector<std::string>& filters,
                             const std::string& filter_name) {
  std::string normalized = filter_name;
  std::replace(normalized.begin(), normalized.end(), '|', ',');
  std::replace(normalized.begin(), normalized.end(), ';', ',');

  std::stringstream ss(normalized);
  std::string item;
  while (std::getline(ss, item, ',')) {
    item = upper_copy(trim_copy(item));
    if (item.empty()) continue;
    if (std::find(filters.begin(), filters.end(), item) == filters.end()) {
      filters.push_back(item);
    }
  }
}

std::vector<std::string>
filter_candidates_from_names(const std::vector<std::string> &filter_names) {
  std::vector<std::string> filters;
  for (const std::string& filter_name : filter_names) {
    append_filter_candidate(filters, filter_name);
  }

  if (filters.empty()) {
    throw std::invalid_argument("No wavelet filter candidates were provided");
  }
  return filters;
}

std::vector<std::string>
filter_candidates_from_name(const std::string &filter_name) {
  return filter_candidates_from_names(std::vector<std::string>{filter_name});
}

struct WaveletFilterCandidate {
  std::string name;
  arma::field<arma::vec> filters;
  int length;
  int J_eff;
  int filter_width;
};

arma::mat shift_circular(const arma::mat& M, int s) {
  if (M.n_cols == 0) return M;

  const int n = static_cast<int>(M.n_cols);
  const arma::uword offset = static_cast<arma::uword>((s % n + n) % n);
  if (offset == 0) return M;

  return arma::circshift(M, -static_cast<arma::sword>(offset), 1);
}

inline bool modwpt_node_uses_scaling_filter(int node) {
  // Percival-Walden sequency order: low pass for n mod 4 in {0, 3}.
  const int remainder = node % 4;
  return remainder == 0 || remainder == 3;
}

int compute_J_max(std::size_t N, std::size_t L) {
  if (N <= 1 || L <= 1) throw std::invalid_argument("N and L must be > 1");
  double Jd = std::log(((static_cast<double>(N) - 1.0) /
                        (static_cast<double>(L) - 1.0)) + 1.0) / std::log(2.0);
  return static_cast<int>(std::floor(Jd));
}

void modwt_forward(
  const arma::mat& Vin,
  int j,
  const arma::vec& h,
  const arma::vec& g,
  arma::mat& Wj,
  arma::mat& Vj
) {
  const int step = 1 << (j - 1);

  Wj.zeros(Vin.n_rows, Vin.n_cols);
  Vj.zeros(Vin.n_rows, Vin.n_cols);
  for (arma::uword m = 0; m < h.n_elem; ++m) {
    const arma::mat shifted = shift_circular(Vin, -static_cast<int>(m) * step);
    Wj += h(m) * shifted;
    Vj += g(m) * shifted;
  }
}

arma::mat modwt_backward(
  const arma::mat& Wj,
  const arma::mat& Vj,
  int j,
  const arma::vec& h,
  const arma::vec& g
) {
  const int step = 1 << (j - 1);
  arma::mat Vprev(Wj.n_rows, Wj.n_cols, arma::fill::zeros);

  for (arma::uword m = 0; m < h.n_elem; ++m) {
    const int shift = static_cast<int>(m) * step;
    Vprev += h(m) * shift_circular(Wj, shift) +
      g(m) * shift_circular(Vj, shift);
  }
  return Vprev;
}

MODWPTree modwpt_decompose(
  const arma::mat& X,
  int J,
  const arma::vec& h,
  const arma::vec& g
) {
  if (J <= 0) throw std::invalid_argument("J must be positive");

  MODWPTree tree;
  tree.levels.resize(J + 1);
  tree.levels[0].resize(1);
  tree.levels[0][0] = X;

  for (int lvl = 1; lvl <= J; ++lvl) {
    const int num_nodes = 1 << lvl;
    tree.levels[lvl].resize(num_nodes);
    for (int parent = 0; parent < (num_nodes >> 1); ++parent) {
      const arma::mat& VinMat = tree.levels[lvl - 1][parent];
      arma::mat low;
      arma::mat high;
      modwt_forward(VinMat, lvl, h, g, high, low);
      const int first_child = 2 * parent;
      const int second_child = first_child + 1;
      if (modwpt_node_uses_scaling_filter(first_child)) {
        tree.levels[lvl][first_child] = low;
        tree.levels[lvl][second_child] = high;
      } else {
        tree.levels[lvl][first_child] = high;
        tree.levels[lvl][second_child] = low;
      }
    }
  }
  return tree;
}

bool node_is_descendant(
  int child_level,
  int child_node,
  int parent_level,
  int parent_node
) {
  if (child_level <= parent_level) return false;
  return (child_node >> (child_level - parent_level)) == parent_node;
}

arma::imat canonicalize_packet_nodes(const arma::imat& nodes) {
  if (nodes.n_cols < 2) return arma::imat(0, 2, arma::fill::zeros);

  std::vector<std::pair<int, int>> sorted;
  sorted.reserve(nodes.n_rows);
  for (arma::uword i = 0; i < nodes.n_rows; ++i) {
    const int lvl = nodes(i, 0);
    const int k = nodes(i, 1);
    if (lvl < 0 || k < 0 || k >= (1 << lvl)) {
      throw std::invalid_argument(
          "seasonal_nodes contains an invalid MODWPT node");
    }
    sorted.push_back({lvl, k});
  }

  std::sort(sorted.begin(), sorted.end());
  sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());

  std::vector<std::pair<int, int>> out;
  for (const auto& node : sorted) {
    bool covered_by_parent = false;
    for (const auto& kept : out) {
      if (node_is_descendant(node.first, node.second, kept.first,
                             kept.second)) {
        covered_by_parent = true;
        break;
      }
    }
    if (!covered_by_parent) out.push_back(node);
  }

  arma::imat ans(static_cast<arma::uword>(out.size()), 2, arma::fill::zeros);
  for (arma::uword i = 0; i < ans.n_rows; ++i) {
    ans(i, 0) = out[i].first;
    ans(i, 1) = out[i].second;
  }
  return ans;
}

int max_packet_node_level(const arma::imat& nodes) {
  if (nodes.n_cols < 2 || nodes.n_rows == 0) return 0;
  return nodes.col(0).max();
}

arma::imat select_packet_node_candidates_for_period(
  double period,
  double dt,
  int max_level
) {
  if (max_level < 1) throw std::invalid_argument("max_level must be >= 1");
  if (dt <= 0.0) throw std::invalid_argument("dt must be positive");
  if (period <= 0.0) throw std::invalid_argument("period must be positive");

  const double nyq = 1.0 / (2.0 * dt);
  const double eps = std::max(1e-14, nyq * 1e-12);
  const double target_f = 1.0 / period;
  if (target_f > nyq + eps) {
    throw std::invalid_argument(
        "period is shorter than the Nyquist-resolvable period 2 * dt");
  }

  const int bins = 1 << max_level;
  const arma::vec centers =
    (arma::regspace<arma::vec>(0, bins - 1) + 0.5) *
    (nyq / static_cast<double>(bins));
  const arma::vec distances = arma::abs(centers - target_f);
  const double best_dist = distances.min();
  const double tie_tol = std::max(1e-14, target_f * 1e-12);
  const arma::uvec closest = arma::find(
    arma::abs(distances - best_dist) <= tie_tol);

  arma::imat raw(closest.n_elem, 2, arma::fill::zeros);
  raw.col(0).fill(max_level);
  raw.col(1) = arma::conv_to<arma::ivec>::from(closest);
  return canonicalize_packet_nodes(raw);
}

arma::field<arma::imat> make_objective_node_candidates(
  double period,
  double dt,
  int J_eff
) {
  arma::field<arma::imat> node_candidates;

  // Conservative search: nearest packet center plus any exact-distance tie.
  arma::imat candidate_nodes =
      select_packet_node_candidates_for_period(period, dt, J_eff);
  node_candidates.set_size(candidate_nodes.n_rows);
  for (arma::uword i = 0; i < candidate_nodes.n_rows; ++i) {
    node_candidates(i) = candidate_nodes.rows(i, i);
  }

  return node_candidates;
}

arma::mat modwpt_reconstruct_single_node(
  const arma::mat& coef,
  int level,
  int node,
  const arma::vec& h,
  const arma::vec& g
) {
  arma::mat current = coef;
  int cur_node = node;

  for (int lvl = level; lvl >= 1; --lvl) {
    arma::mat low(current.n_rows, current.n_cols, arma::fill::zeros);
    arma::mat high(current.n_rows, current.n_cols, arma::fill::zeros);
    if (modwpt_node_uses_scaling_filter(cur_node)) {
      low = current;
    } else {
      high = current;
    }

    current = modwt_backward(high, low, lvl, h, g);
    cur_node /= 2;
  }
  return current;
}

arma::mat modwpt_reconstruct_mixed_nodes(
  const MODWPTree& wt,
  const arma::imat& keep_nodes,
  const arma::vec& h,
  const arma::vec& g
) {
  const int J = static_cast<int>(wt.levels.size()) - 1;
  arma::mat out(
    wt.levels[0][0].n_rows,
    wt.levels[0][0].n_cols,
    arma::fill::zeros
  );
  if (keep_nodes.n_cols < 2 || keep_nodes.n_rows == 0) return out;

  const arma::imat nodes = canonicalize_packet_nodes(keep_nodes);

  for (arma::uword i = 0; i < nodes.n_rows; ++i) {
    const int lvl = nodes(i, 0);
    const int k = nodes(i, 1);
    if (lvl < 0 || lvl > J || k < 0 || k >= (1 << lvl)) {
      throw std::invalid_argument("keep_nodes out of range");
    }

    arma::mat coef = wt.levels[lvl][k];
    out += modwpt_reconstruct_single_node(coef, lvl, k, h, g);
  }
  return out;
}

arma::mat pad_reflect(const arma::mat& X, arma::uword pad) {
  const arma::uword n = X.n_cols;
  if (pad == 0) return X;
  if (n < 2 || pad >= n) {
    throw std::invalid_argument(
        "pad_reflect: pad must be smaller than the number of columns");
  }

  return arma::join_rows(
    arma::fliplr(X.cols(1, pad)),
    X,
    arma::fliplr(X.cols(n - pad - 1, n - 2)));
}

int boundary_width(int L, int J) {
  return ((1 << J) - 1) * (L - 1);
}

int effective_filter_width(int L, int J) {
  return ((1 << J) - 1) * (L - 1) + 1;
}

int choose_pad(
  int n,
  int L,
  int J,
  double safety = 1.25
) {
  int b = boundary_width(L, J);
  int pad = std::max(1, static_cast<int>(std::ceil(safety * b)));
  if (pad >= n) pad = n - 1;
  return pad;
}

MODWPTDeseasonResult
deseasonalise_modwpt(
  arma::mat &X,
  const arma::uword period,
  const std::vector<std::string> &filter_names,
  int J,
  double penalty_lambda,
  double acf_min,
  double dt
) {
  MODWPTDeseasonResult result;
  result.X = X;
  result.acf_original = arma::datum::nan;
  result.acf_adjusted = arma::datum::nan;
  result.acf_gain = arma::datum::nan;
  result.used_adjusted = false;
  result.seasonal_nodes = arma::imat(0, 2, arma::fill::zeros);
  result.selected_filter = "";
  result.filter_width_penalty = arma::datum::nan;
  result.objective_score = -arma::datum::inf;
  result.skipped_low_acf = false;

  if (period <= 0.0) throw std::invalid_argument("period must be positive");
  if (dt <= 0.0) throw std::invalid_argument("dt must be positive");
  if ((1.0 / period) > (1.0 / (2.0 * dt))) {
    throw std::invalid_argument(
        "period is shorter than the Nyquist-resolvable period 2 * dt");
  }
  if (penalty_lambda < 0.0) {
    throw std::invalid_argument("penalty_lambda must be non-negative");
  }
  if (!std::isfinite(acf_min) || acf_min < 0.0 || acf_min > 1.0) {
    throw std::invalid_argument("acf_min must be finite and between 0 and 1");
  }

  const int N = static_cast<int>(X.n_cols);
  result.acf_original = seasonal_structure_score(X, period);
  if (std::isfinite(result.acf_original) && result.acf_original < acf_min) {
    result.acf_adjusted = result.acf_original;
    result.acf_gain = 0.0;
    result.filter_width_penalty = 0.0;
    result.objective_score = 0.0;
    result.skipped_low_acf = true;
    return result;
  }

  std::vector<WaveletFilterCandidate> filter_candidates;
  int min_filter_width = std::numeric_limits<int>::max();
  bool skipped_candidate_filter = false;
  for (const std::string &candidate_name :
       filter_candidates_from_names(filter_names)) {
    WaveletFilterCandidate candidate;
    candidate.name = candidate_name;
    candidate.filters = wavelet_filters_modwt(candidate_name);
    candidate.length = static_cast<int>(candidate.filters(0).n_elem);

    const int J_max = compute_J_max(N, candidate.length);
    candidate.J_eff = J;
    if (candidate.J_eff <= 0) {
      candidate.J_eff = J_max;
    }
    if (candidate.J_eff > J_max) candidate.J_eff = J_max;
    if (candidate.J_eff < 1) {
      skipped_candidate_filter = true;
      continue;
    }

    candidate.filter_width =
        effective_filter_width(candidate.length, candidate.J_eff);
    filter_candidates.push_back(candidate);
    min_filter_width = std::min(min_filter_width, candidate.filter_width);
  }

  if (filter_candidates.empty()) {
    if (skipped_candidate_filter) {
      throw std::invalid_argument("No wavelet filter candidate is valid for "
                                  "the requested series length");
    }
    throw std::invalid_argument("No wavelet filter candidate was evaluated");
  }

  double best_objective = -arma::datum::inf;
  bool found_candidate = false;
  bool found_accepted_candidate = false;

  for (const WaveletFilterCandidate& filter_candidate : filter_candidates) {
    const arma::vec& h = filter_candidate.filters(0);
    const arma::vec& g = filter_candidate.filters(1);
    const int L = filter_candidate.length;
    const int J_eff = filter_candidate.J_eff;

    arma::field<arma::imat> node_candidates =
      make_objective_node_candidates(period, dt, J_eff);

    const arma::uword pad = static_cast<arma::uword>(choose_pad(N, L, J_eff));
    arma::mat Xp = pad_reflect(X, pad);
    MODWPTree wt = modwpt_decompose(Xp, J_eff, h, g);
    const double support_penalty =
      (static_cast<double>(filter_candidate.filter_width) -
       static_cast<double>(min_filter_width)) /
      static_cast<double>(N);

    for (arma::uword i = 0; i < node_candidates.n_elem; ++i) {
      arma::mat seasonal_i = modwpt_reconstruct_mixed_nodes(
        wt, node_candidates(i), h, g);
      arma::mat Yp_i = Xp - seasonal_i;
      arma::mat Y_i = Yp_i.cols(pad, pad + X.n_cols - 1);
      const double acf_adjusted_i = seasonal_structure_score(Y_i, period);
      const double acf_gain_i = result.acf_original - acf_adjusted_i;
      const double objective_i =
        acf_gain_i - penalty_lambda * support_penalty;
      const bool accepted_i =
        std::isfinite(objective_i) && objective_i > 0.0;
      const double objective_tol = std::isfinite(best_objective) ?
        std::max(1e-14, std::abs(best_objective) * 1e-12) : 1e-14;

      if (!found_candidate ||
          (accepted_i && !found_accepted_candidate) ||
          (!std::isfinite(best_objective) && std::isfinite(objective_i)) ||
          (accepted_i == found_accepted_candidate &&
           objective_i > best_objective + objective_tol)) {
        best_objective = objective_i;
        result.X = accepted_i ? Y_i : X;
        result.acf_adjusted = acf_adjusted_i;
        result.used_adjusted = accepted_i;
        result.seasonal_nodes = node_candidates(i);
        result.selected_filter = filter_candidate.name;
        result.acf_gain = acf_gain_i;
        result.filter_width_penalty = support_penalty;
        result.objective_score = objective_i;
        found_candidate = true;
        found_accepted_candidate = found_accepted_candidate || accepted_i;
      }
    }
  }

  if (!found_candidate) {
    if (skipped_candidate_filter) {
      throw std::invalid_argument("No wavelet filter candidate is valid for "
                                  "the requested series length");
    }
    throw std::invalid_argument("No wavelet filter candidate was evaluated");
  }

  return result;
}
