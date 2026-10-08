#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <cmath>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

namespace {

const arma::uword S2_BANDS = 10;
const arma::uword LANDSAT_BANDS = 6;
const int SENSOR_S2 = 0;
const int SENSOR_LANDSAT = 1;

void sanitise_index_column(arma::subview_col<double>& output) {
  output.replace(arma::datum::inf, arma::datum::nan);
  output.replace(-arma::datum::inf, arma::datum::nan);
}

void write_s2_index_column(
  const arma::mat& input,
  const arma::uword source,
  const arma::uword index_id,
  const double scale,
  arma::mat& output,
  const arma::uword destination
) {
  arma::subview_col<double> result = output.col(destination);
  const arma::uword blue = source;
  const arma::uword green = source + 1;
  const arma::uword red = source + 2;
  const arma::uword re1 = source + 3;
  const arma::uword re2 = source + 4;
  const arma::uword re3 = source + 5;
  const arma::uword nir1 = source + 6;
  const arma::uword nir2 = source + 7;
  const arma::uword swir1 = source + 8;
  const arma::uword swir2 = source + 9;

  switch (index_id) {
    case 0:
      result = (input.col(nir1) - input.col(red)) /
        (input.col(nir1) + input.col(red));
      break;
    case 1:
      result = (input.col(nir2) - input.col(re1)) /
        (input.col(nir2) + input.col(re1));
      break;
    case 2:
      result = (input.col(nir2) - input.col(re2)) /
        (input.col(nir2) + input.col(re2));
      break;
    case 3:
      result = (input.col(nir2) - input.col(re3)) /
        (input.col(nir2) + input.col(re3));
      break;
    case 4:
      result = (input.col(nir2) - input.col(swir1)) /
        (input.col(nir2) + input.col(swir1));
      break;
    case 5:
      result = (input.col(nir2) - input.col(swir2)) /
        (input.col(nir2) + input.col(swir2));
      break;
    case 6:
      result = input.col(swir1) / input.col(nir2);
      break;
    case 7:
      result = scale * (
        0.3510 * input.col(blue) + 0.3813 * input.col(green) +
        0.3437 * input.col(red) + 0.7196 * input.col(nir1) +
        0.2396 * input.col(swir1) + 0.1949 * input.col(swir2)
      );
      break;
    case 8:
      result = scale * (
        -0.3599 * input.col(blue) - 0.3533 * input.col(green) -
        0.4734 * input.col(red) + 0.6633 * input.col(nir1) +
        0.0087 * input.col(swir1) - 0.2856 * input.col(swir2)
      );
      break;
    case 9:
      result = scale * (
        0.2578 * input.col(blue) + 0.2305 * input.col(green) +
        0.0883 * input.col(red) + 0.1071 * input.col(nir1) -
        0.7611 * input.col(swir1) - 0.5308 * input.col(swir2)
      );
      break;
    case 10: {
      arma::vec tcb =
        0.3510 * input.col(blue) + 0.3813 * input.col(green) +
        0.3437 * input.col(red) + 0.7196 * input.col(nir1) +
        0.2396 * input.col(swir1) + 0.1949 * input.col(swir2);
      arma::vec tcg =
        -0.3599 * input.col(blue) - 0.3533 * input.col(green) -
        0.4734 * input.col(red) + 0.6633 * input.col(nir1) +
        0.0087 * input.col(swir1) - 0.2856 * input.col(swir2);
      result = arma::atan(tcg / tcb);
      break;
    }
    case 11:
      result = input.col(swir1) / (
        input.col(nir2) + ((input.col(swir2) - input.col(nir2)) /
        (2190.0 - 865.0)) * (1610.0 - 865.0)
      );
      break;
    case 12:
      result = input.col(re1) / (
        input.col(red) + ((input.col(re2) - input.col(red)) /
        (740.0 - 665.0)) * (705.0 - 665.0)
      );
      break;
    case 13:
      result = scale * arma::sqrt(
        arma::square(input.col(red)) + arma::square(input.col(swir2))
      );
      break;
    case 14:
      result = scale * (input.col(re3) - input.col(red)) /
        (input.col(re1) / input.col(re2));
      break;
    case 15:
      result = (input.col(nir2) -
        (input.col(swir1) - input.col(swir2))) /
        (input.col(nir2) +
        (input.col(swir1) - input.col(swir2)));
      break;
  }
  sanitise_index_column(result);
}

void write_landsat_index_column(
  const arma::mat& input,
  const arma::uword source,
  const arma::uword index_id,
  const double scale,
  arma::mat& output,
  const arma::uword destination
) {
  arma::subview_col<double> result = output.col(destination);
  const arma::uword blue = source;
  const arma::uword green = source + 1;
  const arma::uword red = source + 2;
  const arma::uword nir = source + 3;
  const arma::uword swir1 = source + 4;
  const arma::uword swir2 = source + 5;

  switch (index_id) {
    case 0:
      result = (input.col(nir) - input.col(red)) /
        (input.col(nir) + input.col(red));
      break;
    case 4:
      result = (input.col(nir) - input.col(swir1)) /
        (input.col(nir) + input.col(swir1));
      break;
    case 5:
      result = (input.col(nir) - input.col(swir2)) /
        (input.col(nir) + input.col(swir2));
      break;
    case 6:
      result = input.col(swir1) / input.col(nir);
      break;
    case 7:
      result = scale * (
        0.3037 * input.col(blue) + 0.2793 * input.col(green) +
        0.4743 * input.col(red) + 0.5585 * input.col(nir) +
        0.5082 * input.col(swir1) + 0.1863 * input.col(swir2)
      );
      break;
    case 8:
      result = scale * (
        -0.2848 * input.col(blue) - 0.2435 * input.col(green) -
        0.5436 * input.col(red) + 0.7243 * input.col(nir) +
        0.0840 * input.col(swir1) - 0.1800 * input.col(swir2)
      );
      break;
    case 9:
      result = scale * (
        0.1509 * input.col(blue) + 0.1973 * input.col(green) +
        0.3279 * input.col(red) + 0.3406 * input.col(nir) -
        0.7112 * input.col(swir1) - 0.4572 * input.col(swir2)
      );
      break;
    case 10: {
      arma::vec tcb =
        0.3037 * input.col(blue) + 0.2793 * input.col(green) +
        0.4743 * input.col(red) + 0.5585 * input.col(nir) +
        0.5082 * input.col(swir1) + 0.1863 * input.col(swir2);
      arma::vec tcg =
        -0.2848 * input.col(blue) - 0.2435 * input.col(green) -
        0.5436 * input.col(red) + 0.7243 * input.col(nir) +
        0.0840 * input.col(swir1) - 0.1800 * input.col(swir2);
      result = arma::atan(tcg / tcb);
      break;
    }
    case 13:
      result = scale * arma::sqrt(
        arma::square(input.col(red)) + arma::square(input.col(swir2))
      );
      break;
    case 15:
      result = (input.col(nir) -
        (input.col(swir1) - input.col(swir2))) /
        (input.col(nir) +
        (input.col(swir1) + input.col(swir2)));
      break;
  }
  sanitise_index_column(result);
}

void apply_raster_mask_columns(arma::mat& values,
                               const arma::vec& mask_values) {
  
  const arma::uvec masked_rows = arma::find_nan(mask_values);
  if (masked_rows.is_empty()) return;

  values.each_col([&masked_rows](arma::vec& column) {
    column.elem(masked_rows).fill(arma::datum::nan);
  });
}

// Weighted medoid implementation
arma::vec medoid(const arma::mat& x, const arma::vec& weights) {
  const arma::uword dimensions = x.n_rows;
  const arma::uword points = x.n_cols;
  if (points == 0) return arma::vec();
  if (points == 1) return x.col(0);

  arma::vec dimensional_median(dimensions);
  for (arma::uword row = 0; row < dimensions; ++row) {
    dimensional_median(row) = arma::median(x.row(row));
  }

  double minimum_distance = std::numeric_limits<double>::max();
  arma::uword initial_medoid = 0;
  for (arma::uword point = 0; point < points; ++point) {
    const double distance = arma::norm(x.col(point) - dimensional_median);
    if (distance < minimum_distance) {
      minimum_distance = distance;
      initial_medoid = point;
    }
  }

  double minimum_objective = 0.0;
  for (arma::uword other = 0; other < points; ++other) {
    minimum_objective += weights(other) *
      arma::norm(x.col(initial_medoid) - x.col(other));
  }
  arma::uword final_medoid = initial_medoid;

  for (arma::uword point = 0; point < points; ++point) {
    if (point == initial_medoid) continue;
    double objective = 0.0;
    for (arma::uword other = 0; other < points; ++other) {
      objective += weights(other) *
        arma::norm(x.col(point) - x.col(other));
    }
    if (objective < minimum_objective) {
      minimum_objective = objective;
      final_medoid = point;
    }
  }
  return x.col(final_medoid);
}

// Vardi-Zhang implementation of the geometric median
arma::vec geometric_median(
  const arma::mat& x,
  const arma::vec& weights,
  double tolerance = 1e-8,
  arma::uword max_iterations = 10000
) {
  const arma::uword dimensions = x.n_rows;
  const arma::uword points = x.n_cols;
  if (points == 0) return arma::vec(dimensions, arma::fill::zeros);

  arma::vec current(dimensions), next(dimensions);
  arma::vec distances(points);
  arma::vec inverse(points);
  arma::vec inverse_excluding_coincident(points);
  arma::vec numerator(dimensions);

  const double weight_sum = arma::sum(weights);
  current = x * (weights / weight_sum);
  const double tiny = std::max(1e-16, tolerance * 1e-8);
  const double scale = arma::max(arma::vecnorm(x, 2, 0));
  const double duplicate_tolerance = std::max(1e-12, 1e-12 * scale);
  const double safe_distance = std::max(tiny, duplicate_tolerance);

  arma::uword iteration = 0;
  double normalized_difference = 1.0;
  while (iteration < max_iterations &&
         normalized_difference > tolerance) {
    for (arma::uword point = 0; point < points; ++point) {
      distances(point) = arma::norm(x.col(point) - current);
    }

    const arma::uvec coincident = arma::find(distances <= safe_distance);
    if (!coincident.is_empty()) {
      const arma::uword point = coincident(0);
      inverse_excluding_coincident.zeros();
      const arma::uvec other = arma::find(distances > safe_distance);
      if (!other.is_empty()) {
        inverse_excluding_coincident.elem(other) =
          weights.elem(other) / distances.elem(other);
      }

      numerator = x * inverse_excluding_coincident;
      const double denominator = arma::sum(inverse_excluding_coincident);
      const arma::vec coincident_point = x.col(point);
      const arma::vec residual =
        numerator - coincident_point * denominator;
      const double residual_norm = arma::norm(residual);
      if (residual_norm <= weights(point)) {
        next = coincident_point;
      } else {
        const double safe_denominator = denominator <= tiny ?
          tiny : denominator;
        const double step =
          (residual_norm - weights(point)) / safe_denominator;
        next = coincident_point + step * residual / residual_norm;
      }
    } else {
      const arma::uvec other = arma::find(distances > safe_distance);
      inverse.zeros();
      if (!other.is_empty()) {
        inverse.elem(other) = weights.elem(other) / distances.elem(other);
      }
      const double denominator = arma::sum(inverse);
      numerator = x * inverse;
      const double safe_denominator = denominator <= tiny ?
        tiny : denominator;
      next = numerator / safe_denominator;
    }

    normalized_difference = arma::norm(next - current) /
      std::sqrt(static_cast<double>(dimensions));
    current = next;
    ++iteration;
  }
  return current;
}

arma::vec composite_location(
  const arma::mat& x,
  const arma::vec& weights,
  int method
) {
  return method == 0 ? medoid(x, weights) : geometric_median(x, weights);
}

arma::uword representative_observation(
  const arma::mat& x,
  const arma::vec& location,
  const arma::uvec& days,
  arma::uword target_day
) {
  arma::mat differences = x;
  differences.each_col() -= location;
  const arma::vec spectral_distances = arma::vecnorm(differences, 2, 0).t();
  arma::uvec candidates = arma::find(
    spectral_distances == spectral_distances.min()
  );
  if (candidates.n_elem == 1) return candidates(0);

  const arma::vec temporal_distances = arma::abs(
    arma::conv_to<arma::vec>::from(days) - static_cast<double>(target_day)
  );
  const arma::vec candidate_temporal = temporal_distances.elem(candidates);
  candidates = candidates.elem(arma::find(
    candidate_temporal == candidate_temporal.min()
  ));
  if (candidates.n_elem == 1) return candidates(0);

  return candidates(days.elem(candidates).index_min());
}

inline double sigmoid_steepness(
  double minimum,
  double median,
  double epsilon = 1e-3
) {
  if (std::abs(minimum - median) < 1e-15) return 1.0;
  return std::log(1.0 / (1.0 - epsilon) - 1.0) /
    (minimum - median);
}

arma::vec spectral_weights(const arma::vec& reference, const arma::mat& x) {
  arma::mat differences = x;
  differences.each_col() -= reference;
  const arma::vec euclidean = arma::vecnorm(differences, 2, 0).t();
  const arma::vec point_norms = arma::vecnorm(x, 2, 0).t();
  const double reference_norm = arma::norm(reference, 2);
  const arma::vec denominators = reference_norm * point_norms;
  arma::vec angle(x.n_cols, arma::fill::zeros);
  const arma::uvec nonzero = arma::find(denominators > 0.0);
  if (!nonzero.is_empty()) {
    arma::vec cosine = (x.cols(nonzero).t() * reference) /
      denominators.elem(nonzero);
    cosine = arma::clamp(cosine, -1.0, 1.0);
    angle.elem(nonzero) = arma::acos(cosine);
  }
  const double median_euclidean = arma::median(euclidean);
  const double median_angle = arma::median(angle);
  const double a_e = sigmoid_steepness(euclidean.min(), median_euclidean);
  const double a_a = sigmoid_steepness(angle.min(), median_angle);
  return 1.0 / (1.0 + arma::exp(a_e * (euclidean - median_euclidean))) +
    1.0 / (1.0 + arma::exp(a_a * (angle - median_angle)));
}

arma::vec observation_weights(
  const arma::vec& distances,
  const arma::vec& reference,
  const arma::mat& x,
  double requested_distance
) {
  arma::vec distance_weight(distances.n_elem, arma::fill::ones);
  if (requested_distance > 0.0) {
    distance_weight = 1.0 / (1.0 + arma::exp((distances -
      requested_distance / 2.0) * (-10.0 / requested_distance)));
  }
  arma::vec weights =
      arma::exp(distance_weight + spectral_weights(reference, x));
  const double total = arma::accu(weights);
  if (!std::isfinite(total) || total == 0.0) weights.ones();
  else weights /= total;
  return weights;
}

arma::uvec balanced_assignments(arma::uword observations, arma::uword bins) {
  arma::uvec assignment(observations);
  if (observations == 0) return assignment;

  const arma::uword base = observations / bins;
  const arma::uword remainder = observations % bins;
  const arma::uword extended = (base + 1) * remainder;

  if (extended > 0) {
    assignment.head(extended) =
      arma::regspace<arma::uvec>(0, extended - 1) / (base + 1);
  }
  if (extended < observations) {
    const arma::uword remaining = observations - extended;
    assignment.tail(remaining) = remainder +
      arma::regspace<arma::uvec>(0, remaining - 1) / base;
  }
  return assignment;
}

arma::vec composite_cell(
  arma::rowvec values,
  arma::uword spectral_bands,
  const arma::uvec& acquisition_days,
  const arma::uvec& acquisition_doys,
  const arma::uvec& acquisition_years,
  const arma::uvec& acquisition_sensors,
  const arma::umat& intervals,
  const arma::uvec& target_days,
  const arma::uvec& target_doys,
  arma::uword bins,
  int method,
  bool hybrid_binning,
  bool use_weights,
  double requested_distance,
  double hybrid_max_days,
  const unsigned int qai_mask
) {
  const arma::uword input_layers = spectral_bands + (use_weights ? 2 : 1);
  const arma::uword quality_layer = input_layers - 1;
  const arma::uword output_layers = spectral_bands + 6;
  const arma::uword observations = values.n_elem / input_layers;
  arma::mat input(values.memptr(), input_layers, observations, false, true);
  
  arma::vec output(output_layers * bins, arma::fill::value(arma::datum::nan));
  
  const arma::rowvec quality_values = input.row(quality_layer);
  arma::uvec valid_quality(observations, arma::fill::zeros);
  const arma::uvec finite_quality = arma::find_finite(quality_values);
  if (!finite_quality.is_empty()) {
    arma::uvec quality_bits = arma::conv_to<arma::uvec>::from(
      quality_values.elem(finite_quality)
    );
    quality_bits.transform([qai_mask](arma::uword value) {
      return (value & qai_mask) == 0u ? 1u : 0u;
    });
    valid_quality.elem(finite_quality) = quality_bits;
  }
  
  const arma::mat spectral_input = input.rows(0, spectral_bands - 1);
  arma::uvec finite_spectra(observations, arma::fill::ones);
  const arma::uvec nonfinite_spectra = arma::find_nonfinite(spectral_input);
  if (!nonfinite_spectra.is_empty()) {
    const arma::uvec invalid_columns = arma::unique(
      nonfinite_spectra / spectral_bands
    );
    finite_spectra.elem(invalid_columns).zeros();
  }
  
  const arma::uvec positions = arma::find(valid_quality % finite_spectra);
  if (positions.is_empty()) return output;
  
  const arma::mat spectra = spectral_input.cols(positions);
  arma::vec distances(positions.n_elem, arma::fill::zeros);
  if (use_weights) {
    const arma::rowvec distance_input = input.row(spectral_bands);
    distances = distance_input.elem(positions);
  }
  const arma::uvec days = acquisition_days.elem(positions);
  arma::vec unit_weights(spectra.n_cols, arma::fill::ones);
  const arma::vec reference = composite_location(spectra, unit_weights, method);
  const arma::vec medoid_reference = method == 0 ? reference :
    medoid(spectra, unit_weights);

  arma::uvec assignment(
      spectra.n_cols,
      arma::fill::value(std::numeric_limits<arma::uword>::max()));
  if (!hybrid_binning) {
    for (arma::uword bin = 0; bin < bins; ++bin) {
      assignment.elem(arma::find((days >= intervals(bin, 0)) &&
        (days <= intervals(bin, 1)))).fill(bin);
    }
  } else {
    assignment = balanced_assignments(spectra.n_cols, bins);
    const arma::vec centres = arma::mean(
      arma::conv_to<arma::mat>::from(intervals), 1);
    arma::mat temporal_distances(spectra.n_cols, bins);
    temporal_distances.each_col() = arma::conv_to<arma::vec>::from(days);
    temporal_distances.each_row() -= centres.t();
    temporal_distances = arma::abs(temporal_distances);
    
    const arma::uvec observation_indices =
      arma::regspace<arma::uvec>(0, spectra.n_cols - 1);
    const arma::uvec assigned_indices =
      assignment * spectra.n_cols + observation_indices;
    const arma::vec assigned_distances =
      temporal_distances.elem(assigned_indices);
    const arma::uvec move = arma::find(assigned_distances > hybrid_max_days);
    if (!move.is_empty()) {
      const arma::uvec closest = arma::index_min(temporal_distances, 1);
      assignment.elem(move) = closest.elem(move);
    }
  }
  
  for (arma::uword bin = 0; bin < bins; ++bin) {
    const arma::uvec selected = arma::find(assignment == bin);
    if (selected.is_empty()) continue;
    const arma::uword offset = bin * output_layers;
    const arma::mat bin_spectra = spectra.cols(selected);
    arma::vec weights(selected.n_elem, arma::fill::ones);
    if (use_weights) {
      weights = observation_weights(distances.elem(selected), reference,
                                    bin_spectra, requested_distance);
    }
    const arma::vec location = composite_location(bin_spectra, weights, method);
    output.subvec(offset, offset + spectral_bands - 1) = location;
    
    arma::vec medoid_weights(selected.n_elem, arma::fill::ones);
    if (use_weights) {
      medoid_weights = method == 0 ? weights : observation_weights(
        distances.elem(selected), medoid_reference, bin_spectra,
        requested_distance);
    }
    const arma::vec information_medoid = method == 0 ? location :
      medoid(bin_spectra, medoid_weights);
    
    const arma::uvec bin_days = days.elem(selected);
    const arma::uword representative = representative_observation(
      bin_spectra, information_medoid, bin_days, target_days(bin));
    const arma::uword selected_position = selected(representative);
    const arma::uword acquisition_position = positions(selected_position);
    const arma::uword information_offset = offset + spectral_bands;
    output(information_offset) = input(quality_layer,
           acquisition_position);                              // QAI
    output(information_offset + 1) = selected.n_elem;          // NOBS
    output(information_offset + 2) = acquisition_doys(
      acquisition_position);                                   // DOY
    output(information_offset + 3) = acquisition_years(
      acquisition_position);                                   // YEAR
    output(information_offset + 4) = static_cast<double>(
      acquisition_doys(acquisition_position)) -
        static_cast<double>(target_doys(bin));                 // D_TDOY
    output(information_offset + 5) = acquisition_sensors(
      acquisition_position);                                   // SENSOR
  }
  return output;
}

} // namespace

// [[Rcpp::export]]
arma::vec mask_chunk_values_cpp(
  arma::vec& values,
  const arma::vec& mask_values,
  const arma::uword input_rows,
  const arma::uword input_columns
) {
  if (values.n_elem != input_rows * input_columns ||
      mask_values.n_elem != input_rows) {
    Rcpp::stop("raster-mask chunk dimensions do not match input values");
  }

  arma::vec output_values = values;
  arma::mat output(
    output_values.memptr(),
    input_rows,
    input_columns,
    false,
    true
  );
  apply_raster_mask_columns(output, mask_values);
  return output_values;
}

// [[Rcpp::export]]
arma::vec force_map_chunk_values_cpp(
  arma::vec& values,
  const arma::vec& mask_values,
  const arma::uword input_rows,
  const arma::uword input_columns,
  const arma::uword time_steps,
  const arma::uvec& selected_band_positions,
  const arma::uvec& index_ids,
  const int sensor,
  const bool rescale_reflectance,
  const int num_threads
) {
#ifndef _OPENMP
  (void)num_threads;
#endif

  if (sensor != SENSOR_S2 && sensor != SENSOR_LANDSAT) {
    Rcpp::stop("sensor must identify Sentinel-2 or Landsat");
  }
  if (index_ids.n_elem > 0 && index_ids.max() > 15) {
    Rcpp::stop("spectral index IDs must be between 0 and 15");
  }

  const arma::uword spectral_bands =
    sensor == SENSOR_S2 ? S2_BANDS : LANDSAT_BANDS;
  if (input_columns != spectral_bands * time_steps ||
      values.n_elem != input_rows * input_columns) {
    Rcpp::stop("mapping chunk has invalid FORCE reflectance dimensions");
  }
  if (selected_band_positions.n_elem > 0 &&
      (selected_band_positions.min() < 1 ||
       selected_band_positions.max() > spectral_bands)) {
    Rcpp::stop("selected band positions are outside the FORCE band set");
  }

  const bool use_mask = mask_values.n_elem > 0;
  if (use_mask && mask_values.n_elem != input_rows) {
    Rcpp::stop("raster-mask chunk dimensions do not match input values");
  }

  const arma::uword variables_per_step =
    selected_band_positions.n_elem + index_ids.n_elem;
  arma::vec output_values(
    input_rows * variables_per_step * time_steps,
    arma::fill::value(arma::datum::nan)
  );
  const arma::mat input(
    values.memptr(),
    input_rows,
    input_columns,
    false,
    true
  );
  arma::mat output(
    output_values.memptr(),
    input_rows,
    variables_per_step * time_steps,
    false,
    true
  );
  const double scale = rescale_reflectance ? 0.0001 : 1.0;
  const arma::uword jobs = variables_per_step * time_steps;

#ifdef _OPENMP
  #pragma omp parallel for schedule(static) num_threads(num_threads)
#endif
  for (arma::uword job = 0; job < jobs; ++job) {
    const arma::uword destination = job;
    const arma::uword step = destination / variables_per_step;
    const arma::uword variable = destination % variables_per_step;
    const arma::uword source = step * spectral_bands;

    if (variable < selected_band_positions.n_elem) {
      output.col(destination) = input.col(
        source + selected_band_positions(variable) - 1
      );
      continue;
    }

    const arma::uword index = variable - selected_band_positions.n_elem;
    if (sensor == SENSOR_S2) {
      write_s2_index_column(
        input, source, index_ids(index), scale, output, destination
      );
    } else {
      write_landsat_index_column(
        input, source, index_ids(index), scale, output, destination
      );
    }
  }

  if (use_mask) {
    apply_raster_mask_columns(output, mask_values);
  }
  return output_values;
}

arma::mat force_composite_chunk_impl(
  arma::vec& values,
  const arma::uword input_rows,
  const arma::uword input_columns,
  const arma::uword spectral_bands,
  const arma::uvec& acquisition_days,
  const arma::uvec& acquisition_doys,
  const arma::uvec& acquisition_years,
  const arma::uvec& acquisition_sensors,
  const arma::umat& intervals,
  const arma::uvec& target_days,
  const arma::uvec& target_doys,
  const arma::uword bins,
  const int method,
  const bool hybrid_binning,
  const bool use_weights,
  const double requested_distance,
  const double hybrid_max_days,
  const int num_threads,
  const unsigned int qai_mask
) {
#ifndef _OPENMP
  (void)num_threads;
#endif

  if (spectral_bands != S2_BANDS && spectral_bands != LANDSAT_BANDS) {
    Rcpp::stop(
        "composite input must contain ten Sentinel-2 or six Landsat bands");
  }
  if (bins == 0) {
    Rcpp::stop("composite output must contain at least one bin");
  }
  if (qai_mask > 65535u) {
    Rcpp::stop("qai_mask must fit within 16 QAI bits");
  }
  const arma::uword input_layers = spectral_bands + (use_weights ? 2 : 1);
  if (input_columns != input_layers * acquisition_days.n_elem) {
    Rcpp::stop(
        "composite input has an invalid number of layers per acquisition");
  }
  if (acquisition_days.n_elem != acquisition_doys.n_elem ||
      acquisition_days.n_elem != acquisition_years.n_elem ||
      acquisition_days.n_elem != acquisition_sensors.n_elem) {
    Rcpp::stop("acquisition metadata must have one value per observation");
  }
  if (intervals.n_rows != bins || target_days.n_elem != bins ||
      target_doys.n_elem != bins) {
    Rcpp::stop("bin metadata must have one value per output bin");
  }
  const arma::mat input(
    values.memptr(),
    input_rows,
    input_columns,
    false,
    true
  );
  arma::mat output(input_rows, (spectral_bands + 6) * bins);

#ifdef _OPENMP
  #pragma omp parallel for schedule(static) num_threads(num_threads)
#endif
  for (arma::uword row = 0; row < input_rows; ++row) {
    output.row(row) = composite_cell(input.row(row), spectral_bands,
      acquisition_days, acquisition_doys, acquisition_years,
      acquisition_sensors, intervals, target_days, target_doys, bins,
      method, hybrid_binning, use_weights,
      requested_distance, hybrid_max_days, qai_mask).t();
  }
  return output;
}

// [[Rcpp::export]]
arma::mat force_composite_chunk_cpp(
  arma::vec& values,
  const arma::uword input_rows,
  const arma::uword input_columns,
  const arma::uword spectral_bands,
  const arma::uvec& acquisition_days,
  const arma::uvec& acquisition_doys,
  const arma::uvec& acquisition_years,
  const arma::uvec& acquisition_sensors,
  const arma::umat& intervals,
  const arma::uvec& target_days,
  const arma::uvec& target_doys,
  const arma::uword bins,
  const int method,
  const bool hybrid_binning,
  const bool use_weights,
  const double requested_distance,
  const double hybrid_max_days,
  const int num_threads,
  const unsigned int qai_mask = 799
) {
  return force_composite_chunk_impl(
    values, input_rows, input_columns, spectral_bands,
    acquisition_days, acquisition_doys, acquisition_years,
    acquisition_sensors, intervals, target_days, target_doys, bins,
    method, hybrid_binning, use_weights, requested_distance,
    hybrid_max_days, num_threads, qai_mask
  );
}
