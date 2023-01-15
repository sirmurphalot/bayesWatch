#include <omp.h>
#include <math.h>
#include <RcppArmadillo.h>
#include <boost/math/special_functions/erf.hpp>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

double log_mvn_pdf(const arma::vec& data_value, const arma::vec& conditional_mean, 
                   const arma::mat& conditional_covariance, const arma::mat& conditional_precision, const int p) {
  arma::vec deviation              = data_value - conditional_mean;
  arma::rowvec deviation_transpose = deviation.t();
  arma::mat exp_term               = deviation_transpose * conditional_precision * deviation;
  double double_exp_term           = (double)arma::as_scalar(exp_term);
  double double_det_term, sign;
  arma::log_det(double_det_term,sign, conditional_covariance);
  double log_pdf                   = - ( (double)p * 0.5 ) * log(2.0 * M_PI) * (-0.5) * double_det_term - 0.5 * double_exp_term;
  return(log_pdf);
}

arma::vec draw_values_parallel(const arma::vec& data_vector, const arma::mat& precision_matrix, const arma::vec& mu_vector, 
                               const arma::vec& lower_bounds, const arma::vec& upper_bounds, 
                               const arma::vec& lower_bound_is_equal, const arma::vec& upper_bound_is_equal, 
                               const arma::vec& row_is_missing, const arma::vec& is_continuous, int p, 
                               const arma::vec& raw_data_vector) {
  arma::vec data_vector_copy = data_vector;
  double sqrt_2              = sqrt(2.0);
  double zero                = 0.0;
  // Determine which values are missing for this row, subset the appropriate precision matrix, and redraw them:
  int total_number_missing = (int)accu(row_is_missing), count_redraw = 0, count_nonredraw = 0, 
    total_number_nonmissing = p - total_number_missing, discrete_val;
  if(total_number_missing > 0){
    arma::uvec location_of_missings(total_number_missing);
    arma::uvec location_of_nonmissings(total_number_nonmissing);
    for( int index = 0; index < p; index++ ){
      if( (int)row_is_missing.at(index) ){
        location_of_missings(count_redraw)       = index;
        count_redraw                            += 1;
      } else {
        location_of_nonmissings(count_nonredraw) = index;
        count_nonredraw                         += 1;
      }
    }
    // To draw the latent data, condition on all other variables, and draw values from the corresponding conditional MVN distribution.
    arma::mat submatrix_unknown            = precision_matrix(location_of_missings, location_of_missings);
    arma::mat submatrix_unknown_known      = precision_matrix(location_of_missings, location_of_nonmissings);
    arma::vec unknown_mu                   = mu_vector(location_of_missings);
    arma::vec known_mu                     = mu_vector(location_of_nonmissings);
    arma::mat conditional_covariance       = arma::inv(submatrix_unknown);
    arma::vec conditional_mean             = unknown_mu - conditional_covariance * submatrix_unknown_known * (data_vector(location_of_nonmissings) - known_mu);
    arma::vec std_norm_vals                = arma::randn<arma::vec>(unknown_mu.n_elem);
    arma::mat chol_decomp                  = arma::chol(conditional_covariance , "lower");
    arma::vec latent_vars_raw              = conditional_mean + chol_decomp * std_norm_vals;
    // log_proposal_density                  += log_mvn_pdf(latent_vars_raw, conditional_mean, conditional_covariance, submatrix_unknown, p);
    data_vector_copy(location_of_missings) = latent_vars_raw;
  }
  
  // Redraw all the discrete values and censored values (at boundaries) using truncated normals:
  arma::uvec nondiscrete_data_indicator(p);
  arma::uvec discrete_data_indicator(p);
  for( int index = 0; index < p; index++){
    nondiscrete_data_indicator(index) = 1;
    discrete_data_indicator(index)    = 0;
  }
  double conditional_covariance_scalar = 0.0, conditional_mean_scalar = 0.0, rand_value = 0.0, 
    lower_tail_probability = 0.0, upper_tail_probability = 0.0, rand_unif = 0.0, cdf_value = 0.0, 
    quantile_input = 0.0, epsilon = 1e-5;
  for(int index = 0; index < p; index++){
    
    // ---------------- If the variable is discrete. ---------------- //
    if( ((int)is_continuous.at(index) == 0) & ((int)row_is_missing.at(index) == 0) ){
      nondiscrete_data_indicator(index) = 0;
      discrete_data_indicator(index)    = 1;
      
      arma::uvec discrete_data    = arma::find(discrete_data_indicator);
      arma::uvec nondiscrete_data = arma::find(nondiscrete_data_indicator);
      
      // Gather all information for the conditional multivariate normal distribution.
      // To draw the latent data, condition on all other variables, and draw values from the corresponding conditional MVN distribution.
      double submatrix_unknown               = precision_matrix.at(index, index);
      arma::rowvec submatrix_12_full_row     = precision_matrix.row(index);
      arma::rowvec submatrix_unknown_known   = submatrix_12_full_row(nondiscrete_data).t();
      double unknown_mu                      = mu_vector.at(index);
      arma::vec known_mu                     = mu_vector(nondiscrete_data);
      double conditional_covariance_scalar   = pow(precision_matrix.at(index, index), -0.5);
      
      double conditional_mean_scalar         = unknown_mu - pow(submatrix_unknown, -1.0 ) * arma::as_scalar((submatrix_unknown_known * (data_vector_copy(nondiscrete_data) - known_mu)) );

      double cdf_value                       = arma::normcdf(zero, conditional_mean_scalar, conditional_covariance_scalar);
      rand_unif                              = arma::randu<double>();
      
      // If the value is ZERO, draw from a truncated normal beneath zero.
      discrete_val = (int)raw_data_vector.at(index);
      if( discrete_val == 0){
        lower_tail_probability  = (cdf_value) * rand_unif;
        quantile_input          = 2.0*(lower_tail_probability)-1.0;
        if(lower_tail_probability <= epsilon){
          rand_value            = -4.2649;
        } else if(lower_tail_probability >= 1 - epsilon ){
          rand_value            = 4.2649;
        } else {
          rand_value            = sqrt_2 * boost::math::erf_inv(quantile_input);
        }
        // rand_value              = sqrt_2 * boost::math::erf_inv(quantile_input);
        
        // If the value is ONE, draw from a truncated normal above zero
      } else {
        upper_tail_probability  = cdf_value + (1.0 - cdf_value) * rand_unif;
        quantile_input          = 2.0 * (upper_tail_probability) - 1.0;
        if(upper_tail_probability <= epsilon){
          rand_value            = -4.2649;
        } else if( upper_tail_probability >= 1 - epsilon ){
          rand_value            = 4.2649;
        } else {
          rand_value            = sqrt_2 * boost::math::erf_inv(quantile_input);
        }
        // rand_value              = sqrt_2 * boost::math::erf_inv(quantile_input);
      }
      
      // Use this random value to draw from the conditional, truncated normal distribution.
      data_vector_copy(index)       = rand_value*conditional_covariance_scalar + conditional_mean_scalar;

      nondiscrete_data_indicator(index) = 1;
      discrete_data_indicator(index)    = 0;
      
      // ---------------- If the variable is continuous, with a (realized) lower bound. ---------------- //
    } else if( ( (int)is_continuous.at(index) ) & ( (int)lower_bound_is_equal.at(index) ) ) {
      discrete_data_indicator(index) = 0;
      
      conditional_covariance_scalar = pow(precision_matrix.at(index, index), -0.5);
      arma::vec mean_temp_1         = (data_vector_copy % discrete_data_indicator - 
        mu_vector % discrete_data_indicator);
      arma::rowvec mean_temp_2      =  (precision_matrix.row(index) % discrete_data_indicator.t());
      conditional_mean_scalar       = mu_vector.at(index) + conditional_covariance_scalar * 
        arma::as_scalar( mean_temp_2 * mean_temp_1 ) ;
      
      rand_unif               = arma::randu<double>();
      cdf_value               = arma::normcdf(lower_bounds(index), conditional_mean_scalar, conditional_covariance_scalar);
      lower_tail_probability  = (cdf_value) * rand_unif;
      
      // To handle extreme edge cases:
      quantile_input          = 2.0*(lower_tail_probability)-1.0;
      if(lower_tail_probability <= epsilon){
        rand_value            = -4.2649;
      } else if(lower_tail_probability >= 1 - epsilon ){
        rand_value            = 4.2649;
      } else {
        rand_value            = sqrt_2 * boost::math::erf_inv(quantile_input);
      }
      
      //
      data_vector_copy(index) = conditional_mean_scalar + conditional_covariance_scalar * rand_value;
      
      // log_proposal_density   += arma::log_normpdf(data_vector_copy(index), conditional_mean_scalar, conditional_covariance_scalar) -log(cdf_value);
      
      discrete_data_indicator(index) = 1;
      
      // ---------------- If the variable is continuous, with a (realized) upper bound. ---------------- //
    } else if( ( (int)is_continuous.at(index) ) & ( (int)upper_bound_is_equal.at(index) ) ) {
      discrete_data_indicator(index) = 0;
      
      conditional_covariance_scalar = pow(precision_matrix.at(index, index), -0.5);
      arma::vec mean_temp_1         = (data_vector_copy % discrete_data_indicator - 
        mu_vector % discrete_data_indicator);
      arma::rowvec mean_temp_2      =  (precision_matrix.row(index) % discrete_data_indicator.t());
      conditional_mean_scalar       = mu_vector.at(index) + conditional_covariance_scalar * 
        arma::as_scalar( mean_temp_2 * mean_temp_1 ) ;
      
      rand_unif               = arma::randu<double>();
      cdf_value               = arma::normcdf(upper_bounds(index), conditional_mean_scalar, conditional_covariance_scalar);
      upper_tail_probability  = (1.0 - cdf_value) * rand_unif;
      // To handle extreme edge cases:
      quantile_input          = 2.0 * (1.0 - upper_tail_probability) - 1.0;
      if(upper_tail_probability <= epsilon){
        rand_value            = -4.2649;
      } else if( upper_tail_probability >= 1 - epsilon ){
        rand_value            = 4.2649;
      } else {
        rand_value            = sqrt_2 * boost::math::erf_inv(quantile_input);
      }
      //
      data_vector_copy(index) = conditional_mean_scalar + conditional_covariance_scalar * rand_value;
      
      // log_proposal_density   += arma::log_normpdf(data_vector_copy(index), conditional_mean_scalar, conditional_covariance_scalar) -log(1.0 - cdf_value);
      discrete_data_indicator(index) = 1;
    }
  }
  // Return the redrawn row.
  return(data_vector_copy);
}


//' Redraw latent data.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
arma::mat redraw_Z_arma(const arma::mat& current_data, const arma::mat& current_precision, const arma::vec& current_mu, int p,
                        const arma::vec& lower_bounds, const arma::vec& upper_bounds, 
                        const arma::mat& lower_bound_is_equal, const arma::mat& upper_bound_is_equal, const arma::mat& is_missing, 
                        const arma::vec& is_continuous, const arma::mat& raw_data, int cores) {
//  omp_set_num_threads( cores );
  int n                             = current_data.n_rows;
  arma::mat new_data                = current_data;
  // arma::vec log_proposal_density_values(n);

   #pragma omp parallel
   {
     #pragma omp for
      for( int index = 0; index < n; index++){
        arma::vec grab_col                 = new_data.row(index).t();
        arma::vec lower_bound_is_equal_row = lower_bound_is_equal.row(index).t();
        arma::vec upper_bound_is_equal_row = upper_bound_is_equal.row(index).t();
        arma::vec row_is_missing           = is_missing.row(index).t();
        arma::vec row_raw_data             = raw_data.row(index).t();
        
        // double temp_log_proposal_density   = 0.0;
        arma::vec drawn_data_temp          = draw_values_parallel(grab_col, current_precision, current_mu, 
                                                                  lower_bounds, upper_bounds, 
                                                                  lower_bound_is_equal_row, upper_bound_is_equal_row, 
                                                                  row_is_missing, is_continuous, p, row_raw_data);
        new_data.row(index)                = drawn_data_temp.t();
        // log_proposal_density_values(index) = temp_log_proposal_density;
      }
   }
  //  double log_proposal_density = arma::accu(log_proposal_density_values);
  
  // check that things work, then code up the MH step.
  
  return(new_data);
}

//' Redraw latent data for missing values only (not including for discrete values).
//' 
//' @noRd
//' 
// [[Rcpp::export]]
List redraw_Z_arma_justmissings(const arma::mat& current_data, const arma::mat& current_precision, const arma::vec& current_mu, int p,
                                const arma::vec& lower_bounds, const arma::vec& upper_bounds, 
                                const arma::mat& lower_bound_is_equal, const arma::mat& upper_bound_is_equal, const arma::mat& is_missing, 
                                const arma::vec& is_continuous, int cores) {
  //  omp_set_num_threads( cores );
  int n                             = current_data.n_rows;
  arma::mat new_data                = current_data;
  arma::vec log_proposal_density_values(n);
  List return_items;
   #pragma omp parallel
   {
   #pragma omp for
  for( int index = 0; index < n; index++){
    arma::vec grab_col                 = new_data.row(index).t();
    arma::vec data_vector_copy         = grab_col;
    arma::vec lower_bound_is_equal_row = lower_bound_is_equal.row(index).t();
    arma::vec upper_bound_is_equal_row = upper_bound_is_equal.row(index).t();
    arma::vec row_is_missing           = is_missing.row(index).t();
    
    double temp_log_proposal_density   = 0.0;
    double log_proposal_density   = 0.0;
    
    
    int total_number_missing = (int)accu(row_is_missing), count_redraw = 0, count_nonredraw = 0, total_number_nonmissing = p - total_number_missing;
    if(total_number_missing > 0){
      arma::uvec location_of_missings(total_number_missing);
      arma::uvec location_of_nonmissings(total_number_nonmissing);
      for( int index = 0; index < p; index++ ){
        if( (int)row_is_missing.at(index) ){
          location_of_missings(count_redraw)       = index;
          count_redraw                            += 1;
        } else {
          location_of_nonmissings(count_nonredraw) = index;
          count_nonredraw                         += 1;
        }
      }
      arma::mat submatrix_unknown            = current_precision(location_of_missings, location_of_missings);
      arma::mat submatrix_unknown_known      = current_precision(location_of_missings, location_of_nonmissings);
      arma::vec unknown_mu                   = current_mu(location_of_missings);
      arma::vec known_mu                     = current_mu(location_of_nonmissings);
      arma::mat conditional_covariance       = arma::inv(submatrix_unknown);
      arma::vec conditional_mean             = unknown_mu - conditional_covariance * submatrix_unknown_known * (grab_col(location_of_nonmissings) - known_mu);
      arma::vec std_norm_vals                = arma::randn<arma::vec>(unknown_mu.n_elem);
      arma::mat chol_decomp                  = arma::chol(conditional_covariance , "lower");
      arma::vec latent_vars_raw              = conditional_mean + chol_decomp * std_norm_vals;
      log_proposal_density                  += log_mvn_pdf(latent_vars_raw, conditional_mean, conditional_covariance, submatrix_unknown, p);
      data_vector_copy(location_of_missings) = latent_vars_raw;
    }
    
    new_data.row(index)                = data_vector_copy.t();
    log_proposal_density_values(index) = temp_log_proposal_density;
  }
   }
  double log_proposal_density = arma::accu(log_proposal_density_values);
  return_items["new_data"] = new_data;
  return_items["log_proposal_density"] = log_proposal_density;
  return return_items;
}


//' Fast evaluation of multivariate normal density of latent data values corresponding to missing values.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
List get_justmissings_density(const arma::mat& current_data, const arma::mat& current_precision, const arma::vec& current_mu, int p,
                                const arma::vec& lower_bounds, const arma::vec& upper_bounds, 
                                const arma::mat& lower_bound_is_equal, const arma::mat& upper_bound_is_equal, const arma::mat& is_missing, 
                                const arma::vec& is_continuous, int cores) {
  //  omp_set_num_threads( cores );
  int n                             = current_data.n_rows;
  arma::mat new_data                = current_data;
  arma::vec log_proposal_density_values(n);
  List return_items;
  #pragma omp parallel
  {
  #pragma omp for
    for( int index = 0; index < n; index++){
      arma::vec grab_col                 = new_data.row(index).t();
      arma::vec data_vector_copy         = grab_col;
      arma::vec lower_bound_is_equal_row = lower_bound_is_equal.row(index).t();
      arma::vec upper_bound_is_equal_row = upper_bound_is_equal.row(index).t();
      arma::vec row_is_missing           = is_missing.row(index).t();
      
      double temp_log_proposal_density   = 0.0;
      double log_proposal_density   = 0.0;
      
      
      int total_number_missing = (int)accu(row_is_missing), count_redraw = 0, count_nonredraw = 0, total_number_nonmissing = p - total_number_missing;
      if(total_number_missing > 0){
        arma::uvec location_of_missings(total_number_missing);
        arma::uvec location_of_nonmissings(total_number_nonmissing);
        for( int index = 0; index < p; index++ ){
          if( (int)row_is_missing.at(index) ){
            location_of_missings(count_redraw)       = index;
            count_redraw                            += 1;
          } else {
            location_of_nonmissings(count_nonredraw) = index;
            count_nonredraw                         += 1;
          }
        }
        arma::mat submatrix_unknown            = current_precision(location_of_missings, location_of_missings);
        arma::mat submatrix_unknown_known      = current_precision(location_of_missings, location_of_nonmissings);
        arma::vec unknown_mu                   = current_mu(location_of_missings);
        arma::vec known_mu                     = current_mu(location_of_nonmissings);
        arma::mat conditional_covariance       = arma::inv(submatrix_unknown);
        arma::vec conditional_mean             = unknown_mu - conditional_covariance * submatrix_unknown_known * (grab_col(location_of_nonmissings) - known_mu);
        arma::vec std_norm_vals                = arma::randn<arma::vec>(unknown_mu.n_elem);
        arma::mat chol_decomp                  = arma::chol(conditional_covariance , "lower");
        arma::vec latent_vars_raw              = data_vector_copy(location_of_missings);
        log_proposal_density                  += log_mvn_pdf(latent_vars_raw, conditional_mean, conditional_covariance, submatrix_unknown, p);
        data_vector_copy(location_of_missings) = latent_vars_raw;
      }
      
      new_data.row(index)                = data_vector_copy.t();
      log_proposal_density_values(index) = temp_log_proposal_density;
    }
  }
  double log_proposal_density = arma::accu(log_proposal_density_values);
  return_items["new_data"] = new_data;
  return_items["log_proposal_density"] = log_proposal_density;
  return return_items;
}



