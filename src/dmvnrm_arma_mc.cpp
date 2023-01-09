// [[Rcpp::depends("RcppArmadillo")]]

#include "dmvnrm_arma_mc.h"

static double const log2pi = std::log(2.0 * M_PI);

void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
      x[j] = tmp;
  }
}

// [[Rcpp::export]]
double log_dmvnrm_arma_regular(arma::mat const &data_x,  
                                arma::rowvec const &mean,  
                                arma::mat const &prec) {  
  int n         = data_x.n_rows;
  int p         = data_x.n_cols;
  double log_ll = 0.0, log_det_term, sign;
  arma::rowvec data_shift;
  arma::log_det(log_det_term, sign, prec);
  log_ll       += -0.5 * (p*n) * log(2*M_PI);
  log_ll       += 0.5 * (n) * log_det_term;
  for(int i = 0; i < n; i++){
    data_shift  = data_x.row(i) - mean;
    log_ll     += -0.5 * arma::as_scalar( data_shift * prec * (data_shift.t()) );
  }
  return log_ll;
}

