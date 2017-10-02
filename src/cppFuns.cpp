#include <RcppArmadillo.h>
using namespace Rcpp;

// Rcpp functions used to run a fast Gibbs sampler.
//
//



//' @export
// [[Rcpp::export]]
NumericVector rnormCpp(int n_) {
  //double alph = Rcpp::as<double>(alpha) ;
  NumericVector norms = rnorm(n_) ;

  return norms ;
}

