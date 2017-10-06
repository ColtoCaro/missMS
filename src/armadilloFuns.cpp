#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// Armadillo functions used to run a fast Gibbs sampler.


namespace arf{

// Approximate the error function (used in the cdf of a normal)
  double erf(double x){
    if(x >= 0){
    double y = 1.0 / (1.0 + 0.3275911 * x) ;
    return 1 - (((((
      + 1.061405429 * y
      - 1.453152027) * y
      + 1.421413741) * y
      - 0.284496736) * y
      + 0.254829592) * y)
      * exp(-x * x) ;
    }else{
      double z = -1 * x ;
      double y = 1.0 / (1.0 + 0.3275911 * z) ;
      return -1 * (1 - (((((
          + 1.061405429 * y
                       - 1.453152027) * y
                       + 1.421413741) * y
                       - 0.284496736) * y
                       + 0.254829592) * y)
                       * exp(-z * z) ) ;
    }
  } // End erf()

//fast computation of the cdf of a normal(mu, sigma)
  double cdfn(double x, double mu, double sigma){
    return 0.5 * (1 + arf::erf((x - mu) / (sigma * sqrt(2.0)))) ;
  } //End cdfn

//Generate uniform(0,1) variate
NumericVector runifArm(int n_) {
  arma::vec unis = arma::randu(n_) ;

  return wrap(unis) ;
}

// Generate random normal(0,1) variate
NumericVector rnormArm(int n_) {
  arma::vec norms = arma::randn(n_) ;

  return wrap(norms) ;
}

// Generate a random number from an extended skew normal distribution
NumericVector resn(int &n_, NumericVector &xi,
                   NumericVector &omega, NumericVector &alpha,
                   double &tau) {
  RNGScope rngScope ;
  //double alph = Rcpp::as<double>(alpha) ;
  NumericVector delta = alpha / sqrt(1 + alpha * alpha) ;
  NumericVector lb(n_) ;
  for(int i = 0; i < n_; i++){
    lb[i] = arf::cdfn(-1 * tau, 0, 1) ;
  }

  NumericVector uni = arf::runifArm(n_) * (1 - lb) + lb ;
  NumericVector truncN = qnorm(uni) ;
  NumericVector z_ = delta * truncN + sqrt(1 - (delta * delta)) *
    arf::rnormArm(n_) ;
  NumericVector esn = xi + omega * z_ ;

  return esn;
}



} // end arf namespace




// Run a custom Gibbs sampler for the SMP model
// [[Rcpp::export()]]
void gibbsCpp(List y_list,
              arma::mat y_miss,
              arma::mat r_obs,
              List matList,
              List pointers,
              arma::mat intercepts,
              arma::mat fcs,
              arma::mat peps,
              arma::mat miss_a,
              arma::mat miss_b,
              arma::mat sigma,
              arma::mat tau_int,
              arma::mat tau_fc,
              arma::mat tau_pep){

  Rcpp::Rcout << y_miss(4, 0) << std::endl  ;

}

