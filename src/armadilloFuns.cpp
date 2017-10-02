#include <RcppArmadillo.h>
#include <iostream>
using namespace Rcpp;

// Armadillo functions used to run a fast Gibbs sampler.
//
//


// [[Rcpp::depends(RcppArmadillo)]]

namespace arf{

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
    }

    if(x < 0){
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
  }

  double cdfn(double x, double mu, double sigma){
    return 0.5 * (1 + arf::erf((x - mu) / (sigma * sqrt(2.0)))) ;
  }

NumericVector runifArm(int n_) {
  arma::vec unis = arma::randu(n_) ;

  return wrap(unis) ;
}

NumericVector rnormArm(int n_) {
  arma::vec norms = arma::randn(n_) ;

  return wrap(norms) ;
}

} // end arf namespace

arma::vec resnArm(int n_, arma::vec xi,
                  arma::vec omega, arma::vec alpha,
                   NumericVector tau) {
  //double alph = Rcpp::as<double>(alpha) ;
  arma::vec delta = alpha / sqrt(1 + alpha % alpha) ;
  arma::vec lb = pnorm(-1 * tau) ;
  arma::vec uni = arma::randu(n_) % (1 - lb) + lb ;
  NumericVector unir = wrap(uni) ;
  arma::vec truncN = qnorm(unir) ;
  arma::vec z_ = delta % truncN + sqrt(1 - (delta % delta)) % arma::randn(n_) ;
  arma::vec esn = xi + omega % z_ ;

  return esn;
}


//' Generate a random number from an extended skew normal distribution
//'
//' @param n_ The number of random variates to create
//' @param xi esn parameter
//' @param omega esn parameter
//' @param alpha esn parameter
//' @param tau esn parameter
//'
//' @export
// [[Rcpp::export]]
NumericVector resn(int n_, NumericVector xi,
                   NumericVector omega, NumericVector alpha,
                   NumericVector tau) {
  //double alph = Rcpp::as<double>(alpha) ;
  NumericVector delta = alpha / sqrt(1 + alpha * alpha) ;
  NumericVector lb = pnorm(-1 * tau) ;
  NumericVector uni = runif(n_) * (1 - lb) + lb ;
  NumericVector truncN = qnorm(uni) ;
  NumericVector z_ = delta * truncN + sqrt(1 - (delta * delta)) * rnorm(n_) ;
  NumericVector esn = xi + omega * z_ ;

  return esn;
}

//' @export
// [[Rcpp::export]]
NumericVector resn2(int &n_, NumericVector &xi,
                    NumericVector &omega, NumericVector &alpha,
                    double &tau) {
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

// [[Rcpp::export]]
NumericVector testfn(NumericVector x){
  int N_ = x.size() ;
  NumericVector ps(N_) ;
  for (int i = 0; i < N_; i++){
    ps[i] = arf::cdfn(x[i], 0, 1) ;
  }
  return ps ;
}

