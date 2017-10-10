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
double resn(const double xi,
                   const double omega, const double alpha,
                   const double tau) {
  RNGScope rngScope ;
  double delta = alpha / sqrt(1 + alpha * alpha) ;
  double lb = arf::cdfn(-1 * tau, 0, 1) ;

  arma::vec uniVec(1); uniVec.randu() ;
  double uni = uniVec(0) * (1 - lb) + lb ;
  NumericVector nvUni(1); nvUni[0] = uni ;
  NumericVector truncN_(1); truncN_ = qnorm(nvUni) ;
  double truncN = truncN_(0) ;

  arma::vec normVec(1); normVec.randn() ;
  double norm = normVec(0) ;
  double z_ = delta * truncN + sqrt(1 - (delta * delta)) * norm ;
  double esn = xi + omega * z_ ;

  return esn;
}

// Function to update missing values and mean parameters
arma::vec updateBlock(int index, int iter, NumericMatrix yMat_,
                 arma::mat &yMiss,
                 NumericMatrix xMat_, NumericMatrix pointers,
                 arma::mat &intercepts, arma::mat fcs,
                 arma::mat &peps, const double miss_a,
                 const double miss_b, const double sigma,
                 const double tau_int, const double tau_fc,
                 const double tau_pep){
  arma::mat yMat(yMat_.begin(), yMat_.nrow(), yMat_.ncol(), false) ;
  arma::mat xMat(xMat_.begin(), xMat_.nrow(), xMat_.ncol(), false) ;

  //get parameter vector
  int nPars = xMat.n_cols ;
  arma::vec theta(nPars); theta.zeros() ;
  for(int p = 0; p < nPars; p++){
    if(pointers(p, 0) == 1){
      theta(p) = intercepts(pointers(p, 1) - 1, iter) ;
    }else{
      if(pointers(p, 0) == 2){
        theta(p) = fcs(pointers(p, 1) - 1, iter) ;
      }else{
        theta(p) = peps(pointers(p, 1) - 1, iter) ;
      }
    }
  }

  //compute esn sampling parameters
  arma::vec xi = xMat * theta ;
  double omega = sqrt(sigma) ;
  double alpha = -1 * miss_b * omega ;

  //get outcome vector
  int nObs = yMat.n_rows ;
  arma::vec y_(nObs); y_.zeros() ;
  for(int j = 0; j < nObs; j++){

    if(yMat(j, 0) == 0){
      double temptau =
        (-1 * miss_a - miss_b * xi(j)) / sqrt(1 + miss_b * miss_b * sigma) ;

      y_(j) = arf::resn(xi(j), omega, alpha, temptau) ;
      yMiss(yMat(j, 1) - 1, iter + 1) = y_(j) ;
    }else{
      y_(j) = yMat(j, 0) ;
    }
  }

  return(y_) ;
  } // end updateBlock

} // end arf namespace




// Run a custom Gibbs sampler for the SMP model
// [[Rcpp::export()]]
List gibbsCpp(List y_list,
              NumericMatrix y_miss_,
              NumericMatrix r_obs_,
              List matList,
              List pointers,
              NumericMatrix intercepts_,
              NumericMatrix fcs_,
              NumericMatrix peps_,
              NumericMatrix miss_a_,
              NumericMatrix miss_b_,
              NumericMatrix sigma_,
              NumericMatrix tau_int_,
              NumericMatrix tau_fc_,
              NumericMatrix tau_pep_){

  int n_prot = matList.size() ;
  //convert to armadillo objects
  arma::mat y_miss(y_miss_.begin(), y_miss_.nrow(), y_miss_.ncol(), false) ;
  arma::mat r_obs(r_obs_.begin(), r_obs_.nrow(), r_obs_.ncol(), false) ;
  arma::mat intercepts(intercepts_.begin(), intercepts_.nrow(),
                       intercepts_.ncol(), false) ;
  arma::mat fcs(fcs_.begin(), fcs_.nrow(), fcs_.ncol(), false) ;
  arma::mat peps(peps_.begin(), peps_.nrow(), peps_.ncol(), false) ;
  arma::mat miss_a(miss_a_.begin(), miss_a_.nrow(), miss_a_.ncol(), false) ;
  arma::mat miss_b(miss_b_.begin(), miss_b_.nrow(), miss_b_.ncol(), false) ;
  arma::mat sigma(sigma_.begin(), sigma_.nrow(), sigma_.ncol(), false) ;
  arma::mat tau_int(tau_int_.begin(), tau_int_.nrow(), tau_int_.ncol(), false) ;
  arma::mat tau_fc(tau_fc_.begin(), tau_fc_.nrow(), tau_fc_.ncol(), false) ;
  arma::mat tau_pep(tau_pep_.begin(), tau_pep_.nrow(), tau_pep_.ncol(), false) ;

  //update missing values and mean parameters
  NumericMatrix y1 = y_list[0] ;
  NumericMatrix mat1 = matList[0] ;
  NumericMatrix point1 = pointers[0] ;
  arma::vec test = arf::updateBlock(0, 0, y1, y_miss,
                               mat1, point1,
                               intercepts, fcs,
                               peps, miss_a(0,0),
                               miss_b(0,0), sigma(0),
                               tau_int(0), tau_fc(0),
                               tau_pep(0)) ;

  return List::create(Named("test") = test, Named("PointerChange") = y_miss) ;

}

