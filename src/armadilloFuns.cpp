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

// Function for updating Gaussian mean parameters
double sampN(const double sumDiff, const double tau, const double beta,
             const double sigma, const int J){
  double condMean = (sigma * beta + tau * sumDiff) / (sigma + tau * J) ;
  double condSd = sqrt((tau * sigma) / (sigma + tau * J)) ;

  arma::vec normVec(1); normVec.randn() ;
  double norm = normVec(0) ;

  return(condMean + condSd * norm) ;
}

// Function to update missing values and mean parameters
arma::vec updateBlock(int index, int iter, NumericMatrix yMat_,
                 arma::mat &yMiss,
                 NumericMatrix xMat_, NumericMatrix pointers,
                 arma::mat &intercepts, arma::mat &fcs,
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
 // Rcout << "before " << yMat << " theta = " << theta << std::endl ;
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
//  Rcout << "after " << y_ << " theta = " << theta << std::endl ;
  //Now update the mean parameters
  arma::mat xxi(size(xMat)) ; xxi.zeros() ;
  arma::mat groupYs(nObs, nPars) ; groupYs.zeros() ;
  for(int rw = 0; rw < nObs; rw++){
    xxi.row(rw) = xMat.row(rw) * xi(rw) ;
    groupYs.row(rw) = xMat.row(rw) * y_(rw) ;
 }

  arma::mat addBack(nObs, nPars) ; addBack.zeros() ;
  for(int cl = 0; cl < nPars; cl++){
    addBack.col(cl) = xMat.col(cl) * theta(cl) ;
  }

  arma::mat newMat = groupYs - xxi + addBack ;
  arma::rowvec sumDiff = sum(newMat) ;
  arma::rowvec nYs = sum(xMat) ;

  arma::vec theta2(nPars) ; theta2.zeros() ;
  for(int p = 0; p < nPars; p++){
    if(pointers(p, 0) == 1){
      double newInt = sampN(sumDiff(p), tau_int, 0, sigma, nYs(p)) ;
      intercepts(pointers(p, 1) - 1, iter + 1) = newInt ;
      theta2(p) = newInt ;
    }else{
      if(pointers(p, 0) == 2){
        double newFc = sampN(sumDiff(p), tau_fc, 0, sigma, nYs(p)) ;
        fcs(pointers(p, 1) - 1, iter + 1) = newFc ;
        theta2(p) = newFc ;
      }else{
        double newPep = sampN(sumDiff(p), tau_pep, 0, sigma, nYs(p)) ;
        peps(pointers(p, 1) - 1, iter + 1) = newPep ;
        theta2(p) = newPep ;
      }
    }
  }

  arma::vec resid = y_ - xMat * theta2 ;
  return (resid);
  } // end updateBlock

double sampV(arma::vec parVec, double hyp, double parMean){
  int nPars = parVec.size() ;
  double postShape = hyp + nPars / 2 ;
  arma::vec meanVec(nPars) ; meanVec.fill(parMean) ;
  double postScale = 1 / (hyp + sum(square(parVec - parMean)) / 2 ) ;

  Rcout << "gamma mean = " << postShape * postScale  <<
    " gamma Var = " << postShape * postScale * postScale << std::endl;
  double newGamma = rgamma(1, postShape, postScale)(0) ;

  return(1 / newGamma) ;
}

double sampMu(arma::vec parVec, double hypVar, double parVar){
  int nPars = parVec.size() ;
  double postMean = (sum(parVec) / parVar) / ((1 / hypVar) + (nPars / parVar));
  double postVar = 1 / ((1 / hypVar) + (nPars / parVar)) ;

  double newMean = rnorm(1, postMean, sqrt(postVar))(0) ;

  return(newMean) ;
}


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
              NumericMatrix int_mu_,
              NumericMatrix miss_a_,
              NumericMatrix miss_b_,
              NumericMatrix sigma_,
              NumericMatrix tau_int_,
              NumericMatrix tau_fc_,
              NumericMatrix tau_pep_,
              NumericMatrix yVec_,
              Function rProbit){

  int n_prot = matList.size() ;
  //convert to armadillo objects
  arma::mat y_miss(y_miss_.begin(), y_miss_.nrow(), y_miss_.ncol(), false) ;
  arma::mat r_obs(r_obs_.begin(), r_obs_.nrow(), r_obs_.ncol(), false) ;
  arma::mat intercepts(intercepts_.begin(), intercepts_.nrow(),
                       intercepts_.ncol(), false) ;
  arma::mat fcs(fcs_.begin(), fcs_.nrow(), fcs_.ncol(), false) ;
  arma::mat peps(peps_.begin(), peps_.nrow(), peps_.ncol(), false) ;
  arma::mat int_mu(int_mu_.begin(), int_mu_.nrow(), int_mu_.ncol(), false) ;
  arma::mat miss_a(miss_a_.begin(), miss_a_.nrow(), miss_a_.ncol(), false) ;
  arma::mat miss_b(miss_b_.begin(), miss_b_.nrow(), miss_b_.ncol(), false) ;
  arma::mat sigma(sigma_.begin(), sigma_.nrow(), sigma_.ncol(), false) ;
  arma::mat tau_int(tau_int_.begin(), tau_int_.nrow(), tau_int_.ncol(), false) ;
  arma::mat tau_fc(tau_fc_.begin(), tau_fc_.nrow(), tau_fc_.ncol(), false) ;
  arma::mat tau_pep(tau_pep_.begin(), tau_pep_.nrow(), tau_pep_.ncol(), false) ;
  arma::mat yVec(yVec_.begin(), yVec_.nrow(), yVec_.ncol(), false) ;

  //Set missing index before imputing the missing values
  arma::uvec missIndex = find(yVec == 0) ;

  //update missing values and mean parameters
  for(int iter = 0; iter < (tau_pep.n_cols - 1) ; iter++){  // (tau_pep.n_cols - 1)
    arma::vec resids(r_obs.size()) ; resids.zeros() ;
    int residCount = 0 ;
    for(int prot = 0; prot < n_prot; prot++){  //
      NumericMatrix y1 = y_list[prot] ;
      NumericMatrix mat1 = matList[prot] ;
      NumericMatrix point1 = pointers[prot] ;
      arma::vec tempResid = arf::updateBlock(prot, iter, y1, y_miss,
                                 mat1, point1,
                                 intercepts, fcs,
                                 peps, miss_a(iter),
                                 miss_b(iter), sigma(iter),
                                 tau_int(iter), tau_fc(iter),
                                 tau_pep(iter)) ;
      for(int resI = 0; resI < tempResid.size(); resI++){
        resids(resI + residCount) = tempResid(resI) ;
      }
      residCount = residCount + tempResid.size() ;


    } // end protein loop
    Rcout << "Prot loop done "<< std::endl ;

  //update variance components
  tau_int(iter + 1) = arf::sampV(intercepts.col(iter + 1), .001, int_mu(iter)) ;


  tau_fc(iter + 1) = arf::sampV(fcs.col(iter + 1), .001, 0) ;
  tau_pep(iter + 1) = arf::sampV(peps.col(iter + 1), .001, 0) ;

  sigma(iter + 1) = arf::sampV(resids, .001, 0) ;

  Rcout << "VC's updated "<< std::endl ;
  //update mean hyperparameters
  int_mu(iter + 1) = arf::sampMu(intercepts.col(iter + 1),
         10000, tau_int(iter + 1)) ;
  Rcout << "check 1 "<< std::endl ;
  //update missing data parameters
  yVec.elem(missIndex) = y_miss.col(iter + 1) ;

  Rcout << "check 2" << std::endl ;

  NumericVector tempMiss = rProbit(r_obs, yVec) ;
  Rcout << "check 3"<< std::endl ;
  miss_a(iter + 1) = tempMiss(0) ;
  miss_b(iter + 1) = tempMiss(1) ;

  Rcout << "End iteration " << iter << std::endl ;
  } //end iteration loop

  return List::create(Named("fcs") = fcs, Named("intercepts") = intercepts,
                      Named("peps") = peps, Named("int_mu") = int_mu,
                      Named("miss_a") = miss_a, Named("miss_b") = miss_b,
                      Named("sigma") = sigma, Named("tau_int") = tau_int,
                      Named("tau_fc") = tau_fc, Named("tau_pep") = tau_pep,
                      Named("y_miss") = y_miss)  ;

} // End gibbsCpp
