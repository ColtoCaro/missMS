#Primary functions for running the selection model

#' Fit the Selection Model for Proteomics
#'
#' @param dat Correctly formated and normalized data frame
#' @param nCores Number of cores to use
#' @param ndraws Number of draws of the Gibbs Sampler
#'
#' @export

smp <- function(dat, nCores = 1, ndraws = 2000, pop = FALSE){

  if(dat[2, 1] == 1){
    pop <- TRUE
    stop("Sorry, the population level model is currently in development")
  }else{
    pop <- FALSE
  }
  #make sure protein names do not have "_" characters
  dat <- within(dat, Protein <- gsub("_", "-", Protein))

  #transform data and initialize the Gibbs sampler
  readyDat <- transformDat(dat)
  initList <- prepare(readyDat, ndraws, pop) #function returns, in order:
  #y_list, y_miss, r_obs, matList, pointers,
  # fcs, peps, int_mu, miss_a, miss_b,
  #sigma, tau_int, tau_fc, tau_pep, pop_mu, n_used, estimable, resids

  yVec <- readyDat$lintensity
  yVec[is.na(yVec)] <- 0
  #call the C++ Gibbs Sampler
  set.seed(777)
  testRes <- gibbsCpp(initList[[1]],
           as.matrix(initList[[2]]),
           as.matrix(initList[[3]]),
           initList[[4]],
           initList[[5]],
           as.matrix(initList[[6]]),
           as.matrix(initList[[7]]),
           as.matrix(initList[[8]]),
           as.matrix(initList[[9]]),
           as.matrix(initList[[10]]),
           as.matrix(initList[[11]]),
           as.matrix(initList[[12]]),
           as.matrix(initList[[13]]),
           as.matrix(initList[[14]]),
           as.matrix(yVec), rProbit, rsn,
           as.matrix(initList[[18]]))

#5464 might be trouble
  levels(factor(readyDat$protein))[5464]
} #end of smp function


