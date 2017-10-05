#Primary functions for running the selection model

#' Generate a random number for Gaussian mean parameter
#'
#' @param n_ The number of random variates to create
#' @param xi esn parameter
#' @param omega esn parameter
#' @param alpha esn parameter
#' @param tau esn parameter
#'
#' @export

smp <- function(dat, nCores = 1, ndraws = 2000){

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
  #intercepts, fcs, peps, miss_a, miss_a,
  #sigma, tau_int, tau_fc, tau_pep, pop_mu

  #call the C++ Gibbs Sampler



  RES <- list()
  RES[[1]] <- resDf
  RES[[2]] <- ptmDf
  if(resultsOnly){
    RES[[3]] <- NULL
  }else{
    RES[[3]] <- model
  }

  RES
} #end of smp function
