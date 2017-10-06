#Primary functions for running the selection model

#' Fit the Selection Model for Proteomics
#'
#' @param dat Correctly formated and normalized data frame
#' @param nCores Number of cores to use
#' @param ndraws Number of draws of the Gibbs Sampler
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
  #intercepts, fcs, peps, miss_a, miss_b,
  #sigma, tau_int, tau_fc, tau_pep, pop_mu

  #call the C++ Gibbs Sampler
  gibbsCpp(initList[[1]],
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
           as.matrix(initList[[14]]))

   for (i in 1:length(initList)){
   print(typeof(initList[[i]]))
 }

  # RES <- list()
  # RES[[1]] <- resDf
  # RES[[2]] <- ptmDf
  # if(resultsOnly){
  #   RES[[3]] <- NULL
  # }else{
  #   RES[[3]] <- model
  # }
  #
  # RES
} #end of smp function
