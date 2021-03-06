#Primary functions for running the selection model

#' Fit the Selection Model for Proteomics
#'
#' @param dat Correctly formated and normalized data frame
#' @param ndraws Number of draws of the Gibbs Sampler
#' @param burn Number of draws to discard before summarizing the posterior
#' @param fc_prior Explicitly set a prior for the variance of the fold change
#'   parameter.  By default this is set to zero which uses the usual shared
#'   variance component.  Specifying a weak prior will avoid shrinkage in the
#'   estimates.
#'
#' @export

smp <- function(dat, ndraws = 20000, burn = 1000, melted = FALSE, fc_prior = 0){
  if(melted){
    readyDat <- dat
    pop <- FALSE
    }else{
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
    }

  initList <- prepare(readyDat, ndraws, pop) #function returns, in order:
  #y_list, y_miss, r_obs, matList, pointers,
  # fcs, peps, int_mu, miss_a, miss_b,
  #sigma, tau_int, tau_fc, tau_pep, pop_mu, n_used, estimable, resids

  #reset tau_fc if a prior was selected
  if(fc_prior > 0){
    initList[[13]][1, 1] <- fc_prior
    initList[[14]][1, 1] <- fc_prior
  }

  yVec <- readyDat$lintensity
  yVec[is.na(yVec)] <- 0
  #call the C++ Gibbs Sampler
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
                      as.matrix(yVec), rProbit, sn::rsn,
                      as.matrix(initList[[18]]),
                      fc_prior)

  #extract summary information
  postMeans <- apply(testRes[["fcs"]][ , burn:ndraws], 1, mean)
  postVar <- apply(testRes[["fcs"]][ , burn:ndraws], 1, var)

  protNames <- levels(factor(readyDat$protein))
  conditions <- length(levels(factor(readyDat$condID)))
  nameCol <- rep(protNames, each = (conditions - 1))
  condCol <- rep(2:conditions, length(protNames))

  resTable <- data.frame(Protein = nameCol, Condition = condCol,
                         Mean = postMeans, Var = postVar,
                         N_used = initList[[16]],
                         Estimable = initList[[17]] )


  list(resTable, testRes)


} #end of smp function


