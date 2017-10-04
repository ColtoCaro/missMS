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



compBayes <- function(dat,
                      approx = FALSE,
                      resultsOnly = FALSE,
                      pp=.95,
                      nCores = 1,
                      iter = 2000,
                      nullSet = c(-.2, .2),
                      normalize = TRUE
){

  #Put single dataframe into a list so that we will always work with a list of dataframes
  if(is.data.frame(dat)){dat <- list(dat)}

  #test to make sure each list component is a dataframe
  if (length(dat) > 1){
    testDf <- lapply(dat, is.data.frame)
    if(sum(unlist(testDf)) != length(dat)){
      stop("Error: at least one list component is not a dataframe")
    }
    #make sure that each dataframe has the same reference condition
    refList <- lapply(dat, function(x) paste(x[1, "tag1"],
                                             x[2, "tag1"]))
    if(!do.call(all.equal,refList)){
      stop("Error: Plexes have different reference channels")
    }
  }

  #make sure protein names do not have "_" characters
  dat <- lapply(dat, function(x) within(x, Protein <-
                                          gsub("_", "-", Protein)))

  #determine how many redundancies are being used
  maxRedun <- max(unlist(lapply(dat, function(x) max(x[1, "varCat"]))))


  readyDat <- lapply(1:length(dat), function(x)
    transformDat(dat[[x]], plexNumber = x, normalize = normalize))
  oneDat <- do.call(rbind, readyDat)

  #set data variables
  #Do ptms first since it might change the dimensions of the data
  N_ <- nrow(oneDat)
  sumPtm <- sum(unlist(lapply(dat, function(x) (x[3, 1] == 1))))
  if(sumPtm == 0){
    n_p <- 0
    n_ptm <- 0
    ptm <- rep(0,N_)
    ptmPep <- rep(0,N_)
  }else{
    #find and remove ptm data with no corresponding protein data
    globalProts <- unique(oneDat$bioID[oneDat$ptm == 0])
    ptmProts <- unique(oneDat$bioID[oneDat$ptm > 0])
    orphanProts <- setdiff(ptmProts, globalProts)
    if(length(orphanProts > 0)){
      #remove orphan prots from ptm data
      ptmIndex <- which(oneDat$ptm > 0)
      orphanIndex <- which(oneDat$bioID %in% orphanProts)
      bad <- intersect(ptmIndex, orphanIndex)
      oneDat <- oneDat[-bad, ]
      wText <- paste(length(orphanIndex), "PTM data points were removed because they had no corresponding protein level measurements")
      warning(wText)
    }

    nonPtms <- which(oneDat$ptm == 0)
    ptmName <- levels(factor(oneDat[-nonPtms , ]$ptmID))
    n_p <- length(ptmName)
    n_ptm <- length(unique(oneDat[-nonPtms , ]$ptm))
    ptm <- as.integer(oneDat$ptm)
    ptmPep <- rep(0, nrow(oneDat))
    ptmPep[-nonPtms] <- as.integer(factor(oneDat[-nonPtms , ]$ptmID))

  } # end actions for ptm experiments

  oneDat <- oneDat[order(oneDat$condID, oneDat$bioID, oneDat$ptm,
                         oneDat$ptmID),]

  N_ <- nrow(oneDat)
  n_c <- length(unique(oneDat$condID))
  condKey <- data.frame(number = 1:n_c,
                        name = levels(factor(oneDat$condID)))
  condID <- as.integer(factor(oneDat$condID))

  #Create a tag by varCat variable which determines vc's
  if(maxRedun == 0){
    tagID <- as.integer(factor(oneDat$tag_plex))
  }else{
    redunStr <- paste("R", oneDat$varCat, sep = "")
    if(sumPtm > 0){
      redunStr[which(oneDat$ptm > 0)] <- "Rp"
    }
    oneDat$tag_plex <- paste(oneDat$tag_plex, redunStr, sep="")
    tagID <- as.integer(factor(oneDat$tag_plex))
  }


  n_t <- length(unique(oneDat$tag_plex))


  sumBio <- sum(unlist(lapply(dat, function(x) (x[1, 3] ==1 | x[2,1] == 1)  )))
  if(sumBio == 0){
    n_b <- 0
    bioID <- rep(0,N_)
    condToBio <- matrix(0, nrow = n_c, ncol = 1)
    n_nc <- rep(0, n_c)
    max_nc <- 0
  }else{
    n_b <- length(unique(oneDat$bioID))
    bioID <- as.integer(factor(oneDat$bioID))
    #make a mapping for use in a heierarchical model (not yet implemented)
    bioMap <- oneDat$condID[match(levels(factor(oneDat$bioID)),
                                  oneDat$bioID)]
    conditionNumber <- condKey$number[match(bioMap, condKey$name)]
    bioKey <- data.frame(number = 1:n_b, bioID = levels(factor(oneDat$bioID)),
                         condID = bioMap, conditionNumber)
    bioToCond <- bioKey$conditionNumber

    #make an array giving bio positions for each condition

    n_nc <- unlist(lapply(1:n_c, function(x) sum(bioToCond == x)))
    max_nc <- max(n_nc)
    ncMap <- lapply(1:n_c, function(x) which(bioToCond == x))
    condToBio <- matrix(0, nrow = n_c, ncol = max_nc)
    for (i in 1:n_c){
      condToBio[i, 1:n_nc[i]] <- ncMap[[i]]
    }
  }


  sumCov <- sum(unlist(lapply(dat, function(x) x[1, "Covariate"])))
  useCov <- 1*(sumCov > 0)

  if(useCov){
    covariate <- oneDat$covariate/quantile(oneDat$covariate, probs = pp)
  }else{covariate <- oneDat$covariate}

  lr <- oneDat$lr
  if(sum(lr) == 0){stop("Outcomes are all zero. This might be the
                        consequence of normalizing values already less than
                        one")}

  summaryStr <- paste("Estimating ", max(n_b, n_c), " relative protein abundances, and ", n_p, "protein adjusted ptm changes")
  print(summaryStr)

  #local call for testing
  # model <- stan(file="~/Documents/compMS/exec/allModels.stan",
  #             iter = 2000, cores = 4, control = list(adapt_delta = .8))

  sMod <- compMS:::stanmodels$allModels
  if(approx){
    model <- rstan::vb(sMod, cores = nCores)
  }else{
    model <- rstan::sampling(sMod, cores = nCores, iter = iter)
  }

  #create summary table
  if(n_b > 0){
    targetChain <- rstan::extract(model, pars="avgCond")$avgCond
    postMeans <- colMeans(targetChain)
    postVar <- apply(targetChain, 2, var)
    pvals <- pnorm(nullSet[2], postMeans, sqrt(postVar)) -
      pnorm(nullSet[1], postMeans, sqrt(postVar))
    pValue <- 1 - pchisq((postMeans^2)/postVar, 1)
    justProts <- oneDat[oneDat$ptm == 0 , ]
    n_obs <- unlist(by(justProts$lr, justProts$condID, length))
  }else{
    if(useCov == 0){
      targetChain <- rstan::extract(model, pars="beta")$beta
    }else{
      targetChain <- rstan::extract(model, pars="betaP_c")$betaP_c
    }
    postMeans <- colMeans(targetChain)
    postVar <- apply(targetChain, 2, var)
    pvals <- pnorm(nullSet[2], postMeans, sqrt(postVar)) -
      pnorm(nullSet[1], postMeans, sqrt(postVar))
    pValue <- 1 - pchisq((postMeans^2)/postVar, 1)
    justProts <- oneDat[oneDat$ptm == 0, ]
    n_obs <- unlist(by(justProts$lr, justProts$condID, length))
  }

  resDf <- data.frame(name = levels(factor(oneDat$condID)), mean = postMeans,
                      var = postVar, P_null = pvals,
                      pValue, n_obs = as.vector(n_obs),
                      stringsAsFactors = F)

  if(n_p > 0){
    targetChain <- rstan::extract(model, pars="alpha")$alpha
    postMeans <- colMeans(targetChain)
    postVar <- apply(targetChain, 2, var)
    pvals <- pnorm(nullSet[2], postMeans, sqrt(postVar)) -
      pnorm(nullSet[1], postMeans, sqrt(postVar))
    pValue <- 1 - pchisq((postMeans^2)/postVar, 1)
    justPtms <- oneDat[oneDat$ptm > 0 , ]
    n_obs <- unlist(by(justPtms$lr, justPtms$ptmID, length))
    ptmDf <- data.frame(ptmName, mean = postMeans, var = postVar,
                        P_null = pvals, pValue, n_obs = as.vector(n_obs),
                        stringsAsFactors = F)
  }else{
    ptmDf <- NULL
  }


  RES <- list()
  RES[[1]] <- resDf
  RES[[2]] <- ptmDf
  if(resultsOnly){
    RES[[3]] <- NULL
  }else{
    RES[[3]] <- model
  }

  RES
  } #end of compFit function
