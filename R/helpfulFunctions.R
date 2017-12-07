#Helpful internal R functions called from the "main" file


#Function that makes a single groupID from the first 2 rows of df
makeHeader <- function(df, index){

  header <- paste(colnames(df)[index], df[1, index],
                  df[2, index],  sep = "qqqq")
  header
}

#function that takes a dataframe and returns a dataframe with unique ids
transformDat <- function(df){
  print("Transforming data")
  #convert factors to strings
  facIndex <- which(sapply(df, is.factor))
  df[facIndex] <- lapply(df[facIndex], as.character)

  #Zero out unused rows
  if(df[2, 1] == 0){
    df[2, ][] <- 0
  }
  #eliminate rows with zero observations
  missIndex <- apply(df[3:nrow(df), 3:ncol(df)], 1, function(x)
    (sum(is.na(x)) == ncol(df) - 2 ))
  if(sum(missIndex) > 0){
    mI <- which(missIndex == 1) + 2
    df <- df[-mI, ]
  }# end eliminate missing rows

  n_ <- nrow(df)

  jDat <- df[3:(n_), ]
  df[3:n_, ] <- jDat[order(jDat$Protein, jDat$Peptide), ]

  value_index <- grep("Run", colnames(df))

  nMat <- df[3:(n_), value_index]
  lMat <- log2(nMat)

  header <- makeHeader(df[ , value_index])
  colnames(lMat) <- header

  newDf1 <- data.frame(Protein = df[3:n_, ]$Protein,
                       Peptide = df[3:n_, ]$Peptide,                                                    lMat, stringsAsFactors = F)

  melted <- reshape2::melt(newDf1, id.vars = c("Protein", "Peptide"),
                            value.name = "lintensity", variable.name = "header")

  separated <- stringr::str_split_fixed(as.character(melted$header), "qqqq",3)

  #figure out if we are using a bioid from the header or from a column

  finalDat <- data.frame(protein = melted$Protein, peptide = melted$Peptide,
                         condID = separated[, 2],
                         bioID = paste(melted$Protein,
                                       separated[, 3], sep = "_"),
                         lintensity = melted$lintensity, stringsAsFactors = F)

  finalDat <- finalDat[order(finalDat$bioID, finalDat$peptide,
                             finalDat$condID), ]

  finalDat
}#end transformDat()


#function for setting up the data structure and initializing the sampler
prepare <- function(df, ndraws, pop){
  print("Creating parameter matrices")

  #set missing parameters
  missIndex <- which(is.na(df$lintensity))
  missPointer <- rep(0, nrow(df))
  missPointer[missIndex] <- 1:length(missIndex)

  y_miss <- matrix(0, nrow = length(missIndex), ncol = ndraws)
  miss_a <- rep(0, ndraws)
  miss_b <- rep(0, ndraws)
  r_obs <- 1 * (!(is.na(df$lintensity)))
  #make initial guesses
  x0 <- min(df$lintensity[r_obs == 1])
  x100 <- max(df$lintensity[r_obs == 1])
  miss_b[1] <- (qnorm(.8) - qnorm(.01)) / (x100 - x0)
  miss_a[1] <- qnorm(.01) - miss_b[1] * x0
  #make initial guess for missing values.  Don't want complete separation!
  #y_miss[ , 1] <- mean(df$lintensity[r_obs == 1])  #do this later#

  #create list of design matrices
  if(pop == FALSE){
    matList <- by(df, df$protein, FUN = makeX)
  }else{
    matList <- by(df, df$bioID, FUN = makeX)
  }

  n_prot <- length(matList)
  n_parList <- lapply(matList, ncol)
  n_pars <- do.call(sum, n_parList)
  n_cond <- length(unique(df$condID))
  n_pep <- n_pars - ((n_cond - 1) * n_prot)

  #generate parameter matrices
  #intercepts <- matrix(0, nrow = n_prot, ncol = ndraws)
  fcs <- matrix(0, nrow = n_prot * (n_cond - 1), ncol = ndraws)
  n_used <- fcs[ , 1]
  peps <- matrix(0, nrow = n_pep, ncol = ndraws)
  int_mu <- matrix(0, nrow = 1, ncol = ndraws)

  #create a way to map between submatrix calculations and parmater arrays
  #pointer matrix has two columns.  The first represents the type of mean
  #parameter (1=intercept, 2=protein fold-change, 3=peptide effect).
  #The second column determines the position in the respective array.
  pointers <- list()
  int_i <- fc_i <- pep_i <- 0
  for(i in 1:n_prot){
    pointMat <- matrix(0, nrow = n_parList[[i]], ncol = 2)
    for(j in 1:n_parList[[i]]){
      # if(j == 1){
      #   int_i <- int_i + 1
      #   pointMat[j, ] <- c(1, int_i)
      # }
      # if((j > 1) & (j <= n_cond)){
      #   fc_i <- fc_i + 1
      #   pointMat[j, ] <- c(2, fc_i)
      # }
      # if(j > n_cond){
      #   pep_i <- pep_i + 1
      #   pointMat[j, ] <- c(3, pep_i)
      # }
      if(j <= (nrow(pointMat) - n_cond + 1)){
           pep_i <- pep_i + 1
           pointMat[j, ] <- c(3, pep_i)
      }else{
           fc_i <- fc_i + 1
           pointMat[j, ] <- c(2, fc_i)
        }
    }
    pointers[[i]] <- pointMat
  }

  #now put the outcomes into list form
  if(pop == FALSE){
    y_list <- by(data.frame(df$lintensity, missPointer),
                 df$protein, function(x) {
                   x[is.na(x[ , 1]), 1] <- 0
                   return(as.matrix(x))})
  }else{
    y_list <- by(data.frame(df$lintensity, missPointer),
                 df$bioID, function(x) {
                   x[is.na(x[ , 1]), 1] <- 0
                   return(as.matrix(x))})
  }

  #closure function for computing intial parameter estimates
  ols_init <- function(y_, X_, pointers){
    vec <- rep(0, nrow(y_))
    isMiss <- which(y_[ , 1] == 0)
    vec[isMiss] <- NA
    vec[y_[ , 1] != 0] <- y_[which(y_[ , 1] != 0), 1]

    beta <- coef(lm(vec ~ 0 + X_))
    obsUsed <- apply(t(t(solve(t(X_) %*% X_) %*% t(X_)) * vec), 1, sumObserved)

    #intercepts[pointers[1, 2], 1] <<- beta[1]
    index <- which(pointers[ , 1] == 2)
    fcs[pointers[index , 2], 1] <<- beta[index]
    n_used[pointers[index , 2]] <<- obsUsed[index]
    index <- which(pointers[ , 1] == 3)
    peps[pointers[index , 2], 1] <<- beta[index]

    #enter initial missing values
    beta0 <- beta
    beta0[which(is.na(beta0))] <- 0
    xi0 <- X_ %*% beta0
    y_miss[y_[isMiss , 2], 1] <<- xi0[isMiss]

    piter <- get("x", parent.frame())
    if(piter == 500){print("500 proteins initialized")}
    if(piter == 1000){print("1000 proteins initialized")}
    if(piter == 2500){print("2500 proteins initialized")}
    if(piter == 5000){print("5000 proteins initialized")}

    NULL
  } #end ols_init()

  #generate initial mean parameters
  print("Computing inital parameter values")
  invisible(lapply(1:n_prot, function(x)
    ols_init(y_list[[x]], matList[[x]], pointers[[x]])))
  print("done with OLS")

  #set intercept hyper mean
  int_mu[1] <- mean(peps[ , 1])
  #Take care of inestimable proteins
  estimable <- !is.na(fcs[ , 1])
  fcs[which(estimable == FALSE), 1] <- 0

  #Now do the variance components
  #generate parameter matrices
  sigma <- matrix(0, nrow = 1, ncol = ndraws)
  tau_int <- matrix(0, nrow = 1, ncol = ndraws)
  tau_fc <- matrix(0, nrow = 1, ncol = ndraws)
  tau_pep <- matrix(0, nrow = 1, ncol = ndraws)

  sigma[1 , 1] <- .001
  tau_int[1 , 1] <- 2
  tau_fc[1 , 1] <- 1
  tau_pep[1 , 1] <- 1

  #create higher level for population studies
  if(pop == FALSE){pop_mu <- NULL}else{
    pop_mu <- NULL # Work on this later
  }
  #create matrix for residuals
  resids <- matrix(0, nrow = nrow(df), ncol = ndraws)

  list(y_list, y_miss, r_obs, matList, pointers,
       fcs, peps, int_mu, miss_a, miss_b,
       sigma, tau_int, tau_fc, tau_pep, pop_mu, n_used, estimable, resids)
} #end prepare()



#function for creating design matrices
makeX <- function(df){
  multiPep <- (length(unique(df$peptide)) > 1)
    if(multiPep){
      mat <- model.matrix(~ 0 + factor(peptide) + factor(condID), df)
    }else{
      mat <- model.matrix(~ factor(condID), df)
    }
} #end makeX()

#function for adding up the number of observed non-zero values in a row
sumObserved <- function(vec){
  vec <- round(vec, 10)
  vec[vec == 0] <- NA
  sum(!is.na(vec))
}

# ccX <- function(df){
#   missIndex <- which(is.na(df$lintensity))
#   df <- df[-missIndex, ]
#   multiCond <- (length(unique(df$condID)) > 1)
#   if(multiCond){
#     mat <- model.matrix(~ 0 + factor(peptide) + factor(condID), df)
#   }else{
#     mat <- model.matrix(~ 0 + factor(peptide), df)
#   }
# } #end ccX()

#function designed to fit a probit regression and extract relevant info
rProbit <- function(boolVec, covar){
  maxMiss <- max(covar[boolVec == 0])
  minObs <- min(covar[boolVec == 1])
  minMiss <- min(covar[boolVec == 0])
 #print(paste("minMiss = ", minMiss, " maxMiss = ", maxMiss,
 #             "minObs = ", minObs ))
  if(maxMiss < minObs){stop("Complete separation has occured!")}
  mod <- glm(boolVec ~ covar,
                              family = binomial(link = "probit"))
  newMiss <- mvtnorm::rmvnorm(1, mod$coefficients, vcov(mod))
  c(newMiss[1], newMiss[2])
}

#Maybe need the below, not sure yet
#function for extracting the condition number from labels
getCond <- function(strVec, ptm = FALSE){
  sPosition <- regexpr("_", strVec)
  if(ptm){
    subbed <- sub("_", "*", strVec)
    ePosition <- regexpr("_", subbed)
    condNumber <- as.integer(substring(strVec, sPosition +1, ePosition -1))
  }else{
    condNumber <- as.integer(substring(strVec, sPosition +1))
  }

  condNumber
}

#function to reverse strings
strReverse <- function(x){
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
}




