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
  y_miss[ , 1] <- mean(df$lintensity[r_obs == 1])

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
    vec[y_[ , 1] == 0] <- y_miss[y_[ , 2], 1]
    vec[y_[ , 1] != 0] <- y_[which(y_[ , 1] != 0), 1]

    beta <- solve(t(X_) %*% X_) %*% t(X_) %*% vec

    #intercepts[pointers[1, 2], 1] <<- beta[1]
    index <- which(pointers[ , 1] == 2)
    fcs[pointers[index , 2], 1] <<- beta[index]
    index <- which(pointers[ , 1] == 3)
    peps[pointers[index , 2], 1] <<- beta[index]

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

  list(y_list, y_miss, r_obs, matList, pointers,
       fcs, peps, int_mu, miss_a, miss_b,
       sigma, tau_int, tau_fc, tau_pep, pop_mu)
} #end prepare()



#function for creating design matrices
makeX <- function(df){
  multiPep <- (length(unique(df$peptide)) > 1)
    if(multiPep){
      mat <- model.matrix(~ 0 + factor(peptide) + factor(condID), df)
    }else{
      mat <- model.matrix(~ 0 + factor(condID), df)
    }
} #end makeX()

#function designed to fit a probit regression and extract relevant info
rProbit <- function(boolVec, covar){
  maxMiss <- max(covar[boolVec == 0])
  minObs <- min(covar[boolVec == 1])
  minMiss <- min(covar[boolVec == 0])
 print(paste("minMiss = ", minMiss, " maxMiss = ", maxMiss,
              "minObs = ", minObs ))
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





