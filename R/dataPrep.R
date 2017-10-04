#File for re-exploring the dilution data


agc <- (dat$Ion.Injection.Time < max(dat$Ion.Injection.Time))
dat <- data.frame(dat, agc)
dat2 <- dat[-grep("##", dat$Reference),
            c("z", "SrchName", "Reference", "Trimmed.Peptide", "Gene.Symbol",
              "Precursor.Max.Intensity", "agc" )]

df <- data.frame(run = dat2$SrchName, prot = dat2$Reference, pep = dat2$Trimmed.Peptide, gene = dat2$Gene.Symbol, intensity = dat2$Precursor.Max.Intensity,
                 agc = dat2$agc)

ordered <- df[order(df$prot, df$pep, df$run), ]

uid <- paste(ordered$prot, ordered$pep, ordered$run)

onlyOne <- function(df){
  index <- which.max(df$intensity)
  newdf <- df[index, ]
  newdf
}

oneList <- by(ordered, uid, onlyOne)
reduced <- do.call(rbind, oneList)

ironI <- dcast(reduced, prot + pep + gene ~ run, value.var = "intensity")
ironA <- dcast(reduced, prot + pep + gene ~ run, value.var = "agc")
#take a look at the HEK cell summaries
noAGC <- (apply(ironA[ , 4:8], 1, sum, na.rm=T) == 0)
noMiss <- !(is.na(apply(ironI[ , 4:8], 1, sum)))
perfect <- noAGC*noMiss

HEKI <- ironI[which(perfect == T) , 1:8]
cmeans <- apply(HEKI[, 4:8], 2, mean, na.rm = T)
facs <- max(cmeans)/(cmeans * c(100,16,4,1,1))
normHEK <- t(t(ironI[ , 4:8]) * facs)
HEK <- data.frame(ironI[,1:3], normHEK)
resid <- c(mean((HEK[,4]/HEK[,7]),na.rm=T) - 1/100,
           mean((HEK[,5]/HEK[,8]),na.rm=T) - 1/16,
           mean((HEK[,6]/HEK[,7]),na.rm=T) - 1/4,
           mean((HEK[,8]/HEK[,7]),na.rm=T) - 1)
title <- paste("residuals \n", paste(signif(resid, 2), collapse = " "))
plot(resid, ylab = "", main = title)

sds <- c(sd((HEK[,4]/HEK[,7]),na.rm=T),
         sd((HEK[,5]/HEK[,8]),na.rm=T),
         sd((HEK[,6]/HEK[,7]),na.rm=T),
         sd((HEK[,8]/HEK[,7]),na.rm=T))
title <- paste("SDs \n", paste(signif(sds, 2), collapse = " "))
plot(sds, ylab = "",main = title)

#take a look at the SH-SY5Y cell summaries
noAGC <- (apply(ironA[ , 10:14], 1, sum, na.rm=T) == 0)
noMiss <- !(is.na(apply(ironI[ , 10:14], 1, sum)))
perfect <- noAGC*noMiss

matI <- ironI[which(perfect == T) , 10:14]
cmeans <- apply(matI, 2, mean, na.rm = T)
facs <- max(cmeans)/(cmeans * c(100,100,16,4,1))
normdat <- t(t(ironI[ , 10:14]) )
resid <- c(median((normdat[,5]/normdat[,1]),na.rm=T) - 100,
           median((normdat[,5]/normdat[,2]),na.rm=T) - 100,
           median((normdat[,5]/normdat[,3]),na.rm=T) - 16,
           median((normdat[,5]/normdat[,4]),na.rm=T) - 4,
           median((normdat[,2]/normdat[,1]),na.rm=T) - 1)
title <- paste("residuals \n", paste(signif(resid, 2), collapse = " "))
plot(resid, ylab = "", main = title, type ="l", col = "red", ylim = c(-70,30))
lines(mqResid, type = "l", col = "blue")

heasds <- c(sd((HEK[,4]/HEK[,7]),na.rm=T),
            sd((HEK[,5]/HEK[,8]),na.rm=T),
            sd((HEK[,6]/HEK[,7]),na.rm=T),
            sd((HEK[,8]/HEK[,7]),na.rm=T))
title <- paste("SDs \n", paste(signif(sds, 2), collapse = " "))
plot(sds, ylab = "",main = title)


mqDat <- read.table("cleanPeps_AOAS.txt", sep = "\t", header = T)
mqMat <- mqDat[, 4:18]
mqMat[mqMat == 0] <- NA
mqDat <- data.frame(mqDat[, 1:3], mqMat)
mean(mqMat$Intensity.70/mqMat$Intensity.66, na.rm = T)
mean(mqMat$Intensity.70/mqMat$Intensity.67, na.rm = T)
mean(mqMat$Intensity.66/mqMat$Intensity.67, na.rm = T)

mqResid <- c(mean(mqMat$Intensity.56/mqMat$Intensity.52, na.rm = T) - 100,
             mean(mqMat$Intensity.56/mqMat$Intensity.53, na.rm = T) - 16,
             mean(mqMat$Intensity.56/mqMat$Intensity.54, na.rm = T) - 4,
             mean(mqMat$Intensity.55/mqMat$Intensity.56, na.rm = T) - 1)

HEKdat <- data.frame(Protein = paste("HEK", mqDat$Leading.razor.protein),
                     Peptide = mqDat$Sequence, Run1 = mqDat$Intensity.52,
                     Run2 = mqDat$Intensity.53, Run3 = mqDat$Intensity.54,
                     Run4 = mqDat$Intensity.55, Run5 = mqDat$Intensity.56)
missRows <- apply(HEKdat[3:7], 1, function(vec)
  (sum(is.na(vec)) == length(vec)))
HEKdat <- HEKdat[-which(missRows == TRUE), ]
nHEK <- length(unique(HEKdat$Protein))

HLdat <- data.frame(Protein = paste("HeLa", mqDat$Leading.razor.protein),
                    Peptide = mqDat$Sequence, Run1 = mqDat$Intensity.59,
                    Run2 = mqDat$Intensity.74, Run3 = mqDat$Intensity.75,
                    Run4 = mqDat$Intensity.76, Run5 = mqDat$Intensity.77)
missRows <- apply(HLdat[3:7], 1, function(vec)
  (sum(is.na(vec)) == length(vec)))
HLdat <- HLdat[-which(missRows == TRUE), ]
nHL <- length(unique(HLdat$Protein))


SHdat <- data.frame(Protein = paste("SH-SY5Y", mqDat$Leading.razor.protein),
                    Peptide = mqDat$Sequence, Run1 = mqDat$Intensity.66,
                    Run2 = mqDat$Intensity.67, Run3 = mqDat$Intensity.68,
                    Run4 = mqDat$Intensity.69, Run5 = mqDat$Intensity.70)
missRows <- apply(SHdat[3:7], 1, function(vec)
  (sum(is.na(vec)) == length(vec)))
SHdat <- SHdat[-which(missRows == TRUE), ]
nSH <- length(unique(SHdat$Protein))

allDat <- rbind(HEKdat, HLdat, SHdat)
allDat <- allDat[order(allDat$Protein, allDat$Peptide),]
protKey <- data.frame(Protein = levels(factor(allDat$Protein)),
                      Run1 = c(rep(1/100, nHEK),
                               rep(1/100, nHL),
                               rep(1/100, nSH)),
                      Run2 = c(rep(1/16, nHEK),
                               rep(1/16, nHL),
                               rep(1/100, nSH)),
                      Run3 = c(rep(1/4, nHEK),
                               rep(1/4, nHL),
                               rep(1/16, nSH)),
                      Run4 = c(rep(1, nHEK),
                               rep(1/4, nHL),
                               rep(1/4, nSH)),
                      Run5 = c(rep(1, nHEK),
                               rep(1, nHL),
                               rep(1, nSH))
) #end make protKey

set.seed(777)
newOrder <- matrix(NA, nrow = nrow(protKey), ncol = 5)
for (i in 1:nrow(protKey)){
  newOrder[i, ] <- sample(1:5)
}
permuteKey <- t(sapply(1:nrow(protKey), function(x) protKey[x, 1+ newOrder[x,]]))
colnames(permuteKey) <- colnames(protKey[2:6])
permuteKey <- data.frame(Protein = protKey$Protein, permuteKey)
permuteKey <- lapply(permuteKey, unlist)
permuteKey <- as.data.frame(permuteKey)

useKey <- function(subdat){
  lookupID <- grep(subdat[1,1], permuteKey$Protein)
  perm <- newOrder[lookupID, ]
  newDat <- subdat[, c(1:2, 2 + perm)]
  colnames(newDat) <- colnames(subdat)
  newDat
}
#refactor allDat
allDat$Protein <- factor(allDat$Protein)
permuteDat <- by(allDat, allDat$Protein, useKey)
permuteDf <- do.call(rbind, permuteDat)

#add a header
header <- matrix(0, nrow = 2, ncol = ncol(permuteDf))
colnames(header) <- colnames(permuteDf)
sampleDat <- rbind(header, permuteDf)

devtools::use_data(sampleDat)
write.csv(permuteDf, "permuted.csv")
write.csv(permuteKey, "ratioKey.csv")

