lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p[1])
}


# initialize study variables
heritability = 0.3
heritabilityTopProp = 0.7
heritabilityBottomProp = 0
heritabilityRareProp = 1 - heritabilityTopProp - heritabilityBottomProp
numVariantsTop = 200
numVariantsBottom = 0
numVariantsRare = 2
baselinePrevalence = 0.1

numIndividuals = 500000

# define our functions to go back and forth between p and eta

# phi: eta --> p
phi<-function(eta){
  1 / (1 + exp(-eta))
  
}

# psi: p --> eta
psi<-function(p){
  log(p / (1 - p))
  
}

# sign flip 
sign.flip<-function(n,p){
  coins <- rbinom(n,1,p)
  coins[coins == 0] = -1
  return(coins)}


baselineScore = psi(baselinePrevalence)

assignVariantFrequencies<-function(){
  # variant frequencies top 
  variantFrequenciestop<<-rbeta(numVariantsTop, 1.5, 8.5)
  # variant frequencies bottom 
  variantFrequenciesbottom<<-rbeta(numVariantsBottom, 1.5, 8.5)
  # variant frequencies rare
  #variantFrequenciesrare<<-c(.005,.005)
  myVariances = heritabilityRareProp*heritability / numVariantsRare
  variantFrequenciesrare<<-rep((1 - sqrt(1-2*myVariances/4)) / 2, 2)
  
}


assignVariantBetaVariances<-function(){
  # variances top
  variantBetaVariancesTop<<-rep(heritabilityTopProp*heritability / numVariantsTop, numVariantsTop)
  # variances bottom
  variantBetaVariancesBottom<<-rep(heritabilityBottomProp*heritability / numVariantsBottom, numVariantsBottom)
  # variances rare
  variantBetaVariancesrare<<-rep(heritabilityRareProp*heritability / numVariantsRare, numVariantsRare)
  
}

assignVariantBetas<-function(){
  # Betas top
  variantBetasTop<<-sign.flip(numVariantsTop,.5)*sqrt(variantBetaVariancesTop / (2 * variantFrequenciestop * (1 - variantFrequenciestop)))
  #Betas bottom
  variantBetasBottom<<-sign.flip(numVariantsBottom,.5)*sqrt(variantBetaVariancesBottom / (2 * variantFrequenciesbottom * (1 - variantFrequenciesbottom)))
  #Betas rare
  variantBetasRare<<- -sqrt(variantBetaVariancesrare / (2 * variantFrequenciesrare * (1 - variantFrequenciesrare)))
  
}

# how to generate a genotype
generateGenotype<-function(){
  rbinom(numVariantsTop + numVariantsBottom + numVariantsRare, 2, c(variantFrequenciestop, variantFrequenciesbottom, variantFrequenciesrare))
  
}

# how to take PRS of genotype
PRS<-function(genotype){
  sum(genotype * c(variantBetasTop, variantBetasBottom, variantBetasRare))
  
}

# how to get total risk score
riskScore<-function(genotype){
  baselineScore + PRS(genotype)
  
}

createIndividuals<-function(){
  # generate genotypes
  individualGenotypes<<-replicate(numIndividuals, generateGenotype())
  
  # generate our riskScores
  individualScores<<-apply(individualGenotypes, 2, riskScore)
  individualPRSs<<-individualScores - baselineScore
  
  # generate phenotypes
  individualPhenotypes<<-rbinom(numIndividuals, 1, phi(individualScores))
  
}
createIndividualsBX<-function(){
  # generate genotypes
  individualGenotypes<<-replicate(numIndividuals, generateGenotype())
  
  # lets add some correlation for 2 rare variants
  # MAKE THIS MORE GENERALIZED
  for(i in 1:numIndividuals){
    if(individualGenotypes[201,i] > 0 & individualGenotypes[202,i] < 2){
      individualGenotypes[202, i]<<-individualGenotypes[202, i] + 1
      
    }
    
  }
  
  # generate our riskScores
  individualScores<<-apply(individualGenotypes, 2, riskScore)
  individualPRSs<<-individualScores - baselineScore
  
  # generate phenotypes
  individualPhenotypes<<-rbinom(numIndividuals, 1, phi(individualScores))
  
}


# how to estimate PRS of genotype
estimatePRS<-function(genotype){
  sum(genotype[1:numVariantsTop] * betaEstimates)
  
}

createPRSEstimates<-function(){
  # estimate our Beta values
  betaEstimates<<- rnorm(numVariantsTop, variantBetasTop, 0.01)
  # estimate everyone's PRS
  individualPRSEstimates<<- apply(individualGenotypes, 2, estimatePRS)
  
}

# *************** Analysis ***************

# in analysis, we can use \beta_{m+1} through \beta_n as seperate attempts 
# since we only estimate one at a time


# *************** Method A ***************

# First attempt: measure beta by finding disease prevalence in those with and without
# no information needed other than phenotypes and G_{m+1}

# initialize our estimate vector
#unknownEstimatesA = rep(NaN, numVariantsRare)

# how to estimate a beta
estimateA<-function(locus){
  
  # select indices for each number of alleles
  gZero = which(individualGenotypes[locus,] == 0)
  gOne = which(individualGenotypes[locus,] == 1)
  gTwo = which(individualGenotypes[locus,] == 2)
  
  # find disease prevalence in these groups
  zeroAllelePrevalence = sum(individualPhenotypes[gZero]) / length(gZero)
  oneAllelePrevalence = sum(individualPhenotypes[gOne]) / length(gOne)
  twoAllelePrevalence = sum(individualPhenotypes[gTwo]) / length(gTwo)
  
  # use prevalence to estimate beta
  # just use prevalence for 0 and 1 alleles. too few have 2 normally.
  psi(oneAllelePrevalence) - psi(zeroAllelePrevalence)
  
}


# *************** Method B ***************

# select individuals without the disease and compare estimated PRS means
# in those with 1 allele and with 0

# initialize our estimate vector
#unknownEstimatesB = rep(NaN, numVariantsRare*2)


# how to estimate a beta
estimateB<-function(locus){
  
  # select out healthy inds
  healthy<<-which(individualPhenotypes == 0)
  
  # difference should be beta
 fit <- lm(individualPRSEstimates[healthy] ~ individualGenotypes[locus,healthy])
  return(list(-summary(fit)$coefficients[2, "Estimate"], variantBetasRare[locus - numVariantsBottom - numVariantsTop], summary(fit)$coefficients[2, "Pr(>|t|)"]))
  #print(PRSZero - PRSOne)
  
}

studynumber = 0

simulateStudyB<-function(){
  studynumber<<-studynumber + 1
  print(cat("Beginning Study ", studynumber))
  print("Assigning betas...")
  assignVariantFrequencies()
  assignVariantBetaVariances()
  assignVariantBetas()
  
  print("Creating individuals...")
  createIndividuals()
  createPRSEstimates()
  
  estimateB(numVariantsBottom + numVariantsTop + 1)
  
  
}


# # *************** Method C ***************
# 
# # select individuals with the disease and compare estimated PRS means
# # in those with 1 allele and with 0
# # same as B except choosing those with
# 
# # initialize our estimate vector
# unknownEstimatesC = rep(NaN, numVariants - numKnownVariants)
# 
# 
# # how to estimate a beta
# estimateC<-function(locus){
#   
#   # select out healthy inds
#   healthy = which(individualPhenotypes == 1)
#   
#   # select those with 1 and 0 alleles, as those with 2 are too few
#   gZero = which(individualGenotypes[locus, healthy] == 0)
#   gOne = which(individualGenotypes[locus, healthy] == 1)
#   
#   # find PRS estimate means on populations
#   PRSZero = mean(individualPRSEstimates[gZero])
#   PRSOne = mean(individualPRSEstimates[gOne])
#   
#   # difference should be beta
#   PRSOne - PRSZero
#   
# }
# 
# # estimate for each unknown beta
# for(i in (numKnownVariants + 1):numVariants){
#   unknownEstimatesC[i - numKnownVariants] = estimateC(i)
#   
# }

# *************** Method D ***************

# Take the tails and treat them as cases and controls

estimateD<-function(locus){
  
  # select out healthy inds
  healthy<<-which(individualPhenotypes == 0)
  
  # set fraction for tails
  tailsize = 0.10
  
  cutoffs = quantile(individualPRSEstimates[healthy], c(tailsize, 1-tailsize))
  topPRSPeople = healthy[individualPRSEstimates[healthy] > cutoffs[2]]
  bottomPRSPeople = healthy[individualPRSEstimates[healthy] < cutoffs[1]]
  allPeople = c(bottomPRSPeople, topPRSPeople)
  
  # select indices for each number of alleles
  gZero = allPeople[individualGenotypes[locus, allPeople] == 0]
  gOne = allPeople[individualGenotypes[locus, allPeople] == 1]
  gTwo = allPeople[individualGenotypes[locus, allPeople] == 2]
  
  # if we don't have any with gTwo, then we can't really find our beta
  if(length(gTwo) == 0) return("NH2")
  
  # find disease prevalence in these groups
  zeroAllelePrevalence = length(which(gZero %in% topPRSPeople)) / length(gZero)
  oneAllelePrevalence = length(which(gOne %in% topPRSPeople)) / length(gOne)
  twoAllelePrevalence = length(which(gTwo %in% topPRSPeople)) / length(gTwo)
  
  # can't use if twoAllelePrevalence == 0 or 1
  if(twoAllelePrevalence == 0 || twoAllelePrevalence == 1) return("IPI")
  
  # use prevalence to estimate beta (if we can)
  fit<-lm(c(psi(zeroAllelePrevalence), psi(oneAllelePrevalence), psi(twoAllelePrevalence)) ~ c(0, 1, 2))
  return(list(-summary(fit)$coefficients[2, "Estimate"], variantBetasRare[locus - numVariantsBottom - numVariantsTop], summary(fit)$coefficients[2, "Pr(>|t|)"]))
  
}

simulateStudyD<-function(){
  studynumber<<-studynumber + 1
  print(cat("Beginning Study ", studynumber))
  print("Assigning betas...")
  assignVariantFrequencies()
  assignVariantBetaVariances()
  assignVariantBetas()
  
  print("Creating individuals...")
  createIndividuals()
  createPRSEstimates()
  
  estimateD(numVariantsBottom + numVariantsTop + 1)
  
  
}

# *************** Method E ***************

# Same as B, but use inner portions

estimateE<-function(locus){
  
  # select out healthy inds
  healthy<<-which(individualPhenotypes == 0)
  
  # get genotype sets
  gZero = healthy[individualGenotypes[locus, healthy] == 0]
  gOne = healthy[individualGenotypes[locus, healthy] == 1]
  gTwo = healthy[individualGenotypes[locus, healthy] == 2]
  
  # find quantiles for each genotype
  innerPortion = 0.5
  quantilesZero = quantile(individualPRSEstimates[gZero], c((1 - innerPortion) / 2, (1 + innerPortion) / 2))
  quantilesOne = quantile(individualPRSEstimates[gOne], c((1 - innerPortion) / 2, (1 + innerPortion) / 2))
  quantilesTwo = quantile(individualPRSEstimates[gTwo], c((1 - innerPortion) / 2, (1 + innerPortion) / 2))
  
  # find inner portions of PRS spreads for each genotype
  gZeroInner = gZero[individualPRSEstimates[gZero] > quantilesZero[1]]
  gZeroInner = gZeroInner[individualPRSEstimates[gZeroInner] < quantilesZero[2]]
  gOneInner = gOne[individualPRSEstimates[gOne] > quantilesOne[1]]
  gOneInner = gOneInner[individualPRSEstimates[gOneInner] < quantilesOne[2]]
  gTwoInner = gTwo[individualPRSEstimates[gTwo] > quantilesTwo[1]]
  gTwoInner = gTwoInner[individualPRSEstimates[gTwoInner] < quantilesTwo[2]]
  
  # assign weights
  weightZero = 1
  weightOne = 1 #sqrt(length(gZeroInner) / length(gOneInner))
  weightTwo = 1 #sqrt(length(gZeroInner) / length(gTwoInner))
  
  # fitting
  relevantPRSs = c(individualPRSEstimates[gZeroInner], individualPRSEstimates[gOneInner], individualPRSEstimates[gTwoInner])
  relevantGenotypes = c(rep(0, length(gZeroInner)), rep(1, length(gOneInner)), rep(2, length(gTwoInner)))
  relevantWeights = c(rep(weightZero, length(gZeroInner)), rep(weightOne, length(gOneInner)), rep(weightTwo, length(gTwoInner)))
  plot(relevantGenotypes, relevantPRSs)
  fit <- lm(relevantPRSs ~ relevantGenotypes, weights=relevantWeights)
  return(list(-summary(fit)$coefficients[2, "Estimate"], variantBetasRare[locus - numVariantsBottom - numVariantsTop], summary(fit)$coefficients[2, "Pr(>|t|)"]))
  
  
}

simulateStudyE<-function(){
  studynumber<<-studynumber + 1
  print(cat("Beginning Study ", studynumber))
  print("Assigning betas...")
  assignVariantFrequencies()
  assignVariantBetaVariances()
  assignVariantBetas()
  
  print("Creating individuals...")
  createIndividuals()
  createPRSEstimates()
  
  print("Estimating...")
  estimateE(numVariantsBottom + numVariantsTop + 1)
  
  
}

simulateStudyBE<-function(){
  studynumber<<-studynumber + 1
  print(cat("Beginning Study ", studynumber))
  print("Assigning betas...")
  assignVariantFrequencies()
  assignVariantBetaVariances()
  assignVariantBetas()
  
  print("Creating individuals...")
  createIndividuals()
  createPRSEstimates()
  
  print("Estimating...")
  c(unlist(estimateB(numVariantsBottom + numVariantsTop + 1)[1]), unlist(estimateE(numVariantsBottom + numVariantsTop + 1)[1]))
  
  
}

simulateStudyBD<-function(){
  studynumber<<-studynumber + 1
  print(cat("Beginning Study ", studynumber))
  print("Assigning betas...")
  assignVariantFrequencies()
  assignVariantBetaVariances()
  assignVariantBetas()
  
  print("Creating individuals...")
  createIndividuals()
  createPRSEstimates()
  
  print("Estimating...")
  c(unlist(estimateB(numVariantsBottom + numVariantsTop + 1)[1]), unlist(estimateD(numVariantsBottom + numVariantsTop + 1)[1]))
  
  
}

# *************** Method F ***************

# Same as D but with weights

estimateF<-function(locus){
  
  # select out healthy inds
  healthy<<-which(individualPhenotypes == 0)
  
  # set fraction for tails
  tailsize = 0.10
  
  cutoffs = quantile(individualPRSEstimates[healthy], c(tailsize, 1-tailsize))
  topPRSPeople = healthy[individualPRSEstimates[healthy] > cutoffs[2]]
  bottomPRSPeople = healthy[individualPRSEstimates[healthy] < cutoffs[1]]
  allPeople = c(bottomPRSPeople, topPRSPeople)
  
  # select indices for each number of alleles
  gZero = allPeople[individualGenotypes[locus, allPeople] == 0]
  gOne = allPeople[individualGenotypes[locus, allPeople] == 1]
  gTwo = allPeople[individualGenotypes[locus, allPeople] == 2]
  
  # if we don't have any with gTwo, then we can't really find our beta
  if(length(gTwo) == 0) return("NH2")
  
  # find disease prevalence in these groups
  zeroAllelePrevalence = length(which(gZero %in% topPRSPeople)) / length(gZero)
  oneAllelePrevalence = length(which(gOne %in% topPRSPeople)) / length(gOne)
  twoAllelePrevalence = length(which(gTwo %in% topPRSPeople)) / length(gTwo)
  
  # can't use if twoAllelePrevalence == 0 or 1
  if(twoAllelePrevalence == 0 || twoAllelePrevalence == 1) return("IPI")
  
  # attach weights for confidence
  myWeights = c(length(gZero), length(gOne), length(gTwo))
  
  # use prevalence to estimate beta (if we can)
  fit<-lm(c(psi(zeroAllelePrevalence), psi(oneAllelePrevalence), psi(twoAllelePrevalence)) ~ c(0, 1, 2), weights=myWeights)
  return(list(-summary(fit)$coefficients[2, "Estimate"], variantBetasRare[locus - numVariantsBottom - numVariantsTop], summary(fit)$coefficients[2, "Pr(>|t|)"]))
  
}

simulateStudyDF<-function(){
  studynumber<<-studynumber + 1
  print(cat("Beginning Study ", studynumber))
  print("Assigning betas...")
  assignVariantFrequencies()
  assignVariantBetaVariances()
  assignVariantBetas()
  
  print("Creating individuals...")
  createIndividuals()
  createPRSEstimates()
  
  print("Estimating...")
  c(unlist(estimateD(numVariantsBottom + numVariantsTop + 1)[1]), unlist(estimateF(numVariantsBottom + numVariantsTop + 1)[1]))
  
  
}

simulateStudyF<-function(){
  studynumber<<-studynumber + 1
  print(cat("Beginning Study ", studynumber))
  print("Assigning betas...")
  assignVariantFrequencies()
  assignVariantBetaVariances()
  assignVariantBetas()
  
  print("Creating individuals...")
  createIndividuals()
  createPRSEstimates()
  
  print("Estimating...")
  estimator = estimateF(numVariantsBottom + numVariantsTop + 1)
  c(unlist(estimator[1]), unlist(estimator[3]))
  
  
}


simulateStudyAB<-function(){
  studynumber<<-studynumber + 1
  print(cat("Beginning Study ", studynumber))
  print("Assigning betas...")
  assignVariantFrequencies()
  assignVariantBetaVariances()
  assignVariantBetas()
  
  print("Creating individuals...")
  createIndividuals()
  createPRSEstimates()
  
  print("Estimating...")
  c(unlist(estimateA(numVariantsBottom + numVariantsTop + 1)), unlist(estimateB(numVariantsBottom + numVariantsTop + 1)[1]))
  
  
}

simulateStudyAF<-function(){
  studynumber<<-studynumber + 1
  print(cat("Beginning Study ", studynumber))
  print("Assigning betas...")
  assignVariantFrequencies()
  assignVariantBetaVariances()
  assignVariantBetas()
  
  print("Creating individuals...")
  createIndividuals()
  createPRSEstimates()
  
  print("Estimating...")
  c(unlist(estimateA(numVariantsBottom + numVariantsTop + 1)), unlist(estimateF(numVariantsBottom + numVariantsTop + 1)[1]))
  
  
}


simulateStudyBF<-function(){
  studynumber<<-studynumber + 1
  print(cat("Beginning Study ", studynumber))
  print("Assigning betas...")
  assignVariantFrequencies()
  assignVariantBetaVariances()
  assignVariantBetas()
  
  print("Creating individuals...")
  createIndividuals()
  createPRSEstimates()
  
  print("Estimating...")
  c(unlist(estimateB(numVariantsBottom + numVariantsTop + 1)), unlist(estimateF(numVariantsBottom + numVariantsTop + 1)))
  
  
}

simulateStudyABF<-function(){
  studynumber<<-studynumber + 1
  print(cat("Beginning Study ", studynumber))
  print("Assigning betas...")
  assignVariantFrequencies()
  assignVariantBetaVariances()
  assignVariantBetas()
  
  print("Creating individuals...")
  createIndividuals()
  createPRSEstimates()
  
  print("Estimating...")
  c(unlist(estimateA(numVariantsBottom + numVariantsTop + 1)[1]), unlist(estimateB(numVariantsBottom + numVariantsTop + 1)[1]), unlist(estimateF(numVariantsBottom + numVariantsTop + 1)[1]))
  
}

simulatePRSMean<-function(){
  studynumber<<-studynumber + 1
  print(cat("Beginning Study ", studynumber))
  print("Assigning betas...")
  assignVariantFrequencies()
  assignVariantBetaVariances()
  assignVariantBetas()
  
  print("Creating individuals...")
  createIndividuals()
  createPRSEstimates()
  
  print("Estimating...")
  mean(individualPRSs)
  
}

# plotVaryingBeta<-function(){
#   allp = c(unlist(info.h0[3,]), unlist(info.h10[3,]), unlist(info.h25[3,]), unlist(info.h50[3,]), unlist(info.h75[3,]), unlist(info.h90[3,]), unlist(info.h100[3,]))
#   allb = c(unlist(info.h0[1,]), unlist(info.h10[1,]), unlist(info.h25[1,]), unlist(info.h50[1,]), unlist(info.h75[1,]), unlist(info.h90[1,]), unlist(info.h100[1,]))
#   plot(log(allp), allb, col=allcolors, xlab="log p-Value", ylab="Beta Estimate", main="p-Beta Comparison for Heritability Explained (MAF = 0.005)", sub="Fraction Explained = 0%, 10%, 25%, 50%, 75%, 90%, 100%")
#   legend(-17, 0.02, c("0%", "10%", "25%", "50%", "75%", "90%", "100%"), lty=c(1, 1, 1, 1, 1, 1, 1), lwd=c(2, 2, 2, 2, 2, 2, 2), col=c("red", "orange", "yellow", "green", "blue", "purple", "pink"))
#   allcolors = c(rep("red", 100), rep("orange", 100), rep("yellow", 100), rep("green", 100), rep("blue", 100), rep("purple", 100), rep("pink", 100))
# 
# }
# 
# plotVaryingMAF<-function(){
#   allp = c(unlist(infob.h10[3,]), unlist(infob.h25[3,]), unlist(infob.h50[3,]), unlist(infob.h75[3,]), unlist(infob.h90[3,]), unlist(infob.h100[3,]))
#   allb = c(unlist(infob.h10[1,]), unlist(infob.h25[1,]), unlist(infob.h50[1,]), unlist(infob.h75[1,]), unlist(infob.h90[1,]), unlist(infob.h100[1,]))
#   allcolors = c(rep("orange", 100), rep("yellow", 100), rep("green", 100), rep("blue", 100), rep("purple", 100), rep("pink", 100))
#   plot(log(allp), allb, col=allcolors, xlab="log p-Value", ylab="Beta Estimate", main="p-Beta Comparison for Heritability Explained (Beta = -2)", sub="Fraction Explained = 10%, 25%, 50%, 75%, 90%, 100%")
#   legend(-24, 0.005, c("10%", "25%", "50%", "75%", "90%", "100%"), lty=c(1, 1, 1, 1, 1, 1), lwd=c(2, 2, 2, 2, 2, 2), col=c("orange", "yellow", "green", "blue", "purple", "pink"))
# 
# }

# *************** Method BX ***************

# Method B on multiple variants


# how to estimate a beta
estimateBX<-function(loci){
  
  # select out healthy inds
  healthy<<-which(individualPhenotypes == 0)
  
  # aggregate data
  myData = data.frame(t(rbind(individualPRSEstimates[healthy], individualGenotypes[loci, healthy])))
  fit <- lm(myData[,1] ~ ., data=myData[,-1])
  return(cbind(-summary(fit)$coefficients[2:(length(loci) + 1), "Estimate"], summary(fit)$coefficients[2:(length(loci) + 1), "Pr(>|t|)"]))
  #print(PRSZero - PRSOne)
  
}

simulateStudyBX<-function(){
  studynumber<<-studynumber + 1
  print(cat("Beginning Study ", studynumber))
  print("Assigning betas...")
  assignVariantFrequencies()
  assignVariantBetaVariances()
  assignVariantBetas()
  
  print("Creating individuals...")
  createIndividualsBX()
  createPRSEstimates()
  
  print("Estimating...")
  myLoci = numVariantsTop + numVariantsBottom + 1:numVariantsRare
  myEstimates = matrix(estimateBX(myLoci), nrow=2)
  bEstimates = matrix(unlist(sapply(myLoci, estimateB)), nrow=3)[c(1, 3),]
  rbind(t(myEstimates), bEstimates)
  # returns a matrix:
  # first row is betas from BX
  # second row ps from BX
  # third row betas fom B
  # fourth row ps from B
  
}

plotBXBetas<-function(){
  plot(rep(1:10, 4), c(bxinfo[1,1,1:10], bxinfo[1,2,1:10], bxinfo[3,1,1:10], bxinfo[3,2,1:10]), col=c(rep("red", 20), rep("green", 20)))
  
}

plotBXPs<-function(){
  plot(rep(1:10, 4), c(bxinfo[2,1,1:10], bxinfo[2,2,1:10], bxinfo[4,1,1:10], bxinfo[4,2,1:10]), col=c(rep("red", 20), rep("green", 20)))
  
  
}

plotBXlogPs<-function(){
  plot(rep(1:10, 4), log(c(bxinfo[2,1,1:10], bxinfo[2,2,1:10], bxinfo[4,1,1:10], bxinfo[4,2,1:10])), col=c(rep("red", 20), rep("green", 20)))
  
  
}
