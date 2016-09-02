lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p[1])
}


# initialize study variables
heritability = 0.3
heritabilityTopProp = .9
heritabilityBottomProp = 0
heritabilityRareProp = 1 - heritabilityTopProp - heritabilityBottomProp
numVariantsTop = 200
numVariantsBottom = 0
numVariantsRare = 2
numStrongProtVariants = 5
numStrongRiskVariants = 5
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
  variantFrequenciestop<<-rbeta(numVariantsTop, 1.5, 1.5)
  # variant frequencies bottom 
  variantFrequenciesbottom<<-rbeta(numVariantsBottom, 1.5, 1.5)
  # variant frequencies rare
  variantFrequenciesrare<<-c(.005,.005)
  
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

# estimate for each unknown beta
for(i in (numVariantsTop+numVariantsBottom+1):(numVariantsTop+numVariantsBottom+numVariantsRare)){
  unknownEstimatesA[i-numVariantsTop - numVariantsBottom] = estimateA(i)
  
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

simulateStudy<-function(){
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
  tailsize = 0.1
  
  cutoffs = quantile(individualPRSEstimates[healthy], c(tailsize, 1-tailsize))
  topPRSPeople = healthy[which(individualPRSEstimates > cutoffs[2])]
  bottomPRSPeople = healthy[which(individualPRSEstimates < cutoffs[1])]
  allPeople = c(bottomPRSPeople, topPRSPeople)
  
  # select indices for each number of alleles
  gZero = allPeople[which(individualGenotypes[locus,] == 0)]
  gOne = allPeople[which(individualGenotypes[locus,] == 1)]
  gTwo = allPeople[which(individualGenotypes[locus,] == 2)]
  
  # find disease prevalence in these groups
  zeroAllelePrevalence = length(which(gZero %in% topPRSPeople)) / length(gZero)
  oneAllelePrevalence = length(which(gOne %in% topPRSPeople)) / length(gOne)
  twoAllelePrevalence = length(which(gOne %in% topPRSPeople)) / length(gTwo)
  
  # use prevalence to estimate beta
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
  
  estimateC(numVariantsBottom + numVariantsTop + 1)
  
  
}






