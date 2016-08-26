# initialize study variables
heritability = 0.2
numVariants = 200
numKnownVariants = 100

baselinePrevalence = 0.005

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

baselineScore = psi(baselinePrevalence)

# variant rates from beta distribution
variantFractions = rbeta(numVariants, 1, 99)

# variant betas
variantBetaVariances = rep(heritability / numVariants, numVariants)
variantBetas = sqrt(variantBetaVariances / (2 * variantFractions * (1 - variantFractions)))
variantBetas = variantBetas * sample(c(-1, 1), size = numVariants, replace = TRUE, prob = c(0.5, 0.5))

# how to generate a genotype
generateGenotype<-function(){
  rbinom(numVariants, 2, variantFractions)
  
}

# how to take PRS of genotype
PRS<-function(genotype){
  sum(genotype * variantBetas)
  
}

# how to get total risk score
riskScore<-function(genotype){
  baselineScore + PRS(genotype)
  
}

# generate genotypes
individualGenotypes = replicate(numIndividuals, generateGenotype())

# generate our riskScores
individualScores = apply(individualGenotypes, 2, riskScore)
individualPRSs = individualScores - baselineScore

# generate phenotypes
individualPhenotypes = rbinom(numIndividuals, 1, phi(individualScores))



# estimate our Beta values
betaEstimates = rnorm(numKnownVariants, variantBetas[1:numKnownVariants], 0.01)

# how to estimate PRS of genotype
estimatePRS<-function(genotype){
  sum(genotype[1:numKnownVariants] * betaEstimates)
  
}

# estimate everyone's PRS
individualPRSEstimates = apply(individualGenotypes, 2, estimatePRS)

# *************** Analysis ***************

# in analysis, we can use \beta_{m+1} through \beta_n as seperate attempts 
# since we only estimate one at a time


# *************** Method A ***************

# First attempt: measure beta by finding disease prevalence in those with and without
# no information needed other than phenotypes and G_{m+1}

# initialize our estimate vector
unknownEstimatesA = rep(NaN, numVariants - numKnownVariants)

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
for(i in (numKnownVariants + 1):numVariants){
  unknownEstimatesA[i - numKnownVariants] = estimateA(i)
  
}


# *************** Method B ***************

# select individuals without the disease and compare estimated PRS means
# in those with 1 allele and with 0

# initialize our estimate vector
unknownEstimatesB = rep(NaN, numVariants - numKnownVariants)


# how to estimate a beta
estimateB<-function(locus){
  
  # select out healthy inds
  healthy = which(individualPhenotypes == 0)
  
  # select those with 1 and 0 alleles, as those with 2 are too few
  gZero = which(individualGenotypes[locus, healthy] == 0)
  gOne = which(individualGenotypes[locus, healthy] == 1)
  
  # find PRS estimate means on populations
  PRSZero = mean(individualPRSEstimates[gZero])
  PRSOne = mean(individualPRSEstimates[gOne])
  
  # difference should be beta
  PRSOne - PRSZero
  
}

# estimate for each unknown beta
for(i in (numKnownVariants + 1):numVariants){
  unknownEstimatesB[i - numKnownVariants] = estimateB(i)
  
}


# *************** Method C ***************

# select individuals with the disease and compare estimated PRS means
# in those with 1 allele and with 0
# same as B except choosing those with

# initialize our estimate vector
unknownEstimatesC = rep(NaN, numVariants - numKnownVariants)


# how to estimate a beta
estimateC<-function(locus){
  
  # select out healthy inds
  healthy = which(individualPhenotypes == 1)
  
  # select those with 1 and 0 alleles, as those with 2 are too few
  gZero = which(individualGenotypes[locus, healthy] == 0)
  gOne = which(individualGenotypes[locus, healthy] == 1)
  
  # find PRS estimate means on populations
  PRSZero = mean(individualPRSEstimates[gZero])
  PRSOne = mean(individualPRSEstimates[gOne])
  
  # difference should be beta
  PRSOne - PRSZero
  
}

# estimate for each unknown beta
for(i in (numKnownVariants + 1):numVariants){
  unknownEstimatesC[i - numKnownVariants] = estimateC(i)
  
}


