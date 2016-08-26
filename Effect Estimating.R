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





