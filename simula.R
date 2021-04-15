simula <- function(trueOR,n_studies,p,f,n_samples_lower,n_samples_upper,
                   prior.sigma,prior.cor="indep",prior.rho=NA,VarOR=0.01,
                   mcmcIter=500,shotgunIter=30,B=5,repeats=10){
  library(meta)
  simulaIn <- function(trueORi){
    ORs <- rnorm(n_studies,trueORi,VarOR)
    n_samples <- sample(n_samples_lower:n_samples_upper,n_studies)
    dataSim <- dataGenerate(ORs,p,f,n_samples)
    betas <- dataSim$betas
    ses <- dataSim$ses
    metaRe <- metabin(caseExposed,caseUnexposed+caseExposed,
                      controlExposed,controlUnexposed+controlExposed,
                      data=dataSim,sm="OR")
    pvalueFix <- metaRe$pval.fixed  #fixed effect model
    pvalueRan <- metaRe$pval.random  #random effect model
    maxABF <- max(exh.abf(betas,ses,prior.sigma,prior.cor,prior.rho,log=TRUE))
    mcmcABFset <- replicate(repeats,mcmc.abf(betas,ses,prior.sigma,prior.cor,prior.rho,log=TRUE,n.iter=mcmcIter))
    mcmcABF <- mean(mcmcABFset)
    shotgunABFset <- replicate(repeats,shotgun.abf(betas,ses,prior.sigma,prior.cor,prior.rho,log=TRUE,n.iter=shotgunIter,B=B))
    shotgunABF <- mean(shotgunABFset)
    return(c(trueORi,pvalueFix,pvalueRan,maxABF,mcmcABF,shotgunABF))
  }
  dfMatrix <- sapply(trueOR,simulaIn)
  df <- data.frame(trueOR=dfMatrix[1,],pvalueFix=dfMatrix[2,],
                   pvalueRan=dfMatrix[3,],maxABF=dfMatrix[4,],
                   mcmcABF=dfMatrix[5,],shotgunABF=dfMatrix[6,])
  return(df)
}
