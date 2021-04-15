dataGenerate <- function(ORs,p,f,n_samples){
  #a function to generate data from OR to case/control data
  #dominance model:eg.AA+AG vs GG
  #p=P(disease) in population
  #f=frequency of SNP genotype without risk(=frequency of G,>0.5)
  #n_samples=numbers of case/control group(a vector with same size as ORs)   
  library(rootSolve)
  if(length(n_samples)!=length(ORs)){
    n_samples <- rep(n_samples,length.out=length(ORs))
  }
  model <- function(x,parms)c(F1=((1-x[1])*(1+x[2]))/(1-x[1]*(1+x[2]))-parms[1],
                              F2=parms[3]^2*x[1]+(1-parms[3]^2)*x[1]*(1+x[2])-parms[2])
  n_studies <- length(ORs)
  ORhat <- vector(mode="numeric",length=n_studies)
  lnORhat <- vector(mode="numeric",length=n_studies)
  lnORSEhat <- vector(mode="numeric",length=n_studies)
  caseUnexposed <- vector(mode="numeric",length=n_studies)
  caseExposed <- vector(mode="numeric",length=n_studies)
  controlUnexposed <- vector(mode="numeric",length=n_studies)
  controlExposed <- vector(mode="numeric",length=n_studies)
  for(i in 1:n_studies){
    ss <- multiroot(f=model,start=c(0.5,0.5),parms=c(ORs[i],p,f))
    alpha <- ss$root[1]
    theta <- ss$root[2]
    controlGene <- rbinom(n_samples[i],1,1-(1-alpha)*f^2/(1-p))  #GG=0,AA+AG=1
    caseGene <- rbinom(n_samples[i],1,1-f^2*alpha/p)
    controlSum <- table(controlGene)
    caseSum <- table(caseGene)
    ORhat[i] <- (caseSum[2]/caseSum[1])*(controlSum[1]/controlSum[2])
    lnORhat[i] <- log(ORhat[i])
    lnORSEhat[i] <- sqrt(1/caseSum[1]+1/caseSum[2]+1/controlSum[1]+1/controlSum[2])
    caseUnexposed[i] <- caseSum[[1]]
    caseExposed[i] <- caseSum[[2]]
    controlUnexposed[i] <- controlSum[[1]]
    controlExposed[i] <- controlSum[[2]]
  }
  dataSim <- data.frame(ORhat=ORhat,betas=lnORhat,ses=lnORSEhat,
                        caseUnexposed=caseUnexposed,
                        caseExposed=caseExposed,
                        controlUnexposed=controlUnexposed,
                        controlExposed=controlExposed)  #beta=lnOR,se=SE(lnOR)
  return(dataSim)
}
