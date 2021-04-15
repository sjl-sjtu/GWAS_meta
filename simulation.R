trueOR <- union(seq(0.9,1.1,0.02),seq(0.6,1.4,0.1))
L <- length(trueOR)
n_studies <- 20
pvaluelistFix <- vector(mode="numeric",length=L)
pvaluelistRan <- vector(mode="numeric",length=L)
maxABFlist <- vector(mode="numeric",length=L)
mcmcABFlist <- vector(mode="numeric",length=L)
shotgunABFlist <- vector(mode="numeric",length=L)
library(meta)
timestart<-Sys.time()
for(i in 1:L){
  ORs <- rnorm(n_studies,trueOR[i],0.01)
  n_samples <- sample(100:5000,n_studies)
  dataSim <- dataGenerate(ORs,0.05,0.8,n_samples)
  betas <- dataSim$betas
  ses <- dataSim$ses
  maxABFlist[i] <- max(exh.abf(betas,ses,0.5,"correlated",0.7,log=TRUE))
  mcmcABFlist[i] <- mcmc.abf(betas,ses,0.5,"correlated",0.7,log=TRUE,n.iter=500)
  shotgunABFlist[i] <- shotgun.abf(betas,ses,0.5,"correlated",0.7,log=TRUE,n.iter=50,B=5)
  metaRe <- metabin(caseExposed,controlExposed+caseExposed,
                    caseUnexposed,controlUnexposed+caseUnexposed,
                    data=dataSim,sm="OR")
  pvaluelistFix[i] <- metaRe$pval.fixed  #fixed effect model
  pvaluelistRan[i] <- metaRe$pval.random  #random effect model
}
timeend<-Sys.time()
runningtime<-timeend-timestart
dfTotal <- data.frame(trueOR=trueOR,pvalueFix=pvaluelistFix,
                 pvalueRan=pvaluelistRan,maxABF=maxABFlist,
                 mcmcABF=mcmcABFlist,shotgunABF=shotgunABFlist)
library(ggplot2)
boost <- 250
colours <- c("maxABF"="black","mcmcABF"="blue","shotgunABF"="red",
             "pvalueFix"="gold","pvalueRan"="green")
ggplot(dfTotal)+geom_line(aes(trueOR,maxABF,colour="maxABF"))+
  geom_line(aes(trueOR,mcmcABF,colour="mcmcABF"))+
  geom_line(aes(trueOR,shotgunABF,colour="shotgunABF"))+
  geom_line(aes(trueOR,pvalueFix*boost,colour="pvalueFix"))+
  geom_line(aes(trueOR,pvalueRan*boost,colour="pvalueRan"))+
  labs(x="true OR",y="ABF")+
  scale_y_continuous(sec.axis=sec_axis(~./boost,name="p-value"))+
  scale_colour_manual(name="class",values=colours)

