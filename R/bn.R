#' @title Getting suitable genetic IVs through Bayesian network learning
#'
#' @description is used to get the suitable SNPs as instrumental variables (IVs) of specified exposure by Bayesian network (BN) structure learning.
#'
#' @param df a data frame which contains data of SNPs and specified exposure.
#' @param snp a vector of string belonging to colnames of df, which is the name of SNPs included in BN structure learning.
#' @param exposureName a string which is a colname of df corresponding to the exposure studied.
#' @param bn_method method for BN structure learning. Possible values are the function name of structure learning algorithm implemented in bnlearn. Default is "hc".
#' @param cutoff a numeric between 0 to 1. Those SNPs with score larger than "cutoff" will be chosen as IVs. Default is 0.7.
#' @param repeats an integer standing for the times of bootstraping. Default is 100.
#' @param nsam an integer standing for the sample size for bootstraping sampling. Default is 500.
#'
#' @return a list containing:
#'   selectsnp: a vector of string containing the colnames of df corresponding to
#'   selected SNPs.
#'   dfscore: a data frame containing the score calculated for each SNP.
#' @export
#'
#' @examples
#'
#'
bn <- function(df,snp,exposureName,bn_method="hr",cutoff=0.7,repeats=100,nsam=500){
  library("bnlearn")

  learnBN <- function(df,nsam,bn_method){
    n <- nrow(df)
    iSam <- sample(seq(1,n),size = nsam,replace=TRUE)
    dfSam <- df[iSam,]
    if(bn_method=="pc.stable"){
      model <- pc.stable(dfSam)
    }else if(bn_method=="gs"){
      model <- gs(dfSam)
    }else if(bn_method=="iamb"){
      model <- iamb(dfSam)
    }else if(bn_method=="fast.iamb"){
      model <- fast.iamb(dfSam)
    }else if(bn_method=="inter.iamb"){
      model <- inter.iamb(dfSam)
    }else if(bn_method=="iamb.fdr"){
      model <- iamb.fdr(dfSam)
    }else if(bn_method=="hc"){
      model <- hc(dfSam)
    }else if(bn_method=="tabu"){
      model <- tabu(dfSam)
    }else if(bn_method=="mmhc"){
      model <- mmhc(dfSam)
    }else if(bn_method=="rsmax2"){
      model <- rsmax2(dfSam)
    }else if(bn_method=="h2pc"){
      model <- h2pc(dfSam)
    }else if(bn_method=="mmpc"){
      model <- mmpc(dfSam)
    }else if(bn_method=="si.hiton.pc"){
      model <- si.hiton.pc(dfSam)
    }else if(bn_method=="hpc"){
      model <- hpc(dfSam)
    }else if(bn_method=="chow.liu"){
      model <- chow.liu(dfSam)
    }else if(bn_method=="aracne"){
      model <- aracne(dfSam)
    }else{
      return(message("no this bn learning method"))
    }
    dfarc <- data.frame(model$arcs)
    return(dfarc)
  }

  rmBidire <- function(df){
    delL <- c()
    for(i in 1:nrow(df)){
      for(j in i+1:nrow(df)){
        if(all(df[j,]%in%df[i,])){
          delL <- append(delL,j)
        }
      }
    }
    df <- df[-delL,]
    return(df)
  }

  BNbootstrap <- function(df,repeats,nsam,bn_method){
    arcsL <- replicate(repeats,learnBN(df,nsam,bn_method),simplify = FALSE)
    library("plyr")
    arcsL <- do.call(rbind.fill,arcsL)
    arcsL <- rmBidire(arcsL)
    arcsL$from <- as.factor(arcsL$from)
    arcsL$to <-as.factor(arcsL$to)
    arcsL$count <- rep(1,nrow(arcsL))
    dfre <- aggregate(arcsL$count,by=list(arcsL$from,arcsL$to),FUN=sum)
    colnames(dfre) <- c("from","to","count")
    dfre <- arrange(dfre,-count)
    return(dfre)
  }

  getscore <- function(dfre,exposureName,snp,repeats){
    #exposureName is a str, snp is a vector of str.
    score <- rep(0,length(snp))
    for(i in 1:length(snp)){
      sn <- snp[i]
      count1 <- dfre[which(dfre$from==sn&dfre$to==exposureName),"count"]
      count2 <- dfre[which(dfre$from==exposureName&dfre$to==sn),"count"]
      if(length(count1)==0){
        count1 <- 0
      }
      if(length(count2)==0){
        count2 <- 0
      }
      score[i] <- (count1+count2)/repeats
    }
    dfscore <- data.frame(snp=snp,score=score)
    return(dfscore)
  }

  df1 <- df[,c(snp,exposureName)]
  dfsnp <- df[,snp]
  exposure <- df[,exposureName]

  dfre <- BNbootstrap(df1,repeats,nsam,bn_method)
  dfscore <- getscore(dfre,exposureName,snp,repeats)
  selectsnp <- dfscore[which(dfscore$score>=cutoff),"snp"]

  re <- list(IV=selectsnp,score=dfscore)
  return(re)
}
