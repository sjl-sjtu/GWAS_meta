#' @title ABF calculation for GWAS meta-analysis through Subset-Exhaustive
#'
#' @description Subset-Exhaustive to get the ABFs in all subsets for a single SNP and a phenotype in GWAS meta-analysis.
#'
#' @param betas a numeric vector of observed effect sizes of a single SNP in a set of studies. Each element of the vector is assumed to correspond to a study.
#' @param ses a numeric vector of standard errors corresponding to those in betas. It should have the same length as betas.
#' @param prior.sigma the prior on true effect sizes for each SNP in each study. It can be a flat value, set for each study (i.e. a vector whose length is equal to the number of studies in the meta-analysis) or set for each study and SNP (i.e. a matrix of same dimension as betas).
#' @param prior.cor a square matrix whose row and column numbers are the same as the number of studies. Its elements are the pairwise correlations between true effect sizes of the studies. It can take values "indep" (independent effects), "fixed" (fixed effects), "correlated" (correlated effects, which requires the prior.rho parameter to be set), as well as individual matrices. If betas and ses are matrices, the same prior.cor will be applied to every row (representing every SNP).
#' @param prior.rho either a single value or the upper triangle of a correlation matrix for the prior.cor matrix when it is set to "correlated". If this value is set, but prior.cor is not set to "correlated", this parameter will be ignored.
#' @param cryptic.cor a square matrix whose row and coumn numbers are the same as the number of studies. If the studies in the meta-analysis are not independent of each other, it may be necessary to set this parameter so that the covariance in null effects is accounted for.
#' @param log sets whether the answer should be given in log space.
#' @param log10 sets whether the answer should be given in log10 space.
#' @param study.names if set, the output will label the columns with the study names.
#' @param na.rm if there are NAs in the data, these are removed and the calculation is performed with the remaining data. This happens regardless of how this parameter is set. By default, the output will include a column of NAs for the study with the missing data. Changing this parameter to TRUE removes this column.
#' @param tolerance for the ABF calculation, this can be lowered (or raised, if necessary) if the answers are not what was expected. Should probably never be altered, but is there in case it is needed. Default is 1e-1000.
#'
#' @return a data frame containing ABF calculated in all subsets.
#' @export
#'
#' @examples
#'
exh_abf<-function(betas,ses,prior.sigma=0.3,prior.cor="indep",prior.rho=NA,cryptic.cor=NA,log=FALSE,log10=FALSE,study.names=NULL,na.rm=FALSE,tolerance=1e-1000){
  #If betas and ses are data frames, this checks if they can be turned into numeric vectors. Stops the calculation if this is not the case.
  library(MASS)
  library(dplyr)

  if(!class(betas) %in% c("data.frame","numeric")){
    stop("betas should be a numeric vector.")
  }
  if(!class(ses) %in% c("data.frame","numeric")){
    stop("ses should be a numeric vector.")
  }

  if(class(betas)=="data.frame"){
    if(dim(betas)[1]>1){
      stop("betas should be a vector, not multiple rows from a matrix or a data frame.")
    }
    if(all(apply(betas,2,class)=="numeric")){
      betas<-as.numeric(betas)
    } else {
      stop("betas is not numeric.")
    }
  }
  if(class(ses)=="data.frame"){
    if(dim(ses)[1]>1){
      stop("ses should be a vector, not multiple rows from a matrix or a data frame.")
    }
    if(all(apply(ses,2,class)=="numeric")){
      ses<-as.numeric(ses)
    } else {
      stop("ses is not numeric.")
    }
  }
  if(class(prior.sigma)=="data.frame"){
    if(dim(prior.sigma)[1]>1){
      stop("prior.sigma should be a vector or a single value, not multiple rows from a matrix or a data frame.")
    }
    if(all(apply(prior.sigma,2,class)=="numeric")){
      prior.sigma<-as.numeric(prior.sigma)
    } else {
      stop("prior.sigma is not numeric.")
    }
  }
  if(log && log10){
    stop("can give the approximate Bayes factor in log space or log_10 space, but not both at the same time.")
  }
  if(tolerance>sqrt(.Machine$double.eps)){
    warning(paste0("Your tolerance might be too high. The standard value for the internal functions is ",sqrt(.Machine$double.eps),"."))
  }

  if(length(betas)!=length(ses)){
    stop("betas and ses should be the same length.")
  }
  nstudies<-length(betas)

  if(length(study.names)>0){
    if(class(study.names)!="character"){
      stop("study.names needs to be a character vector.")
    }
    if(length(study.names)!=nstudies){
      stop(paste0("study.names needs to be of length ",nstudies,"."))
    }
  }

  ##Check prior.sigma isn't erroneous.
  if(!length(prior.sigma) %in% c(1,nstudies)){
    stop("prior.sigma should either be a single value or a vector whose length is equal to the number of studies.")
  }
  if(any(prior.sigma<=0)){
    stop("all values of prior.sigma must be > 0.")
  }
  if(length(prior.sigma)==1){
    prior.sigma<-rep(prior.sigma,nstudies)
  }

  ##Get the prior correlation matrix.
  if(class(prior.cor)=="character"){
    if(prior.cor=="correlated"){
      if(is.na(prior.rho)){
        stop("prior.rho must be set or prior.cor must be changed to something other than \"correlated.\"")
      }
      if(!is.numeric(prior.rho)){
        stop("prior.rho must be numeric.")
      }
      if(!length(prior.rho) %in% c(1,choose(nstudies,2))){
        stop("prior.rho should be either a flat value or the full upper triangle of the desired correlation matrix.")
      }
      if(length(prior.cor)==1){
        prior.cor.mat<-matrix(prior.rho,nrow=nstudies,ncol=nstudies)
        diag(prior.cor.mat)<-1
      } else {
        prior.cor.mat<-diag(nstudies)
        prior.cor.mat[upper.tri(prior.cor.mat,diag=FALSE)]<-prior.cor
        prior.cor.mat[lower.tri(prior.cor.mat)]<-t(prior.cor.mat)[lower.tri(prior.cor.mat)]
      }
    } else if(prior.cor=="indep"){
      prior.cor.mat<-diag(nstudies)
    } else if(prior.cor=="fixed"){
      prior.cor.mat<-matrix(1,nrow=nstudies,ncol=nstudies)
    }
  } else if(class(prior.cor)!="matrix"){
    stop("prior.cor should be a matrix, or one of the following values: \"indep\", \"fixed\", or \"correlated\"")
  } else {
    if(dim(prior.cor)[1]!=dim(prior.cor)[2]){
      stop("prior.cor is not a square matrix.")
    }
    if(dim(prior.cor)[1]!=nstudies){
      stop("the dimensions of prior.cor do not match the number of studies in the meta-analysis.")
    }
    if(!isSymmetric.matrix(prior.cor)){
      stop("prior.cor is not a symmetric matrix.")
    }
    if(any(svd(prior.cor)$d<0)){
      warning("prior.cor is not positive semidefinite. Errors may arise due to this.")
    }
    if(!all(diag(prior.cor) %in% c(0,1))){
      stop("the diagonal of prior.cor should be 1 for all studies with true effects and 0 elsewhere.")
    }
    if(any(diag(prior.cor)==0)){
      if(all(prior.cor[which(diag(prior.cor==0)),]!=0) || all(prior.cor[,which(diag(prior.cor==0))]!=0)){
        if(length(which(diag(prior.cor==0)))==1){
          stop(paste0("Row ",which(diag(prior.cor==0))," and column ",which(diag(prior.cor==0))," of prior.cor should all be 0, or else prior.cor[",which(diag(prior.cor==0)),",",which(diag(prior.cor==0)),"] should be 1."))
        } else {
          stop(paste0("Rows ",paste(which(diag(prior.cor==0)),collapse=",")," and columns ",paste(which(diag(prior.cor==0)),collapse=",")," of prior.cor should all be 0, or else prior.cor[c(",paste(which(diag(prior.cor==0)),collapse=","),"), c(",paste(which(diag(prior.cor==0)),collapse=","),")] should all be 1."))
        }
      }
    }
    prior.cor.mat<-prior.cor
  }

  ##Get the cryptic correlation matrix
  if(all(is.na(cryptic.cor)) && length(cryptic.cor)==1){
    cryptic.cor.mat<-diag(nstudies)
  } else if(class(cryptic.cor)!="matrix"){
    stop("cryptic.cor should be a square matrix.")
  } else if(dim(cryptic.cor)[1]!=dim(cryptic.cor)[2]){
    stop("cryptic.cor should be a square matrix.")
  } else if(dim(cryptic.cor)[1]!=nstudies){
    stop("the dimensions of cryptic.cor do not match the number of studies in the meta-analysis.")
  } else if(!isSymmetric.matrix(cryptic.cor)){
    stop("cryptic.cor is not a symmetric matrix.")
  } else if(any(eigen(cryptic.cor)$values<0)){
    stop("cryptic.cor is not positive definite.")
  } else if(any(diag(cryptic.cor)!=1)){
    stop("the diagonal of cryptic.cor should be 1 uniformly.")
  } else if(length(cryptic.cor)>1 && any(is.na(cryptic.cor))){
    stop("If cryptic.cor is defined, then no NAs are permitted.")
  } else {
    cryptic.cor.mat<-cryptic.cor
  }



  #Functions takes as input integer (or vector of integers) and integer b (default is b=2)
  #and outputs integer(s) in base b.
  integer.base.b <-function(x, b=2){
    xi <- as.integer(x)
    if(any(is.na(xi) | ((x-xi)!=0)))
      print(list(ERROR="x not integer", x=x))
    N <- length(x)
    xMax <- max(x)
    ndigits <- (floor(logb(xMax, base=2))+1)
    Base.b <- array(NA, dim=c(N, ndigits))
    for(i in 1:ndigits){#i <- 1
      Base.b[, ndigits-i+1] <- (x %% b)
      x <- (x %/% b)
    }
    if(N == 1) Base.b[1, ] else Base.b
  }


  #Calculate V and the null once.
  ind<-intersect(which(!is.na(betas)),which(!is.na(ses)))
  nind<-union(which(is.na(betas)),which(is.na(ses)))
  n<-length(ind)
  if(n<=1){
    return(NA)
  }

  b<-betas[ind]
  se<-ses[ind]
  ccm<-cryptic.cor.mat[ind,ind]
  prior.V.gen<-prior.cor.mat[ind,ind]*matrix(prior.sigma[ind],nrow=n,ncol=n)*t(matrix(prior.sigma[ind],nrow=n,ncol=n))

  V<-diag(se^2)
  for(i in 1:(nrow(V)-1)){
    for(j in (i+1):ncol(V)){
      V[i,j]<-sqrt(V[i,i])*sqrt(V[j,j])*ccm[i,j]
      V[j,i]<-V[i,j]
    }
  }

  null.calc<-as.numeric(determinant(V,logarithm=TRUE)$modulus)

  nmodels<-2^n
  subset<-integer.base.b(0:(nmodels-1))
  ABF <- rep(0,nmodels)
  for(i in 2:nmodels){
    prior.V<-matrix(subset[i,],nrow=n,ncol=n)*t(matrix(subset[i,],nrow=n,ncol=n))*prior.V.gen
    A<-prior.V+V
    invA<-ginv(A,tol=tolerance)
    quad.form<-t(b) %*% (ginv(V,tol=tolerance) - invA) %*% b
    lbf<-(-0.5*(as.numeric(determinant(A,logarithm=TRUE)$modulus)-null.calc))
    ABF[i]<-(lbf + 0.5 * quad.form)
  }

  if(log10){
    ABF<-ABF/log(10)
  }
  if(!log && !log10){
    ABF<-exp(ABF)
  }
  if(!na.rm && n<nstudies){
    subset<-cbind(subset,matrix(NA,nrow=nmodels,ncol=(nstudies-n)))
    subset[,ind]<-subset[,1:n]
    subset[,nind]<-NA
  }
  if(!is.null(study.names)){
    if(na.rm){
      colnames(subset)<-study.names[ind]
    } else {
      colnames(subset)<-study.names
    }
  }

  return(as.data.frame(cbind(subset,ABF)))
}
