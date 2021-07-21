exh.abf<-function(betas,ses,prior.sigma,prior.cor="indep",prior.rho=NA,cryptic.cor=NA,log=FALSE,log10=FALSE,study.names=NULL,na.rm=FALSE,tolerance=1e-1000){
  #betas is a numeric vector of observed effect sizes of a single SNP in a set of studies. Each element of the vector is assumed to correspond to a study.
  #ses is a numeric vector of standard errors corresponding to those in betas. It should have the same length as betas.
  #prior.sigma is the prior on true effect sizes for each SNP in each study. It can be a flat value, set for each study (i.e. a vector whose length is equal to the number of studies in the meta-analysis) or set for each study and SNP (i.e. a matrix of same dimension as betas).
  #prior.cor is a square matrix whose row and column numbers are the same as the number of studies. Its elements are the pairwise correlations between true effect sizes of the studies. It can take values "indep" (independent effects), "fixed" (fixed effects), "correlated" (correlated effects, which requires the prior.rho parameter to be set), as well as individual matrices. If betas and ses are matrices, the same prior.cor will be applied to every row (representing every SNP).
  #prior.rho is either a single value or the upper triangle of a correlation matrix for the prior.cor matrix when it is set to "correlated". If this value is set, but prior.cor is not set to "correlated", this parameter will be ignored.
  #cryptic.cor a square matrix whose row and coumn numbers are the same as the number of studies. If the studies in the meta-analysis are not independent of each other, it may be necessary to set this parameter so that the covariance in null effects is accounted for.
  #log sets whether the answer should be given in log space.
  #log10 sets whether the answer should be given in log10 space.
  #study.names if set, the output will label the columns with the study names.
  #na.rm if there are NAs in the data, these are removed and the calculation is performed with the remaining data. This happens regardless of how this parameter is set. By default, the output will include a column of NAs for the study with the missing data. Changing this parameter to TRUE removes this column.
  #tolerance for the ABF calculation, this can be lowered (or raised, if necessary) if the answers are not what was expected. Should probably never be altered, but is there in case it is needed.
  
  #If betas and ses are data frames, this checks if they can be turned into numeric vectors. Stops the calculation if this is not the case.
  library(MASS)
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

mcmc.abf<-function(betas,ses,prior.sigma,prior.cor="indep",prior.rho=NA,cryptic.cor=NA,log=FALSE,log10=FALSE,na.rm=FALSE,tolerance=1e-1000,n.iter=500){
  #betas is a numeric vector of observed effect sizes of a single SNP in a set of studies. Each element of the vector is assumed to correspond to a study.
  #ses is a numeric vector of standard errors corresponding to those in betas. It should have the same length as betas.
  #prior.sigma is the prior on true effect sizes for each SNP in each study. It can be a flat value, set for each study (i.e. a vector whose length is equal to the number of studies in the meta-analysis) or set for each study and SNP (i.e. a matrix of same dimension as betas).
  #prior.cor is a square matrix whose row and column numbers are the same as the number of studies. Its elements are the pairwise correlations between true effect sizes of the studies. It can take values "indep" (independent effects), "fixed" (fixed effects), "correlated" (correlated effects, which requires the prior.rho parameter to be set), as well as individual matrices. If betas and ses are matrices, the same prior.cor will be applied to every row (representing every SNP).
  #prior.rho is either a single value or the upper triangle of a correlation matrix for the prior.cor matrix when it is set to "correlated". If this value is set, but prior.cor is not set to "correlated", this parameter will be ignored.
  #cryptic.cor a square matrix whose row and coumn numbers are the same as the number of studies. If the studies in the meta-analysis are not independent of each other, it may be necessary to set this parameter so that the covariance in null effects is accounted for.
  #log sets whether the answer should be given in log space.
  #log10 sets whether the answer should be given in log10 space.
  #study.names if set, the output will label the columns with the study names.
  #na.rm if there are NAs in the data, these are removed and the calculation is performed with the remaining data. This happens regardless of how this parameter is set. By default, the output will include a column of NAs for the study with the missing data. Changing this parameter to TRUE removes this column.
  #tolerance for the ABF calculation, this can be lowered (or raised, if necessary) if the answers are not what was expected. Should probably never be altered, but is there in case it is needed.
  
  #If betas and ses are data frames, this checks if they can be turned into numeric vectors. Stops the calculation if this is not the case.
  
  library(MASS)
  
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
  
  getModel<-function(ind){
    x <- rep(0,n)
	  x[ind] <- 1
	  return(x)
  }

  findNeighborPlus<-function(Curr,Non){
    NeighborP <- matrix(nrow=length(Non),ncol=n)
    for(i in 1:length(Non)){
      plu <- Non[i]
      indP <- c(Curr,plu)
      NeighborP[i,] <- getModel(indP)
    }
    return(NeighborP)
  }

  findNeighborMinus<-function(Curr){
    NeighborM <- matrix(nrow=length(Curr),ncol=n)
    for(i in 1:length(Curr)){
      min <- Curr[i]
      indP <- setdiff(Curr,min)
      NeighborM[i,] <- getModel(indP)
    }
    return(NeighborM)
  }  

  findNeighbor<-function(x){
    Curr <- which(x==1)
    Non <- setdiff(seq(1,n),Curr)
    if(length(Curr)==0){
      neighbors <- rbind(findNeighborPlus(Curr,Non),x)
    }else if(length(Non)==0){
      neighbors <- rbind(findNeighborMinus(Curr),x)
    }else{
      neighbors <- rbind(findNeighborPlus(Curr,Non),findNeighborMinus(Curr),x)
    }
    return(neighbors)
  }
  
  getABF<-function(x){
    prior.V<-matrix(x,nrow=n,ncol=n)*t(matrix(x,nrow=n,ncol=n))*prior.V.gen
    A<-prior.V+V
    invA<-ginv(A,tol=tolerance)
    quad.form<-t(b) %*% (ginv(V,tol=tolerance) - invA) %*% b
    lbf<-(-0.5*(as.numeric(determinant(A,logarithm=TRUE)$modulus)-null.calc))
    ABF<-(lbf + 0.5 * quad.form)
    if(ABF<=0){
      ABF<-10^(-20)
    }
    return(ABF)
  }
  
  k <- sample(1:n,1)
  xind <- sample(1:n,k)
  x <- getModel(xind)
  neighbors <- findNeighbor(x)
  ABFcurr <- getABF(x)
  
  modelist <- list(x)
  ABFlist <- c(ABFcurr)
  for(i in 1:n.iter){
    yI <- sample(nrow(neighbors),1,replace=TRUE)
    y <- neighbors[yI,]
    ABFy <- getABF(y)
    h <- min(1,ABFy/ABFcurr)
    r <- runif(1)
    if(r<h){
      x <- y
      ABFcurr <- ABFy
    }
    neighbors <- findNeighbor(x)
    ABFlist <- append(ABFlist,ABFcurr)
    modelist <- append(modelist,x)		
  }
  
  getmode <- function(v){
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  ABFfinall<-getmode(tail(ABFlist,ceiling(n.iter/2)))
  
  if(log10){
    ABFfinall<-ABFfinall/log(10)
  }
  if(!log && !log10){
    ABFfinall<-exp(ABFfinall)
  }
  if(!na.rm && n<nstudies){
    subset<-cbind(subset,matrix(NA,nrow=nmodels,ncol=(nstudies-n)))
    subset[,ind]<-subset[,1:n]
    subset[,nind]<-NA
  }
  
  return(ABFfinall)
}

shotgun.abf<-function(betas,ses,prior.sigma,prior.cor="indep",prior.rho=NA,cryptic.cor=NA,log=FALSE,log10=FALSE,na.rm=FALSE,tolerance=1e-1000,n.iter=50,B=5){
  #betas is a numeric vector of observed effect sizes of a single SNP in a set of studies. Each element of the vector is assumed to correspond to a study.
  #ses is a numeric vector of standard errors corresponding to those in betas. It should have the same length as betas.
  #prior.sigma is the prior on true effect sizes for each SNP in each study. It can be a flat value, set for each study (i.e. a vector whose length is equal to the number of studies in the meta-analysis) or set for each study and SNP (i.e. a matrix of same dimension as betas).
  #prior.cor is a square matrix whose row and column numbers are the same as the number of studies. Its elements are the pairwise correlations between true effect sizes of the studies. It can take values "indep" (independent effects), "fixed" (fixed effects), "correlated" (correlated effects, which requires the prior.rho parameter to be set), as well as individual matrices. If betas and ses are matrices, the same prior.cor will be applied to every row (representing every SNP).
  #prior.rho is either a single value or the upper triangle of a correlation matrix for the prior.cor matrix when it is set to "correlated". If this value is set, but prior.cor is not set to "correlated", this parameter will be ignored.
  #cryptic.cor a square matrix whose row and coumn numbers are the same as the number of studies. If the studies in the meta-analysis are not independent of each other, it may be necessary to set this parameter so that the covariance in null effects is accounted for.
  #log sets whether the answer should be given in log space.
  #log10 sets whether the answer should be given in log10 space.
  #study.names if set, the output will label the columns with the study names.
  #na.rm if there are NAs in the data, these are removed and the calculation is performed with the remaining data. This happens regardless of how this parameter is set. By default, the output will include a column of NAs for the study with the missing data. Changing this parameter to TRUE removes this column.
  #tolerance for the ABF calculation, this can be lowered (or raised, if necessary) if the answers are not what was expected. Should probably never be altered, but is there in case it is needed.
  
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

  getModel<-function(ind){
    x <- rep(0,n)
	  x[ind] <- 1
	  return(x)
  }

  findNeighborPlus<-function(Curr,Non){
    NeighborP <- matrix(nrow=length(Non),ncol=n)
    for(i in 1:length(Non)){
      plu <- Non[i]
      indP <- c(Curr,plu)
      NeighborP[i,] <- getModel(indP)
    }
    return(NeighborP)
  }

  findNeighborMinus<-function(Curr){
    NeighborM <- matrix(nrow=length(Curr),ncol=n)
    for(i in 1:length(Curr)){
      min <- Curr[i]
      indP <- setdiff(Curr,min)
      NeighborM[i,] <- getModel(indP)
    }
    return(NeighborM)
  }

  findNeighborO<-function(Curr,Non){
    NeighborO <- matrix(nrow=length(Curr)*length(Non),ncol=n)
    for(i in 1:length(Curr)){
      used <- Curr[i]
      indTep <- setdiff(Curr,used)
      for(j in 1:length(Non)){
        new <- Non[j]
        indO <- c(indTep,new)
        NeighborO[(i-1)*length(Non)+j,] <- getModel(indO)
      }
    }
    return(NeighborO)
  }  

  findNeighbor<-function(x){
    Curr <- which(x==1)
    Non <- setdiff(seq(1,n),Curr)
    if(length(Curr)==0){
      neighbors <- list(vector(),findNeighborPlus(Curr,Non),vector())
    }else if(length(Non)==0){
      neighbors <- list(vector(),vector(),findNeighborMinus(Curr))
    }else{
      neighbors <- list(findNeighborO(Curr,Non),findNeighborPlus(Curr,Non),findNeighborMinus(Curr))
    }
    return(neighbors)
  }
  
  getABF<-function(x){
    prior.V<-matrix(x,nrow=n,ncol=n)*t(matrix(x,nrow=n,ncol=n))*prior.V.gen
    A<-prior.V+V
    invA<-ginv(A,tol=tolerance)
    quad.form<-t(b) %*% (ginv(V,tol=tolerance) - invA) %*% b
    lbf<-(-0.5*(as.numeric(determinant(A,logarithm=TRUE)$modulus)-null.calc))
    ABF<-(lbf + 0.5 * quad.form)
    if(ABF<=0){
      ABF<-10^(-20)
    }
    return(ABF)
  }
  
  modelToABF<-function(models){
    abf<-apply(models,1,getABF)
    return(abf)
  }
  
  modelToString<-function(models){
    return(apply(models,1,paste,sep="",collapse=""))
  }
  
  k <- sample(1:n,1)
  xind <- sample(1:n,k)
  x <- getModel(xind)
  df1 <- data.frame(model=paste(x,sep="",collapse=""),abf=getABF(x))
  df1$model <- as.character(df1$model)
  dflist <- list(df1)
  
  for(i in 1:n.iter){
    neighbors <- findNeighbor(x)
    neighborO <- neighbors[[1]]
    if(length(neighborO)==0){
      dfO <- data.frame()
      IO <- NA
    }else if(nrow(neighborO)==1){
      dfO <- data.frame(model=modelToString(neighborO),abf=modelToABF(neighborO))
      dfO$model <- as.character(dfO$model)
      IO <- 1
    }else{
      dfO <- data.frame(model=modelToString(neighborO),abf=modelToABF(neighborO))
      dfO$model <- as.character(dfO$model)
      IO <- sample(nrow(dfO),1,prob=dfO$abf/sum(dfO$abf),replace=TRUE)
    }
    
    neighborPlus <- neighbors[[2]]
    if(length(neighborPlus)==0){
      dfPlus <- data.frame()
      IPlus <- NA
    }else if(nrow(neighborPlus)==1){
      dfPlus <- data.frame(model=modelToString(neighborPlus),abf=modelToABF(neighborPlus))
      dfPlus$model <- as.character(dfPlus$model)
      IPlus <- 1
    }else{
      dfPlus <- data.frame(model=modelToString(neighborPlus),abf=modelToABF(neighborPlus))
      dfPlus$model <- as.character(dfPlus$model)
      IPlus <- sample(nrow(dfPlus),1,prob=dfPlus$abf/sum(dfPlus$abf),replace=TRUE)
    }
    
    neighborMinus <- neighbors[[3]]
    if(length(neighborMinus)==0){
      dfMinus <- data.frame()
      IMinus <- NA
    }else if(nrow(neighborMinus)==1){
      dfMinus <- data.frame(model=modelToString(neighborMinus),abf=modelToABF(neighborMinus))
      dfMinus$model <- as.character(dfMinus$model)
      IMinus <- 1
    }else{
      dfMinus <- data.frame(model=modelToString(neighborMinus),abf=modelToABF(neighborMinus))
      dfMinus$model <- as.character(dfMinus$model)
      IMinus <- sample(nrow(dfMinus),1,prob=dfMinus$abf/sum(dfMinus$abf),replace=TRUE)
    }
    
    df1 <- rbind(df1,dfO,dfMinus,dfPlus)
    df1$model <- as.character(df1$model)
    df1 <- df1[!duplicated(df1),]
    df1 <- arrange(df1,abf)
    if(nrow(df1)>B){
      df1<-df1[-c(1:(nrow(df1)-B)),]
    }
    dflist <- c(dflist,list(df1))  
    
    mod3 <- c(dfO[IO,1],dfPlus[IPlus,1],dfMinus[IMinus,1])
    mod3 <- na.omit(mod3)
    abf3 <- c(dfO[IO,2],dfPlus[IPlus,2],dfMinus[IMinus,2])
    abf3 <- na.omit(abf3)
    if(length(mod3)==0){
      break
    }else if(length(mod3)==1){
      xstr <- mod3[1]
      x <- as.numeric(strsplit(xstr,split="")[[1]])
    }else{
      k <- sample(length(mod3),1,prob=abf3/sum(abf3),replace=TRUE)
      xstr <- mod3[k]
      x <- as.numeric(strsplit(xstr,split="")[[1]])
    }
  }
  
  ABFfinall<-df1[nrow(df1),2]
  
  if(log10){
    ABFfinall<-ABFfinall/log(10)
  }
  if(!log && !log10){
    ABFfinall<-exp(ABFfinall)
  }
  
  return(ABFfinall)
}

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

simulaData <- function(trueOR,n_studies,p,f,n_samples_lower,
                       n_samples_upper,VarOR=0.01){
  simDa <- function(trueORi){
    ORs <- rnorm(n_studies,trueORi,VarOR)
    n_samples <- sample(n_samples_lower:n_samples_upper,n_studies)
    dataSim <- dataGenerate(ORs,p,f,n_samples)
    return(dataSim)
  }
  dataList <- lapply(trueOR,simDa)
  return(dataList)
}

library(MASS)
library(dplyr)

dataList <- simulaData(0.8,20,0.05,0.8,100,2000,0.01)

dfs <- dataList[[1]]
betas <- dfs$betas
ses <- dfs$ses
repeats <- 100

timestart <- Sys.time()
maxABF <- max(exh.abf(betas,ses,0.5,"indep",log=TRUE)$ABF)
timeend <- Sys.time()
maxABFtime <- difftime(timeend,timestart,units="min")
print("maxABF")
c(maxABF,maxABFtime)

timestart <- Sys.time()
mcmc1 <-replicate(repeats,mcmc.abf(betas,ses,0.5,"indep",NA,log=TRUE,n.iter=1000))
timeend <- Sys.time()
mcmcABFlist1 <- mean(mcmc1)
mcmcTimelist1 <- (difftime(timeend,timestart,units="min"))/repeats
print("mcmc1")
mcmcABFlist1
mcmcTimelist1

timestart <- Sys.time()
mcmc2 <-replicate(repeats,mcmc.abf(betas,ses,0.5,"indep",NA,log=TRUE,n.iter=5000))
timeend <- Sys.time()
mcmcABFlist2 <- mean(mcmc2)
mcmcTimelist2 <- (difftime(timeend,timestart,units="min"))/repeats
print("mcmc2")
mcmcABFlist2
mcmcTimelist2

timestart <- Sys.time()
mcmc3 <-replicate(repeats,mcmc.abf(betas,ses,0.5,"indep",NA,log=TRUE,n.iter=10000))
timeend <- Sys.time()
mcmcABFlist3 <- mean(mcmc3)
mcmcTimelist3 <- (difftime(timeend,timestart,units="min"))/repeats
print("mcmc3")
mcmcABFlist3
mcmcTimelist3

timestart <- Sys.time()
mcmc4 <-replicate(repeats,mcmc.abf(betas,ses,0.5,"indep",NA,log=TRUE,n.iter=20000))
timeend <- Sys.time()
mcmcABFlist4 <- mean(mcmc4)
mcmcTimelist4 <- (difftime(timeend,timestart,units="min"))/repeats
print("mcmc4")
mcmcABFlist4
mcmcTimelist4

timestart <- Sys.time()
mcmc5 <-replicate(repeats,mcmc.abf(betas,ses,0.5,"indep",NA,log=TRUE,n.iter=50000))
timeend <- Sys.time()
mcmcABFlist5 <- mean(mcmc5)
mcmcTimelist5 <- (difftime(timeend,timestart,units="min"))/repeats
print("mcmc5")
mcmcABFlist5
mcmcTimelist5

timestart <- Sys.time()
shotgun1 <-replicate(repeats,shotgun.abf(betas,ses,0.5,"indep",NA,log=TRUE,n.iter=100,B=5))
timeend <- Sys.time()
shotgunABFlist1 <- mean(shotgun1)
shotgunTimelist1 <- (difftime(timeend,timestart,units="min"))/repeats
print("shotgun1")
shotgunABFlist1
shotgunTimelist1

timestart <- Sys.time()
shotgun2 <-replicate(repeats,shotgun.abf(betas,ses,0.5,"indep",NA,log=TRUE,n.iter=200,B=5))
timeend <- Sys.time()
shotgunABFlist2 <- mean(shotgun2)
shotgunTimelist2 <- (difftime(timeend,timestart,units="min"))/repeats
print("shotgun2")
shotgunABFlist2
shotgunTimelist2

timestart <- Sys.time()
shotgun3 <-replicate(repeats,shotgun.abf(betas,ses,0.5,"indep",NA,log=TRUE,n.iter=500,B=5))
timeend <- Sys.time()
shotgunABFlist3 <- mean(shotgun3)
shotgunTimelist3 <- (difftime(timeend,timestart,units="min"))/repeats
print("shotgun3")
shotgunABFlist3
shotgunTimelist3

timestart <- Sys.time()
shotgun4 <-replicate(repeats,shotgun.abf(betas,ses,0.5,"indep",NA,log=TRUE,n.iter=1000,B=5))
timeend <- Sys.time()
shotgunABFlist4 <- mean(shotgun4)
shotgunTimelist4 <- (difftime(timeend,timestart,units="min"))/repeats
print("shotgun4")
shotgunABFlist4
shotgunTimelist4

timestart <- Sys.time()
shotgun5 <-replicate(repeats,shotgun.abf(betas,ses,0.5,"indep",NA,log=TRUE,n.iter=2000,B=5))
timeend <- Sys.time()
shotgunABFlist5 <- mean(shotgun5)
shotgunTimelist5 <- (difftime(timeend,timestart,units="min"))/repeats
print("shotgun5")
shotgunABFlist5
shotgunTimelist5

dfSum<-data.frame(exhV=maxABF1,exhT=maxABFtime1,mcmcV1=mcmcABFlist1,mcmcT1=mcmcTimelist1,
                  mcmcV2=mcmcABFlist2,mcmcT2=mcmcTimelist2,mcmcV3=mcmcABFlist3,mcmcT3=mcmcTimelist3,
                  mcmcV4=mcmcABFlist4,mcmcT4=mcmcTimelist4,mcmcV5=mcmcABFlist5,mcmcT5=mcmcTimelist5,
                  shotgunV1=shotgunABFlist1,shotgunT1=shotgunTimelist1,shotgunV2=shotgunABFlist2,shotgunT2=shotgunTimelist2,
                  shotgunV3=shotgunABFlist3,shotgunT3=shotgunTimelist3,shotgunV4=shotgunABFlist4,shotgunT4=shotgunTimelist4,
                  shotgunV5=shotgunABFlist5,shotgunT5=shotgunTimelist5)
write.csv(dfSum,"sim3_sum.csv")
dfDetail <- rbind(mcmc1,mcmc2,mcmc3,mcmc4,mcmc5,shotgun1,shotgun2,shotgun3,shotgun4,shotgun5)
write.csv(dfDetail,"sim3_detail.csv")