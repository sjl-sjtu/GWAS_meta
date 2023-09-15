#' Title ABF calculation across multiple variants for GWAS meta-analysis through shotgun stochastic search
#'
#' @description using shotgun stochastic search algorithm to quickly get the optimal ABF for multiple SNPs and a phenotype in GWAS meta-analysis.
#'
#' @param df a data frame containing the names of variants, the effect sizes and standard errors for each variant in each study.
#' @param vname the index of column of variant name in df. Default is 1.
#' @param vbetas a vectors to represent the index of columns of those containing the effect sizes in each study, default is the column 2,4,6,etc.
#' @param vses a vectors to represent the index of columns of those containing the corresponding standard errors in each study, default is the column 3,5,7,etc. The length of vses should be the same as the length of vbetas.
#' @param prior.sigma the prior on true effect sizes for each SNP in each study. It can be a flat value, set for each study (i.e. a vector whose length is equal to the number of studies in the meta-analysis) or set for each study and SNP (i.e. a matrix of same dimension as betas).
#' @param prior.cor a square matrix whose row and column numbers are the same as the number of studies. Its elements are the pairwise correlations between true effect sizes of the studies. It can take values "indep" (independent effects), "fixed" (fixed effects), "correlated" (correlated effects, which requires the prior.rho parameter to be set), as well as individual matrices. If betas and ses are matrices, the same prior.cor will be applied to every row (representing every SNP).
#' @param prior.rho either a single value or the upper triangle of a correlation matrix for the prior.cor matrix when it is set to "correlated". If this value is set, but prior.cor is not set to "correlated", this parameter will be ignored.
#' @param cryptic.cor a square matrix whose row and coumn numbers are the same as the number of studies. If the studies in the meta-analysis are not independent of each other, it may be necessary to set this parameter so that the covariance in null effects is accounted for.
#' @param log sets whether the answer should be given in log space.
#' @param log10 sets whether the answer should be given in log10 space.
#' @param na.rm if there are NAs in the data, these are removed and the calculation is performed with the remaining data. This happens regardless of how this parameter is set. By default, the output will include a column of NAs for the study with the missing data. Changing this parameter to TRUE removes this column.
#' @param tolerance for the ABF calculation, this can be lowered (or raised, if necessary) if the answers are not what was expected. Should probably never be altered, but is there in case it is needed. Default is 1e-1000.
#' @param n.iter the number of iteration for MCMC. Default is 50.
#' @param B the largest number of subsets involved in the optimal set. Default is 5.
#'
#' @return a data frame containing the following columns:
#' \item{SNP}{Variant Name}
#' \item{ABF}{the final ABF value calculated for each variant}
#' \item{model}{a string of 0-1 to represents the subset model selected to calculate.}
#' \item{n_studies}{the number of studies involved for each variant. If needClean is FALSE, the value will be the same for all variants.}
#' \item{studies_involved}{a string of 0-1 to represents the studies involved for each variant. 1 represents the corresponding study is involved.}
#' @export
#'
#' @examples
#' library(GWASmeta)
#' data(multi)
#' re <- multi_shotgun_abf(multi)
#'
multi_shotgun_abf <- function(df,vname=1,vbetas=seq(2,ncol(df),2),vses=seq(3,ncol(df),2),prior.sigma=0.5,
                              prior.cor="indep",prior.rho=NA,cryptic.cor=NA,log=FALSE,log10=FALSE,
                              na.rm=FALSE,tolerance=1e-1000,n.iter=50,B=5){
  library(MASS)
  library(dplyr)
  if(length(vbetas)!=length(vses)){
    return(message("Vectors betas and ses do not have the same length!"))
  }
  if(length(vbetas)==1){
    return(message("Only one study involved!"))
  }
  get_abf <- function(i){
    cali <- as.numeric(is.na(df[i,vbetas]) | is.na(df[i,vses]))
    cali <- 1-cali
    studiesUsed <- paste(cali,collapse="")
    nstudies <- sum(cali==1)
    SNP <- df[i,vname]
    betas <- df[i,vbetas]
    ses <- df[i,vses]
    abfi <- shotgun_abf_model(betas,ses,prior.sigma,prior.cor,prior.rho,
                             cryptic.cor=NA,log,log10,na.rm,tolerance=1e-1000,n.iter,B=5)
    abfvalue <- abfi$ABF
    submodel <- abfi$model
    return(c(SNP,abfvalue,submodel,nstudies,studiesUsed))
  }
  df[df==0] <- NA
  nstud <- length(vbetas)
  for(i in 1:nstud){
    df[which(is.na(df[vbetas[i]])),vses[i]] <- NA
    df[which(is.na(df[vses[i]])),vbetas[i]] <- NA
  }
  re <- sapply(seq(1,nrow(df)),get_abf)
  abf <- data.frame(SNP=re[1,],ABF=re[2,],model=re[3,],n_studies=re[4,],studies_involved=re[5,])
  abf$ABF <- as.numeric(abf$ABF)
  abf$ABF <- round(abf$ABF,4)
  abf <- arrange(abf,desc(ABF))
  return(abf)
}
