#' @title Causal estimation by adaptive Bayesian Mendelian randomization with pleiotropy
#'
#' @description is used to get the causal estimation between exposure and outcome through Bayesian adaptive Mendelian randomization with pleiotropy after IV selection.
#'
#' @param df a data frame which contains data of IVs, specified exposure and outcome.
#' @param selectsnp a vector of string containing the colnames of df corresponding to the IV used in MR.
#' @param exposureName a string which is a colname of df corresponding to the exposure studied.
#' @param outcomeName a string which is a colname of df corresponding to the outcome studied.
#' @param mr_model model for MR. Possible values are "linear" or "logit". Default is "linear".
#' @param init the init value of theta for MCMC estimation. It can be a specific numeric or a string of "median", "egger" and "ivw", which means the initial value of the iteration will be calculated automatically by the above method. Default is "median".
#' @param n.iter an integer standing for the number of iterations. Default is 500.
#'
#' @return a list containing:
#'   thetaList: a vector cantaining the result of MCMC sampling of the causal parameter we want to estimate.
#'   mean: the mean estimate of the causal parameter.
#'   se: the SE of the estimation.
#'   sd: the SD of the estimation.
#'   lower: the lower boundary of the 95% CI of the causal estimation.
#'   upper: the upper boundary of the 95% CI of the causal estimation.
#'   Rhat: a indicator to measure the convergence (at convergence, Rhat <= 1.05).
#' @export
#'
#' @examples
#'
#'
mr <- function(df,selectsnp,exposureName,outcomeName,mr_model="linear",init="median",n.iter=500){
  library(rstan)
  library(MendelianRandomization)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)

  stanmodelcode <-'
/* lg_t.stan */
functions {
// Vector square root
vector vecsqrt(vector x) {
	vector[dims(x)[1]] res;
	for (m in 1:dims(x)[1]){
		res[m] = sqrt(x[m]);
	}
return res; }
}
data {
  int<lower=0> N;
  int<lower=0> J;
  matrix[N,J] Z;
  vector[N] X;
  vector[N] Y;
}
parameters {
  real <lower=0> sigmax;
  real <lower=0> sigmay;
  real <lower=0> sigmaalpha;
  real <lower=0> r1_global;
  real <lower=0> r2_global;
  real mualpha;
  real omegax;
  real omegay;
  real deltax;
  real deltay;
  real theta;
  vector[N] u;
  vector[J] z;
  vector<lower=0>[J] r1_local;
  vector<lower=0>[J] r2_local;
  vector[J] alpha;
}
transformed parameters {
   real<lower=0> tau;
   vector<lower=0> [J] lambda;
   vector[J] beta;
   tau      = r1_global * sqrt(r2_global);
   lambda	  = r1_local .* vecsqrt(r2_local);
   beta	    =  z .* lambda*tau;
}
model {
  X 	~ normal(omegax+Z*alpha+u*deltax, sigmax);
  Y   ~ normal(omegay+Z*beta+X*theta+u*deltay, sigmay);
  u 	~ normal(0,1);

  for(k in 1:J){
    alpha[k] ~ normal(mualpha, sigmaalpha);
  }
// Constructing the prior for the lambda vector
    z ~ normal(0, 1);
    r1_local ~ normal(0.0, 1.0);
    r2_local ~ inv_gamma(0.5, 0.5);
// Constructing the prior for tau
    r1_global ~ normal(0.0, 1.0);
    r2_global ~ inv_gamma(0.5, 0.5);
    }
'

  stanmodelcodeLogit <-'
/* lg_t.stan */
functions {
// Vector square root
vector vecsqrt(vector x) {
	vector[dims(x)[1]] res;
	for (m in 1:dims(x)[1]){
		res[m] = sqrt(x[m]);
	}
return res; }
}
data {
  int<lower=0> N;
  int<lower=0> J;
  matrix[N,J] Z;
  vector[N] X;
  int Y[N];
}
parameters {
  real <lower=0> sigmax;
  real <lower=0> sigmay;
  real <lower=0> sigmaalpha;
  real <lower=0> r1_global;
  real <lower=0> r2_global;
  real mualpha;
  real omegax;
  real omegay;
  real deltax;
  real deltay;
  real theta;
  vector[N] u;
  vector[J] z;
  vector<lower=0>[J] r1_local;
  vector<lower=0>[J] r2_local;
  vector[J] alpha;
}
transformed parameters {
   real<lower=0> tau;
   vector<lower=0> [J] lambda;
   vector[J] beta;
   vector<lower=0>[N] odds;
   vector<lower=0, upper=1>[N] prob;
   tau      = r1_global * sqrt(r2_global);
   lambda	  = r1_local .* vecsqrt(r2_local);
   beta	    =  z .* lambda*tau;
   for (i in 1:N){
       odds[i] = exp(omegay+Z[i]*beta+X[i]*theta+u[i]*deltay);
       prob[i] = odds[i] / (odds[i] + 1);
   }
}
model {
  X 	~ normal(omegax+Z*alpha+u*deltax, sigmax);
  Y   ~ bernoulli(prob);
  u 	~ normal(0,1);
  for(k in 1:J){
    alpha[k] ~ normal(mualpha, sigmaalpha);
  }
// Constructing the prior for the lambda vector
    z ~ normal(0, 1);
    r1_local ~ normal(0.0, 1.0);
    r2_local ~ inv_gamma(0.5, 0.5);
// Constructing the prior for tau
    r1_global ~ normal(0.0, 1.0);
    r2_global ~ inv_gamma(0.5, 0.5);
    }
'

  exposure <- df[,exposureName]
  outcome <- df[,outcomeName]
  s <- selectsnp
  N <- nrow(df)
  J <- length(s)
  X <- array(exposure,dim=N)
  Y <- array(outcome,dim=N)
  Z <- as.matrix(df[,s],dim=c(N,J))
  mydata <- list(N=N,J=J,X=X,Y=Y,Z=Z)

  betaX <- array(NA, dim=J)
  betaY <- array(NA, dim=J)
  sebetaY <- array(NA, dim=J)
  sebetaX <- array(NA, dim=J)
  for(isnp in 1:J){
    regX <- lm(X ~ Z[,isnp])
    regY <- lm(Y ~ Z[,isnp])
    betaX[isnp] <- summary(regX)$coefficients[2,1]
    sebetaX[isnp] <- summary(regX)$coefficients[2,2]
    betaY[isnp] <- summary(regY)$coefficients[2,1]
    sebetaY[isnp] <- summary(regY)$coefficients[2,2]
  }

  oggetto = mr_input(bx = as.numeric(betaX),
                     bxse = as.numeric(sebetaX),
                     by = as.numeric(betaY),
                     byse = as.numeric(sebetaY),
                     correlation = cor(Z),
                     exposure = "X ", outcome = "Y",
                     snps = colnames(Z))
  if(init=="median"){
    if(J<3){
      return(message("median initial needs at least 3 IVs"))
    }
    risultato = mr_allmethods(oggetto, method = "median")
    thetamedianestimate = risultato$Values[3,2]
  }else if(init=="egger"){
    if(J<3){
      return(message("egger initial needs at least 3 IVs"))
    }
    risultato = mr_allmethods(oggetto, method = "egger")
    thetamedianestimate = risultato$Values[7,2]
  }else if(init=="ivw"){
    risultato = mr_allmethods(oggetto, method = "ivw")
    thetamedianestimate = risultato$Values[4,2]
  }else if(class(init)=="numeric"){
    thetamedianestimate = init[1]
  }else{
    return(message("no this method to get initial value!"))
  }

  init_list = list(c1=list(theta=thetamedianestimate,
                           beta=rep(0,J),alpha=betaX,deltax=0,
                           deltay=0,u=rep(0,N)))
  if(mr_model=="linear"){
    fit <- stan(model_code=stanmodelcode, init=init_list, iter=n.iter,
                chains=1, verbose=F,data=mydata)
  }else if(mr_model=="logit"){
    fit <- stan(model_code=stanmodelcodeLogit, init=init_list, iter=n.iter,
                chains=1, verbose=F,data=mydata)
  }else{
    return(message("no this MR model!"))
  }

  theta = extract(fit,pars='theta',permuted=FALSE)
  motheta = monitor(theta,digits_summary=5)
  re = list(thetaList=theta,mean=motheta$mean,se=motheta$se_mean,sd=motheta$sd,
             lower=motheta$`2.5%`,upper=motheta$`97.5%`,Rhat=motheta$Rhat)
  return(re)
}
