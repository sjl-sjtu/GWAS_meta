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
