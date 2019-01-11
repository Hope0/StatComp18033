#' @title MLE using R.
#' @description MLE is a fundamental estimation method of parameter, the function is working in condition that sample is  independent identical distribution (iid)
#' @importFrom stats optim
#' @param n the numbers of data
#' @param data which dim is n*p, p is the length of one sample
#' @param stheta the vector with length p, the initial value of optimization function
#' @param llm the function of likelihood of one sample and the parameter, input should be (parameter,sample)
#' @param is.log logical; if TRUE, the return values of llm are given as log
#' @return the MLE estimation of parameter
#' @examples
#' \dontrun{
#' llm1<- function(data,theta) return(dnorm(data,theta[1],theta[2]))
#' llm2<- function(data,theta) return(dnorm(data,theta[1],theta[2],log=TRUE))
#' MLE(10,X[1:10],c(-1,2),llm1,FALSE)
#' MLE(100,X,c(2,0.5),llm2,TRUE)
#' }
#' @export
MLE <- function(n,data,stheta,llm,is.log=TRUE){
  data <- as.matrix(data)
  LM<-function(theta){
    l <- 0
    for(i in 1:n) l=l+log(llm(data[i,],theta))
    return(-l)
  }
  LLM <-function(theta){
    l <- 0
    for(i in 1:n) l=l+llm(data[i,],theta)
    return(-l)
  }
  if(is.log==FALSE){
    ans <- optim(stheta,LM)
    return(ans$par)
  }
  if(is.log==TRUE){
    ans <- optim(stheta,LLM)
    return(ans$par)
  }
  return("error")
}
