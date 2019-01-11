#' @title create a MCMC chain  using R.
#' @description  Markov chain Monte Carlo (MCMC) methods comprise a class of algorithms for sampling from a probability distribution.
#' @importFrom stats rnorm
#' @importFrom stats dnorm
#' @importFrom stats dexp
#' @importFrom stats rexp
#' @importFrom stats runif
#' @param N the length of the Markov chain
#' @param ntheta the length of parameter
#' @param datas the data and the prior parameter of Bayesian, also the fixed parameters of LM
#' @param LM  the likelihood function, the input should be a sample of markov chain, and data.
#' @param ssigma the fixed parameters of transition function, such as variance of normal transition function, the acquiescent is 1
#' @param xstart the initial value of the chain, the length is equal to ntheta, the acquiescent is 0
#' @param rg the function of transition function. The input should be a sample of markov chain and ssigma, the output is another, one example is one dimensional normal-based on rnorm
#' @param tdf the function of density function or probability of transition function. The input should be a sample of markov chain, ssigma, the output of function rg, one example is one dimensional normal-based on dnorm.
#' @param is.log logical; if TRUE, the return values of LM and tdf p are given as log
#' @return chain:a markov chain and accept:acceptance times
#' @examples
#' \dontrun{
#' X <- rnorm(10,1,2)
#' LMn1 <- function(theta,data){
#'   t <- 1
#'   for (i in 1:length(data)) t=t*dnorm(data[i],theta[1],2)
#'   return(t)
#' }
#' rnorm1 <- function(mean,sigma) return(rnorm(1,mean,sigma))
#' dnorm1 <- function(mchain,r,ssigma) return(dnorm(mchain,mean=r,sd=ssigma))
#' cha <- MCMCchain(1000,1,X,LMn1,ssigma=1,xstart=0,rg=rnorm1,tdf=dnorm1)
#'
#' LMn2 <- function(theta,data){
#'   t <- 0
#'   for (i in 1:length(data)) t=t+dnorm(data[i],theta[1],theta[2],log=T)
#'   return(t)
#' }
#' t1 <- function(i,ssigma){
#'   j=i
#'   j[1] <- rnorm(1,i[1],ssigma[1])
#'   j[2] <- rexp(1,1/i[2])
#'   return(j)
#' }
#' d1<- function(i,j,ssigma){
#'   p1 <- dnorm(j[1],i[1],ssigma[1])
#'   p2 <- dexp(j[2],1/i[2])
#'   return(log(p1*p2))
#' }
#' cha2 <- MCMCchain(5000,2,X,LMn2,ssigma=c(1,2),xstart=c(0,1),rg=t1,tdf=d1,is.log=TRUE)
#' mean(cha2$c[1,][1000:5000])
#' mean(X)
#' mean(cha2$c[2,][1000:5000])
#' var(X)
#' }
#' @export
MCMCchain <- function(N,ntheta,datas,LM,ssigma=1,xstart=0,rg,tdf,is.log=FALSE){
  k=0
  mchain <- matrix(0,ntheta,N)
  mchain[,1] <- xstart
  for(i in 1:(N-1)){
    r <-rg(mchain[,i],ssigma)
    if(is.log==FALSE)a = LM(r,datas) * tdf(mchain[,i],r,ssigma) / LM(mchain[,i],datas) / tdf(r,mchain[,i],ssigma)
    if(is.log==TRUE) a = exp(LM(r,datas) + tdf(mchain[,i],r,ssigma) - LM(mchain[,i],datas) - tdf(r,mchain[,i],ssigma))
    nr <-runif(1)
    if(nr<a) k=k+1
    if(nr<a) mchain[,i+1]=r else mchain[,i+1]=mchain[,i]
  }
  return(list(chain=mchain,accept=k))
}

