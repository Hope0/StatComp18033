## ------------------------------------------------------------------------
set.seed(1)
x <- rnorm(100,0,1)
x
hist(x)


## ------------------------------------------------------------------------
#9.6
#It is obvious that theta belongs to [0,1]
#i chooose the beta(1+a,2-a) as the proposal distrubution
#
#(125,18,20,34)
#
CMD <- function(theta){
  #return the probabilities vector of the corresponding multinomial distribution 
  if(theta >1) return(c(0,0,0,0))
  if(theta <0) return(c(0,0,0,0))
  return(c(0.5+theta/4,(1-theta)/4,(1-theta)/4,theta/4))
}

LMCMD <- function(r,data){
  p <- CMD(r)
  return(p[1]^data[1]*p[2]^data[2]*p[3]^data[3]*p[4]^data[4])
}

MCMC2 <- function(N,data,startx=0.5){
  mchain <- numeric(N) 
  mchain[1] = startx
  k=0
  for(i in 1:(N-1)){
    r <- rbeta(1,1+mchain[i],2-mchain[i])
    a = LMCMD(r,data) * dbeta(mchain[i],1+r,2-r)/ LMCMD(mchain[i],data)/  dbeta(r,1+mchain[i],2-mchain[i])
    nr <-runif(1)
    if(nr<a) mchain[i+1]=r else mchain[i+1]=mchain[i]
  }
  return(mchain)
}


#cov

Gelman.Rubin <- function(psi) { 
  # psi[i,j] is the statistic psi(X[i,1:j]) 
  # for chain in i-th row of X 
  psi <- as.matrix(psi)
  n <- ncol(psi) 
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means 
  B <- n * var(psi.means) #between variance est. 
  psi.w <- apply(psi, 1, "var") #within variances 
  W <- mean(psi.w) #within est. 
  v.hat <- W*(n-1)/n + (B/n) #upper variance est. 
  r.hat <- v.hat / W #G-R statistic 
  return(r.hat) 
}

set.seed(123)
k=3
N=5000
data <- c(125,18,20,34)
startx <- c(0.1,0.5,0.9)
startn <-1000
a <- matrix(0,k,N)
for (j in 1:k){
  mchain2 <- MCMC2(N,data,startx[j])
  a[j,] <- mchain2
}

Gelman.Rubin(a)

n=5
psi <- t(apply(a, 1, cumsum)) 
for (i in 1:nrow(psi)) psi[i,] <- psi[i,] / (1:ncol(psi)) 
print(Gelman.Rubin(psi))
#plot psi for the three chains 
par(mfrow=c(2,2),mar=c(1,1,1,1)) 
for (i in 1:k) plot(psi[i, (n+1):N], type="l", xlab=i, ylab=bquote(psi)) 
par(mfrow=c(1,1)) #restore default
#plot the sequence of R-hat statistics rhat <- rep(0, n) for (j in (b+1):n) rhat[j] <- Gelman.Rubin(psi[,1:j]) plot(rhat[(b+1):n], type="l", xlab="", ylab="R") abline(h=1.1, lty=2)

b <- numeric(N-n+1)
for(i in n:N)
  b[i-n+1]=Gelman.Rubin(a[,1:i])
plot(n:N,b,type="l",ylab="R")
abline(h=1.2,col=3)
abline(h=1.1,col=2)

#top 50 5-50
plot(n:N/100,b,type="l",ylab="R")
abline(h=1.2,col=3)
abline(h=1.1,col=2)
# chain has approximately converged to the target distribution within approximately 5000 iterations 
#
# chain has approximately converged to the target distribution within approximately 5000 iterations 

## ------------------------------------------------------------------------

Sk_1 <- function(a,k){
  q <- sqrt(a^2*(k-1)/(k-a^2)) 
  return (1-pt(q,df=k-1))
}
Sk <- function(a,k){
  q <- sqrt(a^2*k/(k+1-a^2)) 
  return (1-pt(q,df=k))
}
difSK <- function(x,k) { 
  Sk_1(x,k)-Sk(x,k)
}
kset <- c(4:25,100,500,1000)
out <- 1:length(kset)
for (i in 1:length(kset)){
  out[i] <- uniroot( difSK
                     , lower = 0+1e-5, upper = sqrt(kset[i])-1e-5,k=kset[i]) $root
}
out
###########
#not all right 
###########

kset[ abs(out-sqrt(kset)) < sqrt(kset)*0.01]
#It is shown that when k large than 22 ,the root is a wrong ,so we musst change the  
n <- 1:length(kset)
Kwrongnum <- n[abs(out-sqrt(kset)) < sqrt(kset)*0.01]

#based on the curve and the increasing of the answer,I change the upper 
#Example : k=23
k=23
xx <- seq(0.01,sqrt(k)-1e-5,length=1000)
y <- difSK(xx,k)
plot(xx,y,type="l")

#Example : k=1000
k=1000
xx <- seq(0.01,sqrt(k)-1e-5,length=1000)
y <- difSK(xx,k)
plot(xx,y,type="l")

#change upper to 3

for (i in Kwrongnum){
  out[i] <- uniroot( difSK
                     , lower = 0+1e-5, upper =3,k=kset[i]) $root
}
names(out) <- kset

out


## ------------------------------------------------------------------------
iris

## ----echo=FALSE----------------------------------------------------------
names(iris) <- c("SL", "SW", "PL", "PW", "SPP")
mod.iris <- lm(cbind(SL, SW, PL, PW) ~ SPP, data=iris)
e = resid(mod.iris)
y  = fitted(mod.iris)
par(mfrow=c(2,2),mar=c(1,1,1,1)) 
qqnorm(e[,1], pch = 20, main="Normal Probability Plot")
qqline(e[,1])
qqnorm(e[,2], pch = 20, main="Normal Probability Plot")
qqline(e[,2])
qqnorm(e[,3], pch = 20, main="Normal Probability Plot")
qqline(e[,3])
qqnorm(e[,4], pch = 20, main="Normal Probability Plot")
qqline(e[,4])
plot(y[,1],e[,1],col = 2,pch = 20)
plot(y[,2],e[,2],col = 3,pch = 20)
plot(y[,3],e[,3],col = 4,pch = 20)
plot(y[,4],e[,4],col = 5,pch = 20)

## ----echo=FALSE----------------------------------------------------------
X <- c(0:4)
p <- c(0.1,0.2,0.2,0.2,0.3)
csp <-  cumsum(p)
#the inverse transform method
ITMSample <- function(x,csp,N){
  #x???sample csp???Cumulative Distribution Function  N???Number of Sample
  csp<-c(csp,1)
  y<-1:N
  #note theresult
  r<-runif(N)
  #gemerate random numbers
  for(i in 1:N){
    k = r[i]
    j=1 
    while(csp[j]<k) j=j+1
    y[i]=x[j]
  }
  return(y) 
}
set.seed(123)
N<-1000
rus1<-ITMSample(X,csp,N)
tableRE1 <- rbind(table(rus1),p)
tableRE1[1,]<-tableRE1[1,]/N
rownames(tableRE1)<- c('empirical','theoretical')
knitr::kable(tableRE1, caption="Table 1???a relative frequency table using the inverse transform method to generate ")
#make the table

#sample
rus2<-sample(X,N,replace = TRUE,prob=p)
tableRE2 <- rbind(table(rus2),p)
tableRE2[1,]<-tableRE2[1,]/N
rownames(tableRE2)<- c('empirical','theoretical')
knitr::kable(tableRE2, caption="Table 2???a relative frequency table using Sample ")
#make the table


## ----echo=TRUE-----------------------------------------------------------
X <- c(0:4)
p <- c(0.1,0.2,0.2,0.2,0.3)
csp <-  cumsum(p)
#the inverse transform method
ITMSample <- function(x,csp,N){
  #x???sample csp???Cumulative Distribution Function  N???Number of Sample
  csp<-c(csp,1)
  y<-1:N
  #note theresult
  r<-runif(N)
  #gemerate random numbers
  for(i in 1:N){
    k = r[i]
    j=1 
    while(csp[j]<k) j=j+1
    y[i]=x[j]
  }
  return(y) 
}
set.seed(123)
N<-1000
rus1<-ITMSample(X,csp,N)
tableRE1 <- rbind(table(rus1),p)
tableRE1[1,]<-tableRE1[1,]/N
rownames(tableRE1)<- c('empirical','theoretical')
#knitr::kable(tableRE1, caption="Table 1???a relative frequency table using the inverse transform method to generate ")
#make the table

#sample
rus2<-sample(X,N,replace = TRUE,prob=p)
tableRE2 <- rbind(table(rus2),p)
tableRE2[1,]<-tableRE2[1,]/N
rownames(tableRE2)<- c('empirical','theoretical')
#knitr::kable(tableRE2, caption="Table 2???a relative frequency table using Sample ")
#make the table

## ----echo=FALSE----------------------------------------------------------
#find the maximum of density function 
x<-seq(0,1,0.0001)
N<-length(x)
y<-1:N
for(i in 1:N){
  y[i]<-dbeta(x[i],3,2)
}
#function acceptance-rejection method
ARMsample_beta <- function(beta1,beta2,maxd,N){
  #beta1,beta2 the non-negative parameters of the Beta distribution.
  #maxd the c of acceptance-rejection method
  #N number before rejection
  r1<-runif(N)
  #??????????????????
  y <- NULL
  for(i in 1:N){
    if(dbeta(r1[i],beta1,beta2)> maxd*runif(1) ) y<-c(y,r1[i])
    #????????????
  }
  return(y)
}
set.seed(1234)
rus<-ARMsample_beta(3,2,1.8,2000)
#length(rus)
#???????????????37%
rus1<-rus[1:1000]
hist(rus1,main=" histogram of the sample and the theoretical density superimposed",freq=FALSE)
lines(x,y)

## ----echo=TRUE-----------------------------------------------------------
#find the maximum of density function 
x<-seq(0,1,0.0001)
N<-length(x)
y<-1:N
for(i in 1:N){
  y[i]<-dbeta(x[i],3,2)
}
#max(y)
#function acceptance-rejection method
ARMsample_beta <- function(beta1,beta2,maxd,N){
  #beta1,beta2 the non-negative parameters of the Beta distribution.
  #maxd the c of acceptance-rejection method
  #N number before rejection
  r1<-runif(N)
  #??????????????????
  y <- NULL
  for(i in 1:N){
    if(dbeta(r1[i],beta1,beta2)> maxd*runif(1) ) y<-c(y,r1[i])
    #????????????
  }
  return(y)
}
set.seed(1234)
rus<-ARMsample_beta(3,2,1.8,2000)
#length(rus)
#???????????????37%
rus1<-rus[1:1000]
#hist(rus1,main=" histogram of the sample and the theoretical density superimposed",freq=FALSE)
#lines(x,y)

## ------------------------------------------------------------------------
rExp_Gamma <- function(r,beta,N){
  #N number
  y<-1:N
  for (i in 1:N){
    x<- rgamma(1,r,beta)
    y[i]<-rexp(1,x)
  }
  return(y)
}
set.seed(12)
rus<-rExp_Gamma(4,2,1000)
hist(rus,breaks=seq(0,max(rus)+0.5,0.5),main="histogram of the random number")

## ------------------------------------------------------------------------
betacdfMC <- function(x,beta1,beta2,N=10000) {
  #x F(x),x can be a number or a vector
  ##beta1,beta2 the non-negative parameters of the Beta distribution.
  #N numbers of random number
  y<-rbeta(N,3,3)
  X<-matrix(x,length(x),N)
  Y<-matrix(y,length(x),N,byrow=T)
  #
  return(rowSums(Y<=X)/N)
}
set.seed(123)
A1_1<-betacdfMC(seq(0.1,0.9,0.1),3,3)
#A1_1
#100000
A1_2<-betacdfMC(seq(0.1,0.9,0.1),3,3,100000)
#A1_2
#1000000
A1_3<-betacdfMC(seq(0.1,0.9,0.1),3,3,1000000)
#A1_3
#theoretical value
A1_t<-pbeta(seq(0.1,0.9,0.1),3,3)
#A1_t

A1 <- rbind(A1_1,A1_2,A1_3,A1_t)
colnames(A1)<-c(seq(0.1,0.9,0.1))
rownames(A1)<-c(1e+04,1e+05,1e+06,'theoretical')
A1


## ------------------------------------------------------------------------
invcdfRay <- function(u,e) return(sqrt(-2*e*e*log(1-u)))
#e :??
simpleicRa <- function(Sigma,N=1000){
#sigma: the parameter of Rayleigh(??) distribution
#N numbers of random 
#The function creates N random numbers And N antithetic random number
#the first N number is indenpendent ,the last N is antithetic of the top N
u<-runif(N)
y<-rep(0,2*N)
for(i in 1:N){
y[i]<-invcdfRay(u[i],Sigma)
y[i+N]<-invcdfRay(1-u[i],Sigma)
}
return(y)
}
#I create 1000(N) random numbers 1000(n) times
set.seed(123)
n=100
N=10000
sigma=1
x=xnew<-matrix(0,n,2*N)
for(i in 1:n){
x[i,]<-simpleicRa(sigma,N)
}
for(i in 1:n){
xnew[i,]<-simpleicRa(sigma,N)
}
n=n-n%%2
#xin and x1: (X1+X2)/2 and y1 : n estimated value  
xin<-matrix(0,n,N)
for(i in 1:n){
xin[i,] <- (x[i,1:N]+xnew[i,1:N])/2
}
y1<-rowMeans(xin)

#x2 and xan X+X' and y2 : n estimated value 
xanti<-matrix(0,n,N)
for(i in 1:n){
  xanti[i,] <- (x[i,1:N]+x[i,N+1:N])/2
}
y2<-rowMeans(xanti)
var(y1)
var(y2)
var(y2)/var(y1)
#var(y2) < var (y1)


## ----echo=FALSE----------------------------------------------------------
funcg <- function(x) if(x>=1) x*x*exp(-x*x/2)/sqrt(2*pi) else 0
funcg_vec <-function(X) {
  n<-length(X)
  y<-1:n
  for(i in 1:n)
    y[i]<-funcg(X[i])
  return(y)
}

funcf1 <- function(x) if(x>=1) x*exp(-x*x/2)*exp(1/2) else 0
funcf1_vec <-function(X) {
  n<-length(X)
  y<-1:n
  for(i in 1:n)
    y[i]<-funcf1(X[i])
  return(y)
}
funcf2 <- function(x) if(x>=1) exp(-x/2)/2*exp(1/2) else 0
funcf2_vec <-function(X) {
  n<-length(X)
  y<-1:n
  for(i in 1:n)
    y[i]<-funcf2(X[i])
  return(y)
}
x<-seq(1,10,0.001)
y<-funcg_vec(x)
plot(x,y,ylim=c(0,1),type='l')
y2<-funcf2_vec(x)
lines(x,y2,lty=2)
text(2,0.9,'f2',cex=2)
y1<-funcf1_vec(x)
lines(x,y1,lty=10)
text(6,0.1,'f1',cex=2)
text(1.5,0.2,'g(x)',cex=2)
plot(x,y1/y,type='l',ylim=c(0,10))
lines(x,y2/y,lty=10)
text(6,2,'f1/g',cex=2)
text(4,8,'f2/g',cex=2)

## ------------------------------------------------------------------------
funcg <- function(x) if(x>=1) x*x*exp(-x*x/2)/sqrt(2*pi) else 0
funcg_vec <-function(X) {
  n<-length(X)
  y<-1:n
  for(i in 1:n)
    y[i]<-funcg(X[i])
  return(y)
}

funcf1 <- function(x) if(x>=1) x*exp(-x*x/2)*exp(1/2) else 0
funcf1_vec <-function(X) {
  n<-length(X)
  y<-1:n
  for(i in 1:n)
    y[i]<-funcf1(X[i])
  return(y)
}
funcf2 <- function(x) if(x>=1) exp(-x/2)/2*exp(1/2) else 0
funcf2_vec <-function(X) {
  n<-length(X)
  y<-1:n
  for(i in 1:n)
    y[i]<-funcf2(X[i])
  return(y)
}
#x<-seq(1,10,0.001)
#y<-funcg_vec(x)
#plot(x,y,ylim=c(0,1),type='l')
#y2<-funcf2_vec(x)
#lines(x,y2,lty=2)
#text(2,0.9,'f2',cex=2)
#y1<-funcf1_vec(x)
#lines(x,y1,lty=10)
#text(6,0.1,'f1',cex=2)
#text(1.5,0.2,'g(x)',cex=2)
#plot(x,y1/y,type='l',ylim=c(0,10))
#lines(x,y2/y,lty=10)
#text(6,2,'f1/g',cex=2)
#text(4,8,'f2/g',cex=2)

## ------------------------------------------------------------------------
set.seed(123)
N=10000
u1<-runif(N)
cdff1 <- function(x) 1-exp(-x^2/2+1/2)
Invcdff1 <- function(u) sqrt(1-2*log(1-u))
Invcdff1_vec <- function(X) {
  n<-length(X)
  y<-1:n
  for(i in 1:n)
    y[i]<-Invcdff1(X[i])
  return(y)
}
x1 <- Invcdff1_vec(u1)
y1<-sum(funcg_vec(x1)/funcf1_vec(x1))/N
y1
x1_an <- c(Invcdff1_vec(u1),Invcdff1_vec(1-u1))
y1_an <- sum(funcg_vec(x1_an)/funcf1_vec(x1_an))/N/2
y1_an

u2<-runif(N)
cdff2 <- function(x) 1-exp(-x/2+1/2)
Invcdff2 <- function(u) (1-2*log(1-u))
Invcdff2_vec <- function(X) {
  n<-length(X)
  y<-1:n
  for(i in 1:n)
    y[i]<-Invcdff2(X[i])
  return(y)
}
x2 <- Invcdff2_vec(u2)
y2<-sum(funcg_vec(x2)/funcf2_vec(x2))/N
y2
x2_an <- c(Invcdff2_vec(u2),Invcdff2_vec(1-u2))
y2_an <- sum(funcg_vec(x2_an)/funcf2_vec(x2_an))/N/2
y2_an

## ----echo=FALSE----------------------------------------------------------
GiniratioV <- function(X,Xmean=NA){
  #vector X
  if(is.na(Xmean)) Xmean=mean(X)
  sX <- sort(X)
  G=0
  n=length(X)
  for(i in 1:n){
    G=G+(2*i-n-1)*sX[i]
  }
  G=G/n/n/Xmean
  return(G)
}
N=10000
n=1000
set.seed(1234)
# standard lognormal
# random numbers from  standard lognormal and calculate G
# 10000 random numbers and 1000 times

Gln <- replicate(n, expr = { 
  X<-rlnorm(N)
  GiniratioV(X)
})

## ----echo=FALSE----------------------------------------------------------
mean(Gln)

## ----echo=FALSE----------------------------------------------------------
median(Gln)

## ----echo=FALSE----------------------------------------------------------
quantile(Gln,seq(0.1,0.9,0.1))

## ----echo=FALSE----------------------------------------------------------
hist(Gln,freq=F,main="the density histograms of estimated G from standard lognormal")

## ----echo=FALSE----------------------------------------------------------

# uniform 
# random numbers from uniform and calculate G
# 10000 random numbers and 1000 times

Gun <- replicate(n, expr = { 
  X<-runif(N)
  GiniratioV(X)
})

## ----echo=FALSE----------------------------------------------------------
mean(Gun)

## ----echo=FALSE----------------------------------------------------------
median(Gun)

## ----echo=FALSE----------------------------------------------------------
quantile(Gun,seq(0.1,0.9,0.1))

## ----echo=FALSE----------------------------------------------------------
hist(Gun,freq=F,main="the density histograms of estimated G from unifrom")

## ----echo=FALSE----------------------------------------------------------
# Bernoulli(0.1). 
# random numbers from  Bernoulli(0.1).  and calculate G
# 10000 random numbers and 1000 times


GBn <- replicate(n, expr = { 
  X<-sample(c(0,1),N,replace=T,p=c(0.9,0.1))
  GiniratioV(X)
})

## ----echo=FALSE----------------------------------------------------------
mean(GBn)

## ----echo=FALSE----------------------------------------------------------
median(GBn)

## ----echo=FALSE----------------------------------------------------------
quantile(GBn,seq(0.1,0.9,0.1))

## ----echo=FALSE----------------------------------------------------------
hist(GBn,freq=F,main="the density histograms of estimated G from Bernoulli(0.1)")

## ------------------------------------------------------------------------
GiniratioV <- function(X,Xmean=NA){
  #vector X
  if(is.na(Xmean)) Xmean=mean(X)
  sX <- sort(X)
  G=0
  n=length(X)
  for(i in 1:n){
    G=G+(2*i-n-1)*sX[i]
  }
  G=G/n/n/Xmean
  return(G)
}
N=10000
n=1000
set.seed(1234)
# standard lognormal
# random numbers from  standard lognormal and calculate G
# 10000 random numbers and 1000 times

Gln <- replicate(n, expr = { 
  X<-rlnorm(N)
  GiniratioV(X)
})
mean(Gln)
median(Gln)
quantile(Gln,seq(0.1,0.9,0.1))
#hist(Gln,freq=F)

# uniform 
# random numbers from uniform and calculate G
# 10000 random numbers and 1000 times

Gun <- replicate(n, expr = { 
  X<-runif(N)
  GiniratioV(X)
})
mean(Gun)
median(Gun)
quantile(Gun,seq(0.1,0.9,0.1))
#hist(Gun,freq=F)

# Bernoulli(0.1). 
# random numbers from  Bernoulli(0.1).  and calculate G
# 10000 random numbers and 1000 times


GBN <- replicate(n, expr = { 
  X<-sample(c(0,1),N,replace=T,p=c(0.9,0.1))
  GiniratioV(X)
})
mean(GBn)
median(GBn)
quantile(GBn,seq(0.1,0.9,0.1))
#hist(GBn,freq=F)

## ------------------------------------------------------------------------
erfc <- function(s){
  return(2*pnorm(s,0,sqrt(1/2))-1)
} 

varCILN <- function(x,a){
  # a level
  # x data from log-normal
  n <- length(x)
  ns_2 <- var(log(x))*n
  return(c(ns_2/qchisq(1-a/2,n-1),ns_2/qchisq(a/2,n-1)))
} 

GCI <- function(x,a){
  # a level
  # x data from log-normal
  ci <- varCILN(x,a)
  g1 <- erfc(sqrt(ci[1])/2)
  g2 <- erfc(sqrt(ci[2])/2)
  return(c(g1,g2))
}

## ------------------------------------------------------------------------
set.seed(100)
N=100
n=10000
X <- rlnorm(N)
gci <- GCI(X,0.05)
gci
g <- replicate(n, expr = { 
  X <- rlnorm(N)
  erfc(sd(log(X))/2)
})
sum(gci[1]<g & gci[2]>g)/n

## ------------------------------------------------------------------------
a=0.05
ns_2<-N-1
ci <- c(ns_2/qchisq(1-a/2,N-1),ns_2/qchisq(a/2,N-1))
g1 <- erfc(sqrt(ci[1])/2)
g2 <- erfc(sqrt(ci[2])/2)
gci_2 <-c(g1,g2)
sum(gci_2[1]<g & gci_2[2]>g)/n

## ----echo=FALSE----------------------------------------------------------
#cor.test(x,y,method="pearson")  Pearson product moment correlation ??,
#cor.test(x,y,method="kendall")  Kendall???s coe???cient ??
#cor.test(x,y,method="spearman") Spearman???s rank correlation coe???cient ??s

cortestall <- function(x,y){
  t1<-cor.test(x,y,method="pearson")
  t2<-cor.test(x,y,method="kendall")
  t3<-cor.test(x,y,method="spearman") 
  Assc <- c(t1$p.value,t2$p.value,t3$p.value)
  return(Assc)
}
library(MASS)
set.seed(12)
N=100
n=10000
mu <- c(0,0)
p <- 0.4
Var <- 1
sigma <- matrix(c(1,p,p,1),2,2)
testn <- replicate(n, expr = { 
  X <- mvrnorm(N,mu,Var*sigma)
  cortestall(X[,1],X[,2])
})

## ----echo=FALSE----------------------------------------------------------
sum(testn[1,]<=0.05)/n

## ----echo=FALSE----------------------------------------------------------
sum(testn[2,]<=0.05)/n

## ----echo=FALSE----------------------------------------------------------
sum(testn[3,]<=0.05)/n

## ----echo=FALSE,warning=FALSE--------------------------------------------
#install.packages("DirichletReg")
library(DirichletReg)
N=100
n=10000
testn <- replicate(n, expr = { 
  X <- rdirichlet(N,c(3,1))
  Z <- rdirichlet(N,c(1,3))
  Y <- 0.18*X + 0.5*Z
  cortestall(X[,2],Y[,1])
})

## ----echo=FALSE----------------------------------------------------------
sum(testn[1,]<=0.05)/n

## ----echo=FALSE----------------------------------------------------------
sum(testn[2,]<=0.05)/n

## ----echo=FALSE----------------------------------------------------------
sum(testn[3,]<=0.05)/n

## ------------------------------------------------------------------------
#cor.test(x,y,method="pearson")  Pearson product moment correlation ??,
#cor.test(x,y,method="kendall")  Kendall???s coe???cient ??
#cor.test(x,y,method="spearman") Spearman???s rank correlation coe???cient ??s

cortestall <- function(x,y){
  t1<-cor.test(x,y,method="pearson")
  t2<-cor.test(x,y,method="kendall")
  t3<-cor.test(x,y,method="spearman") 
  Assc <- c(t1$p.value,t2$p.value,t3$p.value)
  return(Assc)
}

library(MASS)
set.seed(12)
N=100
n=10000
mu <- c(0,0)
p <- 0.4
Var <- 1
sigma <- matrix(c(1,p,p,1),2,2)
testn <- replicate(n, expr = { 
  X <- mvrnorm(N,mu,Var*sigma)
  cortestall(X[,1],X[,2])
})
sum(testn[1,]<=0.05)/n
sum(testn[2,]<=0.05)/n
sum(testn[3,]<=0.05)/n

#install.packages("DirichletReg")
library(DirichletReg)

N=100
n=10000

testn <- replicate(n, expr = { 
  X <- rdirichlet(N,c(3,1))
  Z <- rdirichlet(N,c(1,3))
  Y <- 0.18*X + 0.5*Z
  cortestall(X[,2],Y[,1])
})
sum(testn[1,]<=0.05)/n
sum(testn[2,]<=0.05)/n
sum(testn[3,]<=0.05)/n

## ------------------------------------------------------------------------

#install.packages("bootstrap")
library(bootstrap)
library(DAAG)
data(law)


#compute the jackknife replicates, leave-one-out estimates 
set.seed(4321)
n <- length(law[,1])
theta.jack <- numeric(n) 
for (i in 1:n) 
theta.jack[i] <- cor(law[-i,])[1,2]
bias <- (n - 1) * (mean(theta.jack) - cor(law)[1,2])
bias
mjack <- mean(theta.jack)
se.jack <- sqrt((n-1)/n*sum((theta.jack-mjack)^2))
se.jack
c(original=cor(law)[1,2],bias=bias,se=se.jack)


## ------------------------------------------------------------------------
#install.packages("boot")
library(boot)
data(aircondit)
set.seed(4321)
#set up the bootstrap 
B <- 2000 #number of replicates 
n <- nrow(aircondit) #sample size 
a <- 0.05 #level
meantime <- numeric(B) #storage for replicates
#bootstrap estimate of standard error of R 

meanairc <- mean(aircondit[,1]) #mean of the original observed sample
for (b in 1:B) { 
#randomly select the indices 
i <- sample(1:n, size = n, replace = TRUE) 
Bairc <- aircondit[i,1]
meantime[b] <- mean(Bairc)
}
mmean <- mean(meantime)  #mean of the bootstrap replicate
SEmeantime <- sqrt(sum((meantime-mmean)^2)/(B-1))

#the standard normal
the_standard_normalCI <- c(meanairc-qnorm(1-a/2)*SEmeantime,meanairc+qnorm(1-a/2)*SEmeantime)
names(the_standard_normalCI) <- c("2.5%","97.5%")

#basic
percentileCI <- quantile(meantime,probs=c(a/2,1-a/2)) 
basic_CI <- c(2*meanairc-percentileCI[2],2*meanairc-percentileCI[1])

#percentile
#percentileCI

#Bca
x <- aircondit[,1]
alpha <- c(a/2,1-a/2)
zalpha <- qnorm(alpha)
z0 <- qnorm(sum(meantime < mmean) / length(meantime))
mean.jack <- numeric(n)
for (i in 1:n) {
mean.jack[i] <- mean(x[-i]) 
} 
L <- mean(mean.jack) - mean.jack 
sa <- sum(L^3)/(6 * sum(L^2)^1.5)
adj.alpha <- pnorm(z0 + (z0+zalpha)/(1-sa*(z0+zalpha))) 
BCa_CI <- quantile(meantime, adj.alpha, type=6)


#the standard normal
the_standard_normalCI
#the standard normal
basic_CI
#percentile
percentileCI
#Bca
BCa_CI

#another 
x=aircondit$hours
meanx<-function(x,i) mean(x[i])
de <- boot(data=x,statistic=meanx, R =B)
ci <- boot.ci(de,type=c("norm","basic","perc","bca"))
ci

## ------------------------------------------------------------------------
attach(scor)
datascores  <-scor
n <- length(datascores[,1])
theta.jack <- numeric(n) 
calthta <- function(X){
#X show be a matrix n*n
eig <- eigen(X)
sort_eig <- sort(eig$values,decreasing=T)
return(sort_eig[1]/sum(sort_eig))
}

for (i in 1:n) theta.jack[i] <- calthta(cov(datascores[-i,]))
bias <- (n - 1) * (mean(theta.jack) - calthta(cov(datascores)))
bias
mjack <- mean(theta.jack)
sc.jack <- sqrt((n-1)/n*sum((theta.jack-mjack)^2))
sc.jack

## ------------------------------------------------------------------------
#install.packages("DAAG")
#install.packages("caTools")
library(DAAG)
library(caTools)
attach(ironslag)
n <- length(magnetic) #in DAAG ironslag 
K=2 # knife number
kjack <- combs(1:n,K) #create the all combinations of 2 elements from 1:n
N <- choose(n,K) #number of K
e1 <- e2 <- e3 <- e4 <- matrix(,N,K)

# for n-fold cross validation 
# fit models on leave-two-out samples 


for (i in 1:N) { 
k <- kjack[i,]
y <- magnetic[-k] 
x <- chemical[-k]
J1 <- lm(y ~ x) 
yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k] 
e1[i,] <- magnetic[k] - yhat1
J2 <- lm(y ~ x + I(x^2)) 
yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] + J2$coef[3] * chemical[k]^2 
e2[i,] <- magnetic[k] - yhat2
J3 <- lm(log(y) ~ x) 
logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k] 
yhat3 <- exp(logyhat3) 
e3[i,] <- magnetic[k] - yhat3
J4 <- lm(log(y) ~ log(x)) 
logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[k]) 
yhat4 <- exp(logyhat4) 
e4[i,] <- magnetic[k] - yhat4
}

c(sum(e1^2),sum(e2^2),sum(e3^2),sum(e4^2))/N/K

## ------------------------------------------------------------------------
#install.packages("RANN")
#install.packages("energy")
#install.packages("devtools")
#install.packages("Ball")
library(RANN)
library(boot) 
library(energy) 
library(Ball) 
library(devtools) 


#8.1
attach(chickwts) 
set.seed(1234)
x <- sort(as.vector(weight[feed == "soybean"])) 
y <- sort(as.vector(weight[feed == "linseed"])) 
cvMTS <- function(X,Y){
#Cram??r???von Mises Two Sample  statistic
#X  Y : vector
XYrank <- rank(c(X,Y))
N <- length(X)
M <- length(Y)
U= N*sum((XYrank[1:N]-1:N)^2)+M*sum((XYrank[(N+1):(N+M)]-1:M)^2)
T=U/M/N/(N+M)-(4*M*N-1)/6/(M+N)
return(T)
}
## start the replicate
R <- 999 #number of replicates 
z <- c(x, y) #pooled sample 
K <- 1:length(z)
T <- numeric(R) #storage for replicates options(warn = -1) 
T0 <- cvMTS(x,y)
for (i in 1:R) { 
#generate indices k for the first sample 
k <- sample(K, size = 14, replace = FALSE) 
x1 <- z[k] 
y1 <- z[-k]
T[i] <- cvMTS(x1,y1)
}
P0 <- mean(c(T0, T) >= T0)   
pvalue <- min(2*P0,2*(1-P0))
pvalue
#pvalue < 0.05
#Thus, none of the replicates are as large as the observed test statistic. 
#Here the sample evidence supports the alternative hypothesis that the distributions di???er. 
detach(chickwts)

## ------------------------------------------------------------------------
#PLUS

#NN
Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]
n2 <- sizes[2]
n <- n1 + n2 
if(is.vector(z))  z <- data.frame(z,0)
z <- z[ix, ]
NN <- nn2(data=z, k=k+1) 
block1 <- NN$nn.idx[1:n1,-1] 
block2 <- NN$nn.idx[(n1+1):n,-1] 
i1 <- sum(block1 < n1 + .5)
i2 <- sum(block2 > n1+.5) 
return((i1 + i2) / (k * n) )
}

NNtest_1 <- function(x,y,NumR=999,k=3){
N <- c(length(x), length(y)) 
z <- c(x,y)
boot.obj <- boot(data = z, statistic = Tn, R = NumR, sim = "permutation", sizes = N,k=k)
ts <- c(boot.obj$t0,boot.obj$t) 
p.value <- mean(ts>=ts[1])
return(p.value)
}

#Energy test 
Energy_test <- function(x,y,NumR=999){
N <- c(length(x), length(y)) 
z <- c(x,y)
boot.obs <- eqdist.etest(z, sizes=N, R=NumR) 
p.value <- boot.obs$p.value
return(p.value)
}


#install_github("Mamba413/Ball", build_vignettes = TRUE)

#
alpha <- 0.05; 

##1
##Unequal variances and equal expectations
set.seed(123)
m <- 1000; #number od tests
k<-3; 
p<-2; 
mu <- 0; 
n1 <- n2 <- 20; 
R<-999;
sigma1 <- 1 
sigma2 <- 2.5
#x N(0,1) 20 y N(0,4) 20
p.values1 <- matrix(NA,m,3) 
for(i in 1:m){ 
x <- rnorm(n1,mu,sigma1)
y <- rnorm(n2,mu,sigma2)
z <- c(x,y) 
p.values1[i,1] <- NNtest_1(x,y,R,k) 
p.values1[i,2] <- Energy_test(x,y,R)
p.values1[i,3] <- bd.test(x,y,R=R,seed=i*123)$p.value 
}
pow1 <- colMeans(p.values1<alpha)
pow1

##2
##Unequal variances and unequal expectations
set.seed(123)
m <- 1000; #number od tests
k<-3; 
p<-2; 
mu1 <- 0; 
mu2 <- 1
n1 <- n2 <- 20; 
R<-999;
sigma1 <-1
sigma2 <-2
R<-100;
#x N(0,1) 20  y N(1,4) 20
p.values2 <- matrix(NA,m,3) 
for(i in 1:m){ 
x <- rnorm(n1,mu1,sigma1)
y <- rnorm(n2,mu2,sigma2)
z <- c(x,y) 
p.values2[i,1] <- NNtest_1(x,y,R,k) 
p.values2[i,2] <- Energy_test(x,y,R)
p.values2[i,3] <- bd.test(x,y,R=R,seed=i*123)$p.value 
}
pow2 <- colMeans(p.values2<alpha)
pow2

##3.1
##Non-normal distributions: t distribution with 1 df (heavy-tailed distribution),bimodel distribution (mixture of two normal distributions)
set.seed(123)
m <- 1000; #number od tests
k<-3; 
n1 <- n2 <- 20; 
df=1
R<-999;
p=0.4
mu1 =-1
mu2 = 1
sigma1 =1
sigma2=1
#x t1 30  y 0.3N(-1,1)+0.7N(1,2) 30
p.values3 <- matrix(NA,m,3) 
for(i in 1:m){ 
x <- rt(n1,df=df)
y <- 0.3 * rnorm(n2,mu1,sigma1) + 0.7*rnorm(n2,mu2,sigma2)
z <- c(x,y) 
p.values3[i,1] <- NNtest_1(x,y,R,k) 
p.values3[i,2] <- Energy_test(x,y,R)
p.values3[i,3] <- bd.test(x,y,R=999,seed=i*123)$p.value 
}
pow3 <- colMeans(p.values3<alpha)
pow3

#4
#Unbalanced samples (say, 1 case versus 10 controls)
set.seed(123)
m <- 1000; #number od tests
k<-3; 
n1 <- 100; 
n2 <- 100/10
R<-999;
mu1 =-1
mu2 = 0
sigma1 =1
sigma2=2
#x N(-1,1) 100  y N(1,2) 10
p.values4 <- matrix(NA,m,3) 
for(i in 1:m){ 
x <- rnorm(n1,mu1,sigma1)
y <- rnorm(n2,mu2,sigma2)
z <- c(x,y) 
p.values4[i,1] <- NNtest_1(x,y,R,k) 
p.values4[i,2] <- Energy_test(x,y,R)
p.values4[i,3] <- bd.test(x,y,R=999,seed=i*123)$p.value 
}
pow4 <- colMeans(p.values4<alpha)
pow4

## ------------------------------------------------------------------------
#9.3
#the proposal distributio:normal
theta=1
eta=0
MCMC <- function(N,theta,eta,ssigma=1,xstart=0,LM){
#  k=0
mchain <- numeric(N) 
mchain[1] <- xstart
for(i in 1:(N-1)){
r <-rnorm(1,mchain[i],ssigma)
a = LM(r,theta,eta) * dnorm(mchain[i],mean=r,sd=ssigma) / LM(mchain[i],theta,eta) / dnorm(r,mean=mchain[i],sd=ssigma)
nr <-runif(1)
#    if(nr<a) k=k+1
if(nr<a) mchain[i+1]=r else mchain[i+1]=mchain[i]
}
#  return(c(k,mchain))
return(mchain)
}
set.seed(123)
sigma = 8
#sigma of the proposal distributio
N=5000
n=1000

LMcauchy <- function(x,theta,eta){
return(1/theta/pi/(1+((x-eta)/theta)^2))
}
theta=1
eta=0

MCchaincauchy <- MCMC(N,1,0,sigma,0,LMcauchy)
#accrptation rate : 30%
#MCchaincauchy[1]
hist(MCchaincauchy[(n+1):N],freq = FALSE,main="MC chain of cauchy")
x <- seq(min(MCchaincauchy[(n+1):N]),max(MCchaincauchy[(n+1):N]),length=1000)
lines(x,LMcauchy(x,1,0))
xx <- seq(0.1,0.9,0.1)  
MC_deciles <- quantile(MCchaincauchy[(n+1):N],xx)
Th_deciles <- qcauchy(xx)
ans <- rbind(MC_deciles,Th_deciles,abs(Th_deciles-MC_deciles))
rownames(ans)[3] <- "differ"
ans
qqplot(qcauchy(ppoints(N-n)),MCchaincauchy[(n+1):N])
qqline(MCchaincauchy[(n+1):N],distribution = qcauchy)

## ------------------------------------------------------------------------
#9.6
#It is obvious that theta belongs to [0,1]
#i chooose the uniform as the prior distrubution
#
#(125,18,20,34)
#

set.seed(214)
CMD <- function(theta){
#return the probabilities vector of the corresponding multinomial distribution 
if(theta >1) return(c(0,0,0,0))
if(theta <0) return(c(0,0,0,0))
return(c(0.5+theta/4,(1-theta)/4,(1-theta)/4,theta/4))
}

LMCMD <- function(r,data){
p <- CMD(r)
return(p[1]^data[1]*p[2]^data[2]*p[3]^data[3]*p[4]^data[4])
}

MCMC2 <- function(N,data,startx=0.5){
mchain <- numeric(N) 
mchain[1] = startx
#k=0
for(i in 1:(N-1)){
r <- runif(1)
a = LMCMD(r,data)  / LMCMD(mchain[i],data)
nr <-runif(1)
if(nr<a) mchain[i+1]=r else mchain[i+1]=mchain[i]
}
return(mchain)
}

data <- c(125,18,20,34)
mchain2 <- MCMC2(N,data,0.5)
# PS: acceptation rate is 16%
achain <- mchain2[1001:N]
hist(achain)
mean(achain)

## ------------------------------------------------------------------------
LMF <-function(y,theta,eta){
  #theta : scale parameter
  #eta : the location parameter
  1/(theta*3.141592653*(1+((y-eta)/theta)^2))
}
# the cauchy pdf

pdf <-function(x,theta,eta,lower.tail=TRUE){
  if(lower.tail) result<-integrate(LMF,lower = -Inf,upper = x,rel.tol=.Machine$double.eps^0.25,theta=theta,eta=eta)
  else result<-integrate(LMF,lower = x,upper = Inf,rel.tol=.Machine$double.eps^0.25,theta=theta,eta=eta)
  return(result$value)
}
pdf(2,2,1,lower.tail = F )
pcauchy(2,location = 1,scale = 2,lower.tail = F)

## ----echo=FALSE----------------------------------------------------------
        dat <- rbind(Genotype=c('AA','BB','OO','AO','BO','AB','Sum'),
             Frequency=c('p2','q2','r2','2pr','2qr','2pq',1),
             Count=c('nAA','nBB','nOO','nAO','nBO','nAB','n'))
    knitr::kable(dat,digits=4,format='html',caption = "Comparation of them",align = "c")
    

## ------------------------------------------------------------------------
#install.packages("nloptr")
library(nloptr)
# mle
logL <- function(x,xpart,n.A,n.B,nOO,nAB) {
  #x[1] p
  #x[2] q
  #xpart[1] p0
  #xpart[2] q0
  r1<-1-sum(xpart)
  nAA<-n.A*xpart[1]^2/(xpart[1]^2+2*xpart[1]*r1)
  nBB<-n.B*xpart[2]^2/(xpart[2]^2+2*xpart[2]*r1)
  r<-1-sum(x)
  return(-2*nAA*log(x[1])-2*nBB*log(x[2])-2*nOO*log(r)-
           (n.A-nAA)*log(2*x[1]*r)-(n.B-nBB)*log(2*x[2]*r)-nAB*log(2*x[1]*x[2]))
}

# constraint function 
limitx <- function(x,xpart,n.A,n.B,nOO,nAB) {
  return(sum(x)-0.999999)
}

opts <- list("algorithm"="NLOPT_LN_COBYLA",
             "xtol_rel"=1.0e-8)
EMmle<-NULL
r<-matrix(0,1,2)
r<-rbind(r,c(0.2,0.2))
# first p0 and q0
j<-2
while (sum((r[j,]-r[j-1,])^2)>1e-10) {
  res <- nloptr( x0=c(0.2,0.2),
                 eval_f=logL,
                 lb = c(0,0), ub = c(1,1), 
                 eval_g_ineq =limitx, 
                 opts = opts, xpart=r[j,],n.A=28,n.B=24,nOO=41,nAB=70 )
  j<-j+1
  r<-rbind(r,res$solution)
  EMmle<-c(EMmle,logL(x=r[j,],xpart=r[j-1,],n.A=28,n.B=24,nOO=41,nAB=70))
}
#Answer:
r[nrow(r),1:2]
#EM:
r
# negetive log likelihood
EMmle

## ------------------------------------------------------------------------
attach(mtcars)
formulas <- list( mpg ~ disp, 
    mpg ~ I(1 / disp), 
    mpg ~ disp + wt,
    mpg ~ I(1 / disp) + wt 
)
#Loops
lmloop <- list()
for (k in 1:length(formulas)) {
  lmloop[[k]] <- lm(formulas[[k]]) 
}
lmloop

#Lapply
lapply(formulas,lm)

#same answers

## ------------------------------------------------------------------------
set.seed(123)
#follow the question:
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ] })

#Loops
lmboots<-list()
for(i in 1:length(bootstraps)){
  lmboots[[i]] <- lm(mpg~disp,data =bootstraps[[i]])
}
lmboots

#Lapply
set.seed(123)
lapply(bootstraps,lm,formula=mpg~disp)
#same seed ,same answer

## ------------------------------------------------------------------------
#From the question: 
rsq <- function(mod) summary.lm(mod)$r.squared
#Q1 Exercises 3
#Loops
rsqloopsQ3 <- NULL
for (i in seq_along(formulas)) {
 rsqloopsQ3[i] <- rsq(lm(formulas[[i]]))
}
rsqloopsQ3
#lapply
lapply(lapply(formulas,lm),rsq)

#Q2 Exercises 4
#Loops
rsqloopsQ4 <- NULL
set.seed(123)
for(i in seq_along(bootstraps)){
 rsqloopsQ4[i] <- rsq(lm(mpg~disp,data =bootstraps[[i]]))
}
rsqloopsQ4

#lapply
set.seed(123)
lapply(lapply(bootstraps,lm,formula=mpg~disp),rsq)

## ------------------------------------------------------------------------
#sapply
trials <- replicate( 100, t.test(rpois(10, 10), rpois(7, 10)), simplify = FALSE )
p_value<-function(mod) mod$p.value
sapply(trials, p_value)
#[[ directly
#??

## ------------------------------------------------------------------------
Mapvapplay<-function (f,n,type, ...) {  
  #n:the length of output 
  f <- match.fun(f)
  fM=Map(f, ...)
  if(type=="numeric")  return(vapply(fM,cbind,numeric(n)))
  else if (type=="character") return(vapply(fM,cbind,character(n)))
  else if (type=="complex") return(vapply(fM,cbind,complex(n)))
  else if (type=="logical") return(vapply(fM,cbind,logical(n)))
}

#examples :Q3
rsq <- function(mod) summary.lm(mod)$r.squared
Mapvapplay(rsq,1,"numeric",lapply(formulas,lm))

#example 2
comparexy <- function(x,y) return(x<y)
x<- list(1:10,2:11,3:12)
y<- list(10:1,11:2,12:3)
Mapvapplay(comparexy,10,"logical",x,y)

## ------------------------------------------------------------------------
library(microbenchmark)
fastchisq.test<-function(x,y){  
  #Pearson's Chi-squared test
  #input two numeric vectors x,y
  #x,y should have same length
  #computes the chi-square test statistic
  X <- rbind(x,y) 
  n<-sum(c(x,y))
  rs <-rowSums(X)
  cs <-colSums(X)
  np <-rs %*% t(cs)/n
  return(sum((X-np)^2/np))
}
x<-c(42,13,24)
y<-c(111,50,50)  
fastchisq.test(x,y)        
chisq.test(rbind(x,y))
microbenchmark(t1=fastchisq.test(x,y),t2=chisq.test(rbind(x,y)))  
#it is obivous that fastchisq.test is faster

## ------------------------------------------------------------------------
fasttable<-function(...,dnn = list.names(...),deparse.level = 1){
  #input: of two integer vectors
  #two integer vectors should have same length
  #output:  
  #sort by quickSort
  list.names <- function(...) {
    l <- as.list(substitute(list(...)))[-1L]
    nm <- names(l)
    fixup <- if (is.null(nm)) 
    seq_along(l)
    else nm == ""
    dep <- vapply(l[fixup], function(x) switch(deparse.level + 
          1, "", if (is.symbol(x)) as.character(x) else "", 
            deparse(x, nlines = 1)[1L]), "")
    if (is.null(nm)) 
      dep
    else {
      nm[fixup] <- dep
      nm
    }
  }
  args <- list(...)
  if (!length(args)) 
    stop("nothing to tabulate")
  if (length(args) == 1L && is.list(args[[1L]])) {
    args <- args[[1L]]
    if (length(dnn) != length(args)) 
    dnn <- if (!is.null(argn <- names(args))) 
    argn
    else paste(dnn[1L], seq_along(args), sep = ".")
  }
  bin <- 0L
  lens <- NULL
  dims <- integer()
  pd <- 1L
  dn <- NULL
  for (a in args) {
    if (is.null(lens))   lens <- length(a)
    else if (length(a) != lens) 
    stop("all arguments must have the same length")
    fact.a <- is.factor(a)
    if (!fact.a) {
      a0 <- a
      a <- factor(a)
    }
    ll <- levels(a)
    a <- as.integer(a)
    nl <- length(ll)
    dims <- c(dims, nl)
    dn <- c(dn, list(ll))
    bin <- bin + pd * (a - 1L)
    pd <- pd * nl
 }
  names(dn) <- dnn
  bin <- bin[!is.na(bin)]
  if (length(bin)) bin <- bin + 1L
  y <- array(tabulate(bin, pd), dims, dimnames = dn)
  class(y) <- "table"
  y
}
# I delete all useless judgement statement

x <- c(rep(1,100),rep(2,100),rep(3,100),rep(4,100),rep(5,100))
y <- c(rep(1:5,100))
fasttable(x,y)
table(x,y)
microbenchmark(t1=fasttable(x,y),t2=table(x,y))    

