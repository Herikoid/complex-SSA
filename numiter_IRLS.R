library(pcaL1)
library(Rssa)
library(matrixcalc)
M <-10
N<-240
iters <- c(1,2,5,10,20,40)

vec.iter1<-rep(0,length(iters))
MSE.IRLS <- rep(0,M)

for(j in (1:length(iters))){
  print(j)
  set.seed(111)
  for(k in 1:M)
  {
    
    sig <- exp((1:N)/N)+sin(2*pi*(1:N)/Per+pi/6)
    rnk<-3
    
    #sig <- exp(4*(1:N)/N)*sin(2*pi*(1:N)/30) # 2 example
    #rnk<-2
    
    
    sig.outl<-sig
    outlier.seq<-sample(1:(N),N*0.05)
    sig.outl[outlier.seq]<-sig.outl[outlier.seq]+5*sig.outl[outlier.seq]
    ser<-sig.outl+rnorm(N)
    #ser<-sig.outl+0.5*exp(4*(1:N)/N)*rnorm(N)
    
    plot(ser,type='l')
    X<-hankel(ser, L=120)
    
    Pr<-IRLS_orig(X,rnk, maxITER=iters[j] , maxiter = 2, eps=1e-10)
    Pr.IRLS<-hankL2(Pr)
    
    #Pr<-IRLS_mod(X,rnk, "loess", maxITER=iters[j] , maxiter = 2, eps=1e-10)
    #Pr.IRLS<-hankL2(Pr)
    
    #MSE
    MSE.IRLS[k] <- mean((sig - Pr.IRLS)^2)
    
  }
  vec.iter1[j]<-sqrt(mean(MSE.IRLS))
}

vec.iter1
plot(iters, vec.iter1, type = 'l', xlab ='Number of iterations', ylab = 'RMSE', lwd=2)
title(main='Number of iterations for IRLS, example 1, N=240')