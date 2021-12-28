library(pcaL1)
library(Rssa)
library(matrixcalc)
M <-10
N<-240
iters <- c(1,2,5,10,20)

vec.iter1<-rep(0,length(iters))
MSE.L1svd <- rep(0,M)

for(j in (1:length(iters))){
  
  set.seed(111)
  for(k in 1:M)
  {
    
    
    sig <- exp((1:N)/N)+sin(2*pi*(1:N)/Per+pi/6)
    rnk<-3
    
    sig.outl<-sig
    outlier.seq<-sample(1:(N),N*0.05)
    sig.outl[outlier.seq]<-sig.outl[outlier.seq]+5*sig.outl[outlier.seq]
    ser<-sig.outl+rnorm(N)
    
    plot(ser,type='l')
    X<-hankel(ser, L=120)
    
    s.L1svd<-l1pca(X,center=FALSE,projections="l1",projDim=rnk,iterations=iters[j],tolerance=1e-10)
    Pr<-s.L1svd$projPoints
    Pr.L1svd<-hankL1(Pr)
    
    #MSE
    MSE.L1svd[k] <- mean((sig - Pr.L1svd)^2)
    
    
  }
  vec.iter1[j]<-sqrt(mean(MSE.L1svd))
}

vec.iter1
plot(iters, vec.iter1, type = 'l', xlab ='Number of iterations', ylab = 'RMSE', lwd=2)
title(main='Number of iterations for l1pca, example 1, L=240')
