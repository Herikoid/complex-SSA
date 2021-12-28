#Weights function for IRLS
weights<-function(R,sigma,a,m,n){
  x<-R/sigma
  W<-matrix(0L, nrow = m, ncol = n)
  W[which(abs(x)<=a,arr.ind=TRUE)]<-(1-(abs(x[which(abs(x)<=a,arr.ind=TRUE)])/a)^2)^2
  return(W)
}

calc.ends<-function(vec,L){
  N<-length(vec)
  res<-rep(NaN, (L+1))
  for (i in (1:(L+1))) {
    slide<-vec[(N-2*L+i):N]
    res[i]<-median(slide)
  }
  return(res)
}