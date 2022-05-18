#L1-hankelization
hankL1<- function(A) {
  N<-nrow(A)+ncol(A)-1
  F<-rep(0,N)
  L.1<-min(nrow(A),ncol(A))
  K.1<-max(nrow(A),ncol(A))
  for(k in 1:N) {
    v<-A[c(row(A) + col(A) - k == 1)]
    F[k]<-median(v)
  }
  return (F)
}

#L2-hankelization
hankL2<- function(A) {
  N<-nrow(A)+ncol(A)-1
  F<-rep(0,N)
  L.1<-min(nrow(A),ncol(A))
  K.1<-max(nrow(A),ncol(A))
  for(k in 1:N) {
    v<-A[c(row(A) + col(A) - k == 1)]
    F[k]<-mean(v)
  }
  return (F)
}

#L2-hankelization with weights 
#(is equal to L2-hankelization if W is hankel matrix)
hankL2_w<- function(A, W) {
  N<-nrow(A)+ncol(A)-1
  F<-rep(0,N)
  L.1<-min(nrow(A),ncol(A))
  K.1<-max(nrow(A),ncol(A))
  for(k in 1:N) {
    v<-A[c(row(A) + col(A) - k == 1)]
    w<-W[c(row(W) + col(W) - k == 1)]
    wv<-w*v
    zero.weights<- w!=rep(0,length(w))
    if (length(which(zero.weights == FALSE))==0){
      F[k]<-sum(wv)/sum(w)
    }
    else {F[k]<-mean(v)}
  }
  return (F)
}