library(QZ, quiet = TRUE)

weights<-function(R,sigma,a,m,n){
  x<-R/sigma
  W<-matrix(0L, nrow = m, ncol = n)
  W[which(abs(x)<=a,arr.ind=TRUE)]<-(1-(abs(x[which(abs(x)<=a,arr.ind=TRUE)])/a)^2)^2
  return(W)
}

CIRLS_mod<-function(M, k, trend.ver='loess', alpha=4.046, eps=1e-7, maxITER=10, maxiter=5){
  m<-nrow(M)
  n<-ncol(M)
  initial<-svd(M, nu=k, nv=k) #initialization U and V using svd
  U<-initial$u
  U <- as.matrix(U, nrow = nrow(U), ncol = k)
  Lambda<-initial$d
  V<-initial$v
  V <- as.matrix(V, nrow = nrow(), ncol = k)
  U<-U%*%diag(Lambda, nrow = k, ncol = k)
  ITER<-0
  iter<-0
  L<-m #window length for loess, lowess
  
  repeat {
    R<-M-U%*%H(V) #residuals matrix
    r<-as.vector(H(R))
    RR<-hankL1(R)
    RR.trmatrix<-hankel(RR,L=L)
    
    if (trend.ver == 'loess') { 
      loessMod30 <- loess(abs(RR) ~ c(1:length(abs(RR))), span=0.35)
      sigma <- predict(loessMod30)}
    
    else if (trend.ver == 'median') {
      sigma<-runmed(abs(RR),L/2+1)
      sigma[(length(RR)-L/4):length(RR)]<-calc.ends(abs(RR),L/4)}
    
    else if (trend.ver == 'lowess') {
      sigma<-lowess(c(1:length(RR)),abs(RR), f=0.35)$y}
    
    sigma.trmatrix<-hankel(sigma,L=L)
    W<-weights(RR.trmatrix,sigma.trmatrix,alpha,m,n) #weights matrix
    
    
    repeat{
      # updating U
      for (i in (1:m)){
        Wi<-diag(W[i,1:ncol(W)])
        mi<-M[i,1:ncol(M)]
        beta <- qr.solve(H(V)%*%Wi%*%V, H(V)%*%Wi%*%t(H(mi)))
        U[i,1:ncol(U)]<-H(beta)
      }
      # updating V
      for (j in (1:n)){
        Wj<-diag(W[1:nrow(W),j])
        mj<-M[1:nrow(M),j]
        beta <- qr.solve(H(U)%*%Wj%*%U, H(U)%*%Wj%*%mj)
        V[j,1:ncol(V)]<-H(beta)
      }
      iter<-iter+1
      temp <- W^{1/2}*(M-U%*%H(V))
      if ( ((abs(sum(diag(temp%*%H(temp)))))<eps) | (iter > maxiter) ) {break}
    }
    
    ITER<-ITER+1
    temp <- W^{1/2}*(M-U%*%H(V))
    if ( ((abs(sum(diag(temp%*%H(temp)))))<eps) | (ITER > maxITER) ) { break}
  }
  M_est<-U%*%H(V)
  return(M_est)
}