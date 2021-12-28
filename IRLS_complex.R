#Original IRLS method
# M --- trajectory matrix
# k --- rank of signal
# maxITER --- number of iterations N_IRLS
# maxiter --- number of iterations N_alpha

library(QZ, quiet = TRUE)

weights<-function(R,sigma,a,m,n){
  x<-R/sigma
  W<-matrix(0L, nrow = m, ncol = n)
  W[which(abs(x)<=a,arr.ind=TRUE)]<-(1-(abs(x[which(abs(x)<=a,arr.ind=TRUE)])/a)^2)^2
  return(W)
}

IRLS_complex<-function(M, k, alpha=4.685, eps=1e-5, maxITER=10, maxiter=5){
  m<-nrow(M)
  n<-ncol(M)
  initial<-svd(M, nu=k, nv=k) #initialization U and V using svd
  U<-initial$u
  U <- as.matrix(U, nrow = nrow(U), ncol = k)
  Lambda<-initial$d
  V<-initial$v
  U<-U%*%diag(Lambda, nrow = k, ncol = k)
  ITER<-0
  iter<-0
  V <- as.matrix(V, nrow = nrow(), ncol = k)
  repeat {
    R<-M-U%*%H(V) #residuals matrix
    r<-as.vector(H(R))
    sigma<-1.4826*median(abs(r-median(abs(r))))
    W<-weights(R,sigma,alpha,m,n)
    #W <- matrix(1, m, n)
    
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