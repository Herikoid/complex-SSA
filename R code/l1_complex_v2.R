library(BBmisc)
library(spatstat)
library(QZ, quiet = TRUE)


weighted_median <- function(x, w) {
  S_2 <- sum(w) / 2
  temp <- data.frame(x_abs = abs(x), w = w, x = x)
  temp <- temp[order(temp$x_abs),]
  k <- 1
  s <- temp$w[k]
  while(s < S_2) {
    k <- k + 1
    s <- s + temp$w[k]
  }
  
  return(temp$x[k])
}

l1_complex<-function(M, k, eps=1e-5, maxiter=10){
  m<-nrow(M)
  n<-ncol(M)
  initial<-svd(M, nu=k, nv=k) #initialization U and V using svd
  U<-initial$u
  U <- as.matrix(U, nrow = nrow(U), ncol = k)
  V<-initial$v
  V <- as.matrix(V, nrow = nrow(), ncol = k)
  V <- normalize(V, method = "scale", margin = 2)
  ITER<-0
  iter<-0
    
  repeat{
    V_old <- V
    
    # updating U
    for (i in (1:m)){
      mi<-M[i,1:ncol(M)]
      for (c in 1:k) {
        vc <- V[1:nrow(V), c]
        beta <- weighted.median(Re(t(H(mi))/vc), abs(vc)) + complex(imaginary = weighted.median(Im(t(H(mi))/vc), abs(vc)))
        U[i, c]<-H(beta) 
      }
    }
    
    # updating V
    for (j in (1:n)){
      mj<-M[1:nrow(M),j]
      for (c in 1:k) {
        uc <- U[1:nrow(U), c]
        beta <- weighted.median(Re(mj/uc), abs(uc)) +  complex(imaginary = weighted.median(Im(mj/uc), abs(uc)))
        V[j, c]<-H(beta) 
      }
    }
    V <- normalize(V, method = "scale", margin = 2)
  
    iter<-iter+1
    if  ((max(abs(Re(V) - Re(V_old)) + abs(Im(V) - Im(V_old)))<eps) | (iter > maxiter) ) {break}
  }
  
  M_est<-U%*%H(V)
  return(M_est)
}