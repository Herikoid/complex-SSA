library(Rssa)
library(PRIMME)
library(QZ, quiet = TRUE)

outl_error.cssa <- function(ser, N){
  s <- ssa(ser, kind = "cssa", svd.method = "primme")
  rec <- reconstruct(s, groups = 1:2)
  res <- rec$F1-sig
}

outl_error.theor <- function(l, N, L, a, k) {
  ## k > K => reverse
  if (k > K) {
    k = N - k + 1
    l = N - l + 1
  }
  ## L <= k <= K
  if (k <= K & k >= L) {
    a = a / L
    if (l <= k & l >= k - L) {
      return(1 / min(c(l, L)) * (l - k + L) * a)
    }
    if (l > k & l < k + L) {
      return(1 / min(c(N - l + 1, L)) * (k + L - l) * a)
    }
    return(0)
  }
  
  a = a / (L * K)
  ## k <= L / 2
  ### k <= K - L
  if (k <= L / 2) {
    if (k <= K - L) {
      if (l <= k ) {
        return((L + K - k)*a)
      }
      if (l < 2*k ) {
        return((1/l)*(L + K - l)* k * a)
      }
      if (l <= L ) {
        return((1/(l))*(L + K - l)*k*a)
      }
      
      if (l < L  + k) {
        return((a/L)*(K*(k + L - l)))
      }
      if (l <= K ) {
        return(0)
      }
      if (l < K + k) {
        return((1 / (N - l + 1)) * (K - l ) * (L - k) * a)
      }
      return(- k * a)
    }
    ### k > K - L
    if (l <= k ) {
      return((L + K - k)*a)
    }
    if (l < 2*k ) {
      return((1/l)*(L + K - l)* k * a)
    }
    if (l <= L ) {
      return((1/(l))*(L + K - l)*k*a)
    }
    
    if (l < K) {
      return((a/L)*(K*(k + L - l)))
    }
    if (l <= L + k ) {
      return(a * (1 / (N - l + 1)) * ((k + K - l)  * (L - k) + (k + L - l)  * (K - k) + (2*k - l) * k))
    }
    if (l < K + k) {
      return((1 / (N - l + 1)) * (K - l ) * (L - k) * a)
    }
    return(- k * a)
  }
  ## k > L / 2
  ### k <= K - L
  if (k <= K - L) {
    if (l <= k ) {
      return((L + K - k)*a)
    }
    if (l < L ) {
      return((1/l)*(L + K - l)* k * a)
    }
    if (l <= 2 * k ) {
      return((a/(L))*(K * (L + k - l)))
    }
    
    if (l < L  + k) {
      return((a/L)*(K*(k + L - l)))
    }
    if (l <= K ) {
      return(0)
    }
    if (l < K + k) {
      return((1 / (N - l + 1)) * (K - l ) * (L - k) * a)
    }
    return(- k * a)
  }
  
  ### K - L < k <= K / 2   
  if (k <= K / 2) {
    if (l <= k ) {
      return((L + K - k)*a)
    }
    if (l < L ) {
      return((1/l)*(L + K - l)* k * a)
    }
    if (l <= 2 * k ) {
      return((a/(L))*(K * (L + k - l)))
    }
    
    if (l < K) {
      return((a/L)*(K*(k + L - l)))
    }
    if (l < L + k ) {
      return(a * (1 / (N - l + 1)) * ((k + K - l)  * (L - k) + (k + L - l)  * (K - k) + (2 * k - l) * k)) ###
    }
    if (l < K + k) {
      return((1 / (N - l + 1)) * (K - l ) * (L - k) * a)
    }
    return(- k * a)
  }
  ### k > K / 2
  if (l <= k ) {
    return((L + K - k)*a)
  }
  if (l < L ) { 
    return((1/l)*(L + K - l)* k * a)
  }
  if (l < K ) {
    return((a/(L))*(K * (L + k - l)))
  }
  
  if (l <= 2 * k) {
    return((a * (1 / (N - l + 1)))*((K - k) * (L - k) + (K - k) * (L - k) + (N - l + 1 - L - K + 2 * k) * (L + K - k)))
  }
  if (l < L + k) {
    return(a * (1 / (N - l + 1)) * (2 * L * K - l * (L + K - k)))###
  }
  if (l < K + k) {
    return((1 / (N - l + 1)) * (K - l ) * (L - k) * a)
  }
  return(- k * a)
}

hankL2 <- function(A) {
  N <- nrow(A) + ncol(A) - 1
  F <- rep(0, N)
  L.1<-min(nrow(A), ncol(A))
  K.1<-max(nrow(A), ncol(A))
  for(k in 1:N) {
    v <- A[c(row(A) + col(A) - k == 1)]
    F[k] <- mean(v)
  }
  return (F)
}

H_1_1 <- function(E, U, V) {
  as.vector((-1) * H(U) %*% E %*% V) * U %*% H(V) + U %*% H(U) %*% E + E %*% V %*% H(V)
}

a <- 10 + 10i

max.diff_n <- vector()
max.er.theor_n <- vector()
max.er.real_n <- vector()

for(N in c(50, 100, 400, 1600)){
  alpha <- 0.5
  L <- floor(alpha * N)
  K <- N - L + 1
  
  k <- L - 1
  
  sig <- rep(1 + 1i, N)
  
  outl <- rep(0, N)
  outl[k] <- a
  ser <- sig + outl
  
  er.theor <- sapply(1:N, function(i) outl_error.theor(i, N, L, a, k))
  er.real <- outl_error.cssa(ser)
  
  max.diff_n <- c(max.diff_n, max(abs(er.theor - er.real)))
  max.er.theor_n <- c(max.er.theor_n, max(abs(er.theor)))
  max.er.real_n <- c(max.er.real_n, max(abs(er.real)))
}

max.diff_n
max.er.theor_n
max.er.real_n

max.diff_n <- vector()
max.er.theor_n <- vector()
max.er.real_n <- vector()

for(N in c(50, 100, 400, 1600)){
  L <- 20
  K <- N - L + 1
  
  k <- L - 1
  
  sig <- rep(1 + 1i, N)
  
  outl <- rep(0, N)
  outl[k] <- a
  ser <- sig + outl
  
  er.theor <- sapply(1:N, function(i) outl_error.theor(i, N, L, a, k))
  er.real <- outl_error.cssa(ser)
  
  max.diff_n <- c(max.diff_n, max(abs(er.theor - er.real)))
  max.er.theor_n <- c(max.er.theor_n, max(abs(er.theor)))
  max.er.real_n <- c(max.er.real_n, max(abs(er.real)))
}

max.diff_n
max.er.theor_n
max.er.real_n

# GRAPHS
er.theor <- vector()

N <- 50
L <- 20
K <- N - L + 1
sig <- rep(1 + 1i, N)
a <- 10 + 10i
for(k in c(1, 5, 10)) {
  
  outl <- rep(0, N)
  outl[k] <- a
  ser <- sig + outl
  
  er.theor <- c(er.theor, sapply(1:N, function(i) outl_error.theor(i, N, L, a, k)))
  er.real <- c(er.real, outl_error.cssa(ser))
}


#pdf("img/const_outl_err_1.pdf", paper = "special", width = 6, height = 4)
plot(x = 1:N, y = Re(er.theor[1:50]), xlab = "", ylab = "", type = "l", ylim = c(min(Re(er.theor)), max(Re(er.theor))))
lines(x = 1:N, y = Re(er.theor[(N+1):(2*N)]), type = "l", col = "blue")
lines(x = 1:N, y = Re(er.theor[(2*N+1):(3*N)]), type = "l", col = "red")
lines(x = c(1, N), y = c(0, 0), type = "l", lty = 2)
legend('topright', c("k = 1", "k = 5", "k = 10"),
       col=c("black", "blue", "red"), lty=1, cex=0.8, lw=c(2, 2))
#dev.off()

er.theor <- vector()

N <- 51
k <- 26
sig <- rep(1 + 1i, N)
a <- 10 + 10i
outl <- rep(0, N)
outl[k] <- a
ser <- sig + outl

for(L in c(10, 20, 25)) {
  K <- N - L + 1
  er.theor <- c(er.theor, sapply(1:N, function(i) outl_error.theor(i, N, L, a, k)))
}


#pdf("img/const_outl_err_2.pdf", paper = "special", width = 6, height = 4)
plot(x = 1:N, y = Re(er.theor[1:N]), xlab = "", ylab = "", type = "l")
lines(x = 1:N, y = Re(er.theor[(N+1):(2*N)]), type = "l", col = "blue")
lines(x = 1:N, y = Re(er.theor[(2*N+1):(3*N)]), type = "l", col = "red")
legend('topright', c("L = 10", "L = 20", "k = 25"),
       col=c("black", "blue", "red"), lty=1, cex=0.8, lw=c(2, 2))
#dev.off()