library(Rssa)
library(PRIMME)
library(QZ, quiet = TRUE)

outl_error.cssa <- function(ser, N, L, sig, r){
  s <- ssa(ser, L = L, kind = "cssa", svd.method = "primme", neig = r+1)
  rec <- reconstruct(s, groups = list(1:r))
  res <- rec$F1-sig
}

outl_error.theor <- function(l, N, L, a, k) {
  K <- N - L + 1
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

for(N in c(100, 1000, 10000,100000)){
  print(N)
  alpha <- 0.5
  L <- floor(alpha * N)
  K <- N - L + 1
  
  k <- L - 1

  sig <- rep(1 + 1i, N)
  r <- 1
  outl <- rep(0, N)
  outl[k] <- a
  ser <- sig + outl
  
  er.theor <- sapply(1:N, function(i) outl_error.theor(i, N, L, a, k))
  er.real <- outl_error.cssa(ser, N, L, sig,r)
  
  max.diff_n <- c(max.diff_n, max(abs(er.theor - er.real)))
  max.er.theor_n <- c(max.er.theor_n, max(abs(er.theor)))
  max.er.real_n <- c(max.er.real_n, max(abs(er.real)))
}

print(max.diff_n)
print(max.er.theor_n)
print(max.er.real_n)
df <- data.frame(realN2N2 = max.er.real_n, diffN2N2 = max.diff_n)

max.diff_n <- vector()
max.er.theor_n <- vector()
max.er.real_n <- vector()

for(N in c(100, 1000, 10000,100000)){
  print(N)
  L <- 20
  K <- N - L + 1
  
  k <- L - 1
  
  sig <- rep(1 + 1i, N)
  r <- 1
  outl <- rep(0, N)
  outl[k] <- a
  ser <- sig + outl
  
  er.theor <- sapply(1:N, function(i) outl_error.theor(i, N, L, a, k))
  er.real <- outl_error.cssa(ser, N, L, sig,r)
  
  max.diff_n <- c(max.diff_n, max(abs(er.theor - er.real)))
  max.er.theor_n <- c(max.er.theor_n, max(abs(er.theor)))
  max.er.real_n <- c(max.er.real_n, max(abs(er.real)))
}

print(max.diff_n)
print(max.er.theor_n)
print(max.er.real_n)

df$real2020 = max.er.real_n
df$diff2020 = max.diff_n

#############
max.diff_n <- vector()
max.er.theor_n <- vector()
max.er.real_n <- vector()

for(N in c(100, 1000, 10000,100000)){
  print(N)
  alpha <- 0.5
  L <- floor(alpha * N)
  K <- N - L + 1
  
  k <- 20
  
  sig <- rep(1 + 1i, N)
  r <- 1
  outl <- rep(0, N)
  outl[k] <- a
  ser <- sig + outl
  
  er.theor <- sapply(1:N, function(i) outl_error.theor(i, N, L, a, k))
  er.real <- outl_error.cssa(ser, N, L, sig,r)
  
  max.diff_n <- c(max.diff_n, max(abs(er.theor - er.real)))
  max.er.theor_n <- c(max.er.theor_n, max(abs(er.theor)))
  max.er.real_n <- c(max.er.real_n, max(abs(er.real)))
}

print(max.diff_n)
print(max.er.theor_n)
print(max.er.real_n)
df$realN220 = max.er.real_n
df$diffN220 = max.diff_n

max.diff_n <- vector()
max.er.theor_n <- vector()
max.er.real_n <- vector()

for(N in c(100, 1000, 10000,100000)){
  print(N)
  L <- 20
  K <- N - L + 1
  
  k <- N/2 - 1
  
  sig <- rep(1 + 1i, N)
  r <- 1
  outl <- rep(0, N)
  outl[k] <- a
  ser <- sig + outl
  
  er.theor <- sapply(1:N, function(i) outl_error.theor(i, N, L, a, k))
  er.real <- outl_error.cssa(ser, N, L, sig,r)
  
  max.diff_n <- c(max.diff_n, max(abs(er.theor - er.real)))
  max.er.theor_n <- c(max.er.theor_n, max(abs(er.theor)))
  max.er.real_n <- c(max.er.real_n, max(abs(er.real)))
}

print(max.diff_n)
print(max.er.theor_n)
print(max.er.real_n)

df$real20N2 = max.er.real_n
df$diff20N2 = max.diff_n

print(df)
library(xtable)
xtable(df)
print(xtable(df, display=c("s", rep("e",8))), math.style.exponents = TRUE)
print(xtable(t(df), display=c("s", rep("e",4))), math.style.exponents = TRUE)
