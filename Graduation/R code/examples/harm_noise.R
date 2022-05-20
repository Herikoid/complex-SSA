library(Rssa)
library(QZ, quiet = TRUE)

N <- 9
L <- 5

sig.n <- function(N) {
  cos(2 * pi * (1:N) / 10) +  1i * cos(2 * pi * (1:N) / 10 + pi/4)
}

sig <-  sig.n(N)

noise <- function(N) {
  0.1 * (rnorm(N, sd = 1) + 1i * rnorm(N, sd = 1))
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

H.1 <- function(E, U, V) {
  Pcol <- U[[1]] %*% H(U[[1]]) +  U[[2]] %*% H(U[[2]]) 
  Prow <- V[[1]] %*% H(V[[1]]) +  V[[2]] %*% H(V[[2]]) 
  Pcol.orth <- diag(rep(1,length(U[[1]]))) - Pcol
  
  Pcol %*% E + Pcol.orth %*% E %*% Prow
}

H.1.1 <- function(E, U, V) {
  Pcol <- U %*% H(U) 
  Prow <- V %*% H(V) 
  Pcol.orth <- diag(rep(1,length(U))) - Pcol
  
  Pcol %*% E + Pcol.orth %*% E %*% Prow
}

s <- ssa(sig, L = L, kind = "cssa", svd.method = "svd")

rec <- reconstruct(s, groups = list(1:2))[[1]]
U <- list(s$U[,1] + 1i * 0, s$U[,2] + 1i * 0)
V <- list(s$V[,1] + 1i * 0, s$V[,2] + 1i * 0)

set.seed(6)
err.1 <- hankL2(H.1(Rssa::hankel(noise(N), L), U, V))

set.seed(6)
s <- ssa(sig + noise(N), L = L, kind = "cssa", svd.method = "svd")
rec <- reconstruct(s, groups = list(1:2))[[1]]

err <- rec - sig

#pdf(".../img/first_vs_full_re.pdf", paper = "special", width = 6, height = 4) 
plot(Re(err.1), type='l', col = "red", xlab = "l", ylab = "error")
lines(Re(err), type = 'l', col = "blue")
legend('topright', c("first error", "full error"),
       col=c("red", "blue"), lty=1, cex=0.8, lw=c(2, 2))
#dev.off()

#table
diff.1 <- vector()
diff.k <- vector()
set.seed(1)
for(N in c(50, 100, 400, 1600)) {
  L = floor(N / 2)
  k = floor(L / 2)
  sig = sig.n(N)
  noise = noise.n(N)
  
  s <- ssa(sig, L = L, kind = "cssa", svd.method = "svd")
  
  rec <- reconstruct(s, groups = list(1:2))[[1]]
  U <- list(s$U[,1] + 1i * 0, s$U[,2] + 1i * 0)
  V <- list(s$V[,1] + 1i * 0, s$V[,2] + 1i * 0)
  err.1 <- hankL2(H.1(Rssa::hankel(noise, L), U, V))
  
  s <- ssa(sig + noise, L = L, kind = "cssa", svd.method = "svd")
  rec <- reconstruct(s, groups = list(1:2))[[1]]
  
  err <- rec - sig
  
  diff <- abs(err.1 - err)
  
  diff.1 <- c(diff.1, diff[1])
  diff.k <- c(diff.k, diff[k])
}

diff.1
diff.k