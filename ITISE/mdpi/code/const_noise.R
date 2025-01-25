library(Rssa)
library(PRIMME)
library(QZ, quiet = TRUE)

outl_error.cssa <- function(ser, N, L, sig, r){
  s <- ssa(ser, L = L, kind = "cssa", svd.method = "primme")
  rec <- reconstruct(s, groups = list(1:r))
  res <- rec$F1-sig
}

N <- 3999
M <- 1000
L <- 1000
sigma <- 1
#outl[k] <- a

sig <- exp(0.00*(1:N))*(rep(1 + 1i, N))
#sig <- exp(0.01*(1:N))*cos(2*pi*(1:N)/12) + 1i*exp(0.01*(1:N))*cos(2*pi*(1:N)/12 + pi/4)
r <- 1
er.cssa0 <- numeric(N)
set.seed(42)
t0 <- Sys.time()
for(i in 1:M) {
  if(i%/%10*10 == i) {
    t1 <- Sys.time()
    print(i)
    print((t1-t0)/i*(M-i))
  }
  outl <- sigma*(rnorm(N) + 1i*rnorm(N))
  ser <- sig + outl
  K <- N - L + 1
  er.cssa0 <- er.cssa0 + abs(outl_error.cssa(ser, N, L, sig, r))^2
  #er.cssa <- c(er.cssa, sapply(1:N, function(i) outl_error.theor(i, N, L, a, k)))
}

er.cssa0 <- er.cssa0/M
print(e1 <- mean(er.cssa0))
 
pdf("const_noise_err_1.pdf", paper = "special", width = 6, height = 4)
plot(x = 1:N, y = er.cssa0, xlab = "", ylab = "", type = "l")
title(paste0("r = ", r))
dev.off()


sig <- exp(0.00*(1:N))*cos(2*pi*(1:N)/12) + 
  1i*exp(0.00*(1:N))*cos(2*pi*(1:N)/15 + pi/4)
r <- 4
er.cssa <- numeric(N)
set.seed(1)
t0 <- Sys.time()
for(i in 1:M) {
  if(i%/%10*10 == i) {
    t1 <- Sys.time()
    print(i)
    print((t1-t0)/i*(M-i))
  }
  outl <- sigma*(rnorm(N) + 1i*rnorm(N))
  ser <- sig + outl
  K <- N - L + 1
  er.cssa <- er.cssa + abs(outl_error.cssa(ser, N, L, sig, r))^2
  #er.cssa <- c(er.cssa, sapply(1:N, function(i) outl_error.theor(i, N, L, a, k)))
}

er.cssa <- er.cssa/M
print(mean(er.cssa))

pdf("const_noise_err_4.pdf", paper = "special", width = 6, height = 4)
plot(x = 1:N, y = er.cssa, xlab = "", ylab = "", type = "l")
title(paste0("r = ", r))
dev.off()

print(mean(er.cssa)/e1)

er.formula <- numeric(N)
for(i in 0:(N)) {
  lambda <- 2*i/N
  alpha <- L/N
  if(lambda <= 1) er.formula[i] <- calculate_Df(alpha, lambda, sigma1, sigma2)
  else er.formula[i] <- calculate_Df(alpha, 2-lambda, sigma1, sigma2)
}
pdf("const_noise_err_14.pdf", paper = "special", width = 6, height = 4)
plot(x = 1:N, y = (er.cssa)/4, xlab = "", ylab = "", type = "l", col = "grey" )
lines(x = 1:N, y = er.formula, xlab = "", ylab = "", type = "l", col= "red", lwd = 2)
lines(x = 1:N, y = (er.cssa0), xlab = "", ylab = "", type = "l", col= "blue")
legend('topright', c("rank = 4, divided by 4", "rank = 1, first-order", "rank = 1, full"),
       col=c("grey", "red", "blue"), lty=1, cex=0.8, lw=c(2, 2, 2))
#title(paste0("r = ", r))
dev.off()

save.image(file = "1000.Rdata")



