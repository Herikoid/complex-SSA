source("hankelization.R")
source("IRLS_complex.R")
source("l1_complex_v2.R")
source("CIRLS_modified.R")

library(Rssa)
library(matrixcalc)
library(cmvnorm)
library(pcaL1)

N<-240
Per<-120

set.seed(7)
sig <- exp(4*(1:N)/N)*exp(2i*pi*(1:N)/30)
rnk<-1
sig.outl<-sig
outlier.seq<-sample(1:(N),N*0.01)
sig.outl[outlier.seq]<-sig.outl[outlier.seq] + 5*sig.outl[outlier.seq] 

ser<-sig.outl + 0.5*exp(4*(1:N)/N)*rcnorm(N)


plot(Re(ser), col='black', type='l')
legend('topleft', c( "series"), col=c("black"), lty=1, cex=0.8, lw=c(2))

plot(Im(ser), col='black', type='l')
legend('topleft', c( "series"), col=c("black"), lty=1, cex=0.8, lw=c(2))


X<-hankel(ser, L=120)

Pr<-IRLS_complex(X, rnk)
Pr0<-hankL2(Pr)

Pr<-l1_complex(X, rnk)
Pr.L1<-hankL1(Pr)

Pr<-CIRLS_mod(X,rnk,'loess') #CIRLS modification (trend extraction with loess)
Pr1<-hankL2(Pr)

Pr<-CIRLS_mod(X,rnk,'median') #CIRLS modification (trend extraction with median)
Pr2<-hankL2(Pr)

Pr<-CIRLS_mod(X,rnk,'lowess') #CIRLS modification (trend extraction with lowess)
Pr3<-hankL2(Pr)


s <- ssa(ser, kind = "cssa")
r <- reconstruct(s, groups = list(Trend = 1:rnk))

plot(Re(ser),col='black',type='l')
lines( Re(sig), type='l', col='yellow', lw=2)
lines( Re(r$Trend),type='l',col='blue',lw=2)
lines( Re(Pr0),type='l',col='red',lw=2)
lines( Re(Pr.L1), type='l', col='green', lw=2)
legend('topleft', c("series", "target", "cssa", "l2-weighted", "l1"),
       col=c("black", "yellow", "blue", "red", "green"), lty=1, cex=0.8, lw=c(2, 2, 2, 2, 2))



plot(Im(ser),col='black',type='l')
lines( Im(sig), type='l', col='yellow', lw=2)
lines( Im(r$Trend),type='l',col='blue',lw=2)
lines( Im(Pr0),type='l',col='red',lw=2)
lines( Im(Pr.L1), type='l', col='green', lw=2)
legend('topleft', c("series", "target", "cssa", "l2-weight", "l1"),
       col=c("black", "yellow", "blue", "red", "green"), lty=1, cex=0.8, lw=c(2, 2, 2, 2, 2))

print("CSSA RMSE:")
print(sqrt(mean((Re(sig) - Re(r$Trend))^2) + mean((Im(sig) - Im(r$Trend))^2)))
print("L1 RMSE:")
print(sqrt(mean((Re(sig) - Re(Pr.L1))^2) + mean((Im(sig) - Im(Pr.L1))^2)))
print("w-L2 RMSE:")
print(sqrt(mean((Re(sig) - Re(Pr0))^2) + mean((Im(sig) - Im(Pr0))^2)))
print("w-L2 loess RMSE:")
print(sqrt(mean((Re(sig) - Re(Pr1))^2) + mean((Im(sig) - Im(Pr1))^2)))
print("w-L2 median RMSE:")
print(sqrt(mean((Re(sig) - Re(Pr2))^2) + mean((Im(sig) - Im(Pr2))^2)))
print("w-L2 lowess RMSE:")
print(sqrt(mean((Re(sig) - Re(Pr3))^2) + mean((Im(sig) - Im(Pr3))^2)))


###############################################################
# Вычисление RMSE и p-value, слишком долго, лучше не запускать #
###############################################################
a1 <- vector()
a2 <- vector()
a3 <- vector()
a4 <- vector()
a5 <- vector()
a6 <- vector()

b1 <- vector()
b2 <- vector()
b3 <- vector()
b4 <- vector()
b5 <- vector()
b6 <- vector()

c1 <- vector()
c2 <- vector()
c3 <- vector()
c4 <- vector()
c5 <- vector()
c6 <- vector()

e1 <- vector()
e2 <- vector()
e3 <- vector()
e4 <- vector()
e5 <- vector()
e6 <- vector()

for (i in 1:30) {
  sig <- exp(4*(1:N)/N)*exp(2i*pi*(1:N)/30)
  rnk<-1
  sig.outl<-sig
  outlier.seq<-sample(1:(N),N*0.05)
  sig.outl[outlier.seq]<-sig.outl[outlier.seq] + 5*sig.outl[outlier.seq] 
  
  ser<-sig.outl + 0.5*exp(4*(1:N)/N)*rcnorm(N)
  
  
  X<-hankel(ser, L=120)
  
  Pr<-IRLS_complex(X, rnk)
  Pr0<-hankL2(Pr)
  
  Pr<-l1_complex(X, rnk) #l1pca
  Pr.L1<-hankL1(Pr)
  
  Pr<-CIRLS_mod(X,rnk,'loess') #CIRLS modification (trend extraction with loess)
  Pr1<-hankL2(Pr)
  
  Pr<-CIRLS_mod(X,rnk,'median') #CIRLS modification (trend extraction with median)
  Pr2<-hankL2(Pr)
  
  Pr<-CIRLS_mod(X,rnk,'lowess') #CIRLS modification (trend extraction with lowess)
  Pr3<-hankL2(Pr)
  
  
  s <- ssa(ser, kind = "cssa")
  r <- reconstruct(s, groups = list(Trend = 1))
  
  a1 <- c(a1, sqrt(mean((Re(sig) - Re(r$Trend))^2) + mean((Im(sig) - Im(r$Trend))^2)))
  a2 <- c(a2, sqrt(mean((Re(sig) - Re(Pr.L1))^2) + mean((Im(sig) - Im(Pr.L1))^2)))
  a3 <- c(a3, sqrt(mean((Re(sig) - Re(Pr0))^2) + mean((Im(sig) - Im(Pr0))^2)))
  a4 <- c(a4, sqrt(mean((Re(sig) - Re(Pr1))^2) + mean((Im(sig) - Im(Pr1))^2)))
  a5 <- c(a5, sqrt(mean((Re(sig) - Re(Pr2))^2) + mean((Im(sig) - Im(Pr2))^2)))
  a6 <- c(a6, sqrt(mean((Re(sig) - Re(Pr3))^2) + mean((Im(sig) - Im(Pr3))^2)))
  
  sig <- exp(4*(1:N)/N)*exp(2i*pi*(1:N)/30)
  rnk<-1
  sig.outl<-sig
  #outlier.seq<-sample(1:(N),N*0.05)
  #sig.outl[outlier.seq]<-sig.outl[outlier.seq] + 5*sig.outl[outlier.seq] 
  
  ser<-sig.outl + 0.5*exp(4*(1:N)/N)*rcnorm(N)
  
  
  X<-hankel(ser, L=120)
  
  Pr<-IRLS_complex(X, rnk)
  Pr0<-hankL2(Pr)
  
  Pr<-l1_complex(X, rnk) #l1pca
  Pr.L1<-hankL1(Pr)
  
  Pr<-CIRLS_mod(X,rnk,'loess') #CIRLS modification (trend extraction with loess)
  Pr1<-hankL2(Pr)
  
  Pr<-CIRLS_mod(X,rnk,'median') #CIRLS modification (trend extraction with median)
  Pr2<-hankL2(Pr)
  
  Pr<-CIRLS_mod(X,rnk,'lowess') #CIRLS modification (trend extraction with lowess)
  Pr3<-hankL2(Pr)
  
  
  s <- ssa(ser, kind = "cssa")
  r <- reconstruct(s, groups = list(Trend = 1))
  
  b1 <- c(b1, sqrt(mean((Re(sig) - Re(r$Trend))^2) + mean((Im(sig) - Im(r$Trend))^2)))
  b2 <- c(b2, sqrt(mean((Re(sig) - Re(Pr.L1))^2) + mean((Im(sig) - Im(Pr.L1))^2)))
  b3 <- c(b3, sqrt(mean((Re(sig) - Re(Pr0))^2) + mean((Im(sig) - Im(Pr0))^2)))
  b4 <- c(b4, sqrt(mean((Re(sig) - Re(Pr1))^2) + mean((Im(sig) - Im(Pr1))^2)))
  b5 <- c(b5, sqrt(mean((Re(sig) - Re(Pr2))^2) + mean((Im(sig) - Im(Pr2))^2)))
  b6 <- c(b6, sqrt(mean((Re(sig) - Re(Pr3))^2) + mean((Im(sig) - Im(Pr3))^2)))
  
  sig <- exp(2i*pi*(1:N)/30)
  rnk<-1
  sig.outl<-sig
  outlier.seq<-sample(1:(N),N*0.05)
  sig.outl[outlier.seq]<-sig.outl[outlier.seq] + 5*sig.outl[outlier.seq] 
  
  ser<-sig.outl + rcnorm(N)
  
  
  X<-hankel(ser, L=120)
  
  Pr<-IRLS_complex(X, rnk)
  Pr0<-hankL2(Pr)
  
  Pr<-l1_complex(X, rnk) #l1pca
  Pr.L1<-hankL1(Pr)
  
  Pr<-CIRLS_mod(X,rnk,'loess') #CIRLS modification (trend extraction with loess)
  Pr1<-hankL2(Pr)
  
  Pr<-CIRLS_mod(X,rnk,'median') #CIRLS modification (trend extraction with median)
  Pr2<-hankL2(Pr)
  
  Pr<-CIRLS_mod(X,rnk,'lowess') #CIRLS modification (trend extraction with lowess)
  Pr3<-hankL2(Pr)
  
  
  s <- ssa(ser, kind = "cssa")
  r <- reconstruct(s, groups = list(Trend = 1))
  
  c1 <- c(c1, sqrt(mean((Re(sig) - Re(r$Trend))^2) + mean((Im(sig) - Im(r$Trend))^2)))
  c2 <- c(c2, sqrt(mean((Re(sig) - Re(Pr.L1))^2) + mean((Im(sig) - Im(Pr.L1))^2)))
  c3 <- c(c3, sqrt(mean((Re(sig) - Re(Pr0))^2) + mean((Im(sig) - Im(Pr0))^2)))
  c4 <- c(c4, sqrt(mean((Re(sig) - Re(Pr1))^2) + mean((Im(sig) - Im(Pr1))^2)))
  c5 <- c(c5, sqrt(mean((Re(sig) - Re(Pr2))^2) + mean((Im(sig) - Im(Pr2))^2)))
  c6 <- c(c6, sqrt(mean((Re(sig) - Re(Pr3))^2) + mean((Im(sig) - Im(Pr3))^2)))
  
  sig <- exp(2i*pi*(1:N)/30)
  rnk<-1
  sig.outl<-sig
  #outlier.seq<-sample(1:(N),N*0.05)
  #sig.outl[outlier.seq]<-sig.outl[outlier.seq] + 5*sig.outl[outlier.seq] 
  
  ser<-sig.outl + rcnorm(N)
  
  
  X<-hankel(ser, L=120)
  
  Pr<-IRLS_complex(X, rnk)
  Pr0<-hankL2(Pr)
  
  Pr<-l1_complex(X, rnk) #l1pca
  Pr.L1<-hankL1(Pr)
  
  Pr<-CIRLS_mod(X,rnk,'loess') #CIRLS modification (trend extraction with loess)
  Pr1<-hankL2(Pr)
  
  Pr<-CIRLS_mod(X,rnk,'median') #CIRLS modification (trend extraction with median)
  Pr2<-hankL2(Pr)
  
  Pr<-CIRLS_mod(X,rnk,'lowess') #CIRLS modification (trend extraction with lowess)
  Pr3<-hankL2(Pr)
  
  
  s <- ssa(ser, kind = "cssa")
  r <- reconstruct(s, groups = list(Trend = 1))
  
  e1 <- c(e1, sqrt(mean((Re(sig) - Re(r$Trend))^2) + mean((Im(sig) - Im(r$Trend))^2)))
  e2 <- c(e2, sqrt(mean((Re(sig) - Re(Pr.L1))^2) + mean((Im(sig) - Im(Pr.L1))^2)))
  e3 <- c(e3, sqrt(mean((Re(sig) - Re(Pr0))^2) + mean((Im(sig) - Im(Pr0))^2)))
  e4 <- c(e4, sqrt(mean((Re(sig) - Re(Pr1))^2) + mean((Im(sig) - Im(Pr1))^2)))
  e5 <- c(e5, sqrt(mean((Re(sig) - Re(Pr2))^2) + mean((Im(sig) - Im(Pr2))^2)))
  e6 <- c(e6, sqrt(mean((Re(sig) - Re(Pr3))^2) + mean((Im(sig) - Im(Pr3))^2)))
}

print(a1/10)
print(a2/10)
print(a3/10)
print(a4/10)
print(a5/10)
print(a6/10)


library(tictoc)


N<-720
Per<-120

set.seed(7)
sig <- exp(4*(1:N)/N)*exp(2i*pi*(1:N)/30)
rnk<-1
sig.outl<-sig
outlier.seq<-sample(1:(N),N*0.01)
sig.outl[outlier.seq]<-sig.outl[outlier.seq] + 5*sig.outl[outlier.seq] 

ser<-sig.outl + 2*exp(4*(1:N)/N)*rcnorm(N)

X<-hankel(ser, L=120)

tic()
s <- ssa(ser, kind = "cssa")
r <- reconstruct(s, groups = list(Trend = 1))
toc()

tic()
Pr<-l1_complex(X, rnk)
toc()

tic()
Pr<-IRLS_complex(X, rnk)
toc()

p_val_rmse <- function(x, y, n = 30) {
  sd_x <- sd(x)
  sd_y <- sd(y)
  t <- sqrt(n) * (mean(x) - mean(y)) / sqrt(sd_x^2 + sd_y^2 - 2 * sd_x * sd_y * cor(x, y))
  if (t > 0) {
    return(2 - 2 *pnorm(t))
  }
  return(2*pnorm(t))
}


f1 <- vector()
f2 <- vector()
f3 <- vector()
f4 <- vector()
f5 <- vector()
f6 <- vector()

ref1 <- vector()
ref2 <- vector()
ref3 <- vector()
ref4 <- vector()
ref5 <- vector()
ref6 <- vector()

imf1 <- vector()
imf2 <- vector()
imf3 <- vector()
imf4 <- vector()
imf5 <- vector()
imf6 <- vector()

g1 <- vector()
g2 <- vector()
g3 <- vector()
g4 <- vector()
g5 <- vector()
g6 <- vector()

gref1 <- vector()
gref2 <- vector()
gref3 <- vector()
gref4 <- vector()
gref5 <- vector()
gref6 <- vector()

gimf1 <- vector()
gimf2 <- vector()
gimf3 <- vector()
gimf4 <- vector()
gimf5 <- vector()
gimf6 <- vector()


N<-240
Per<-120
set.seed(1)

for (i in 1:30) {
  sig <- exp(4*(1:N)/N)*exp(2i*pi*(1:N)/30)
  rnk<-1
  sig.outl<-sig
  sig.outl1<-sig
  sig.outl2<-sig
  outlier.seq<-sample(1:(N),N*0.05)
  sig.outl[outlier.seq]<-sig.outl[outlier.seq] + 4 * sig.outl[outlier.seq]
  
  abs_sum <- sum(abs(sig.outl[outlier.seq]))
  sq_sum <- sum(abs(sig.outl[outlier.seq])^2)
  
  
  sig.outl1[outlier.seq] <- sig.outl1[outlier.seq] + sqrt(25 * Re(sig.outl1[outlier.seq])^2 + 24 * Im(sig.outl1[outlier.seq])^2) - Re(sig.outl1[outlier.seq])
  #sig.outl2[outlier.seq] <- sig.outl2[outlier.seq] + 4 * abs(sig.outl[outlier.seq])
  
  noise <- 0.5*exp(4*(1:N)/N)*rcnorm(N)
  ser<-sig.outl + noise
  ser1 <- sig.outl1 + noise
    

  X<-hankel(ser, L=120)

  Pr<-IRLS_complex(X, rnk)
  Pr0<-hankL2(Pr)

  Pr<-l1_complex(X, rnk) #l1pca
  Pr.L1<-hankL1(Pr)

  Pr<-CIRLS_mod(X,rnk,'loess') #CIRLS modification (trend extraction with loess)
  Pr1<-hankL2(Pr)

  #Pr<-CIRLS_mod(X,rnk,'median') #CIRLS modification (trend extraction with median)
  #Pr2<-hankL2(Pr)

  Pr<-CIRLS_mod(X,rnk,'lowess') #CIRLS modification (trend extraction with lowess)
  Pr3<-hankL2(Pr)

  X<-hankel(ser1, L=120)
  
  Pr<-IRLS_complex(X, rnk)
  Pr01<-hankL2(Pr)
  
  Pr<-l1_complex(X, rnk) #l1pca
  Pr.L11<-hankL1(Pr)
  
  Pr<-CIRLS_mod(X,rnk,'loess') #CIRLS modification (trend extraction with loess)
  Pr11<-hankL2(Pr)
  
  #Pr<-CIRLS_mod(X,rnk,'median') #CIRLS modification (trend extraction with median)
  #Pr21<-hankL2(Pr)
  
  Pr<-CIRLS_mod(X,rnk,'lowess') #CIRLS modification (trend extraction with lowess)
  Pr31<-hankL2(Pr)
  
  s <- ssa(ser, kind = "cssa")
  r <- reconstruct(s, groups = list(Trend = 1))


  f1 <- c(f1, sqrt(mean((Re(sig) - Re(r$Trend))^2) + mean((Im(sig) - Im(r$Trend))^2)))
  f2 <- c(f2, sqrt(mean((Re(sig) - Re(Pr.L1))^2) + mean((Im(sig) - Im(Pr.L1))^2)))
  f3 <- c(f3, sqrt(mean((Re(sig) - Re(Pr0))^2) + mean((Im(sig) - Im(Pr0))^2)))
  f4 <- c(f4, sqrt(mean((Re(sig) - Re(Pr1))^2) + mean((Im(sig) - Im(Pr1))^2)))
  #f5 <- c(f5, sqrt(mean((Re(sig) - Re(Pr2))^2) + mean((Im(sig) - Im(Pr2))^2)))
  f6 <- c(f6, sqrt(mean((Re(sig) - Re(Pr3))^2) + mean((Im(sig) - Im(Pr3))^2)))

  ref1 <- c(ref1, sqrt(mean((Re(sig) - Re(r$Trend))^2)))
  ref2 <- c(ref2, sqrt(mean((Re(sig) - Re(Pr.L1))^2)))
  ref3 <- c(ref3, sqrt(mean((Re(sig) - Re(Pr0))^2)))
  ref4 <- c(ref4, sqrt(mean((Re(sig) - Re(Pr1))^2)))
  #ref5 <- c(ref5, sqrt(mean((Re(sig) - Re(Pr2))^2)))
  ref6 <- c(ref6, sqrt(mean((Re(sig) - Re(Pr3))^2)))

  imf1 <- c(imf1, sqrt(mean((Im(sig) - Im(r$Trend))^2)))
  imf2 <- c(imf2, sqrt(mean((Im(sig) - Im(Pr.L1))^2)))
  imf3 <- c(imf3, sqrt(mean((Im(sig) - Im(Pr0))^2)))
  imf4 <- c(imf4, sqrt(mean((Im(sig) - Im(Pr1))^2)))
  #imf5 <- c(imf5, sqrt(mean((Im(sig) - Im(Pr2))^2)))
  imf6 <- c(imf6, sqrt(mean((Im(sig) - Im(Pr3))^2)))
  
  g1 <- c(g1, sqrt(mean((Re(sig) - Re(r$Trend))^2) + mean((Im(sig) - Im(r$Trend))^2)))
  g2 <- c(g2, sqrt(mean((Re(sig) - Re(Pr.L11))^2) + mean((Im(sig) - Im(Pr.L11))^2)))
  g3 <- c(g3, sqrt(mean((Re(sig) - Re(Pr01))^2) + mean((Im(sig) - Im(Pr01))^2)))
  g4 <- c(g4, sqrt(mean((Re(sig) - Re(Pr11))^2) + mean((Im(sig) - Im(Pr11))^2)))
  #g5 <- c(g5, sqrt(mean((Re(sig) - Re(Pr21))^2) + mean((Im(sig) - Im(Pr21))^2)))
  g6 <- c(g6, sqrt(mean((Re(sig) - Re(Pr31))^2) + mean((Im(sig) - Im(Pr31))^2)))
  
  gref1 <- c(gref1, sqrt(mean((Re(sig) - Re(r$Trend))^2)))
  gref2 <- c(gref2, sqrt(mean((Re(sig) - Re(Pr.L11))^2)))
  gref3 <- c(gref3, sqrt(mean((Re(sig) - Re(Pr01))^2)))
  gref4 <- c(gref4, sqrt(mean((Re(sig) - Re(Pr11))^2)))
  #gref5 <- c(gref5, sqrt(mean((Re(sig) - Re(Pr21))^2)))
  gref6 <- c(gref6, sqrt(mean((Re(sig) - Re(Pr31))^2)))
  
  gimf1 <- c(gimf1, sqrt(mean((Im(sig) - Im(r$Trend))^2)))
  gimf2 <- c(gimf2, sqrt(mean((Im(sig) - Im(Pr.L11))^2)))
  gimf3 <- c(gimf3, sqrt(mean((Im(sig) - Im(Pr01))^2)))
  gimf4 <- c(gimf4, sqrt(mean((Im(sig) - Im(Pr11))^2)))
  #gimf5 <- c(gimf5, sqrt(mean((Im(sig) - Im(Pr21))^2)))
  gimf6 <- c(gimf6, sqrt(mean((Im(sig) - Im(Pr31))^2)))
  
}
