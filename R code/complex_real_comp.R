library(Rssa)
library(matrixcalc)
library(cmvnorm)
library(pcaL1)

N<-50
L <- 20
repeats <- 100
sig <- (cos(2*pi*(1:N)/10) + 1i*cos(2*pi*(1:N)/10 + pi/2))
crnk <- 1
rnk<-2


set.seed(7)

a1 <- vector()
a2 <- vector()
a3 <- vector()
a4 <- vector()
a5 <- vector()

a1.Re <- vector()
a2.Re <- vector()
a3.Re <- vector()
a4.Re <- vector()
a5.Re <- vector()

a1.Im <- vector()
a2.Im <- vector()
a3.Im <- vector()
a4.Im <- vector()
a5.Im <- vector()

for(i in 1:repeats) {
  sig.outl<-sig
  outlier.seq<-sample(1:(N),N*0.05)
  sig.outl[outlier.seq]<-sig.outl[outlier.seq] + 5*sig.outl[outlier.seq] 
  noise <- 0.1 * rcnorm(N)
  
  ser<-sig.outl + noise
  
  X<-hankel(ser, L=L)
  
  Pr<-IRLS_complex(X, crnk)
  Pr0<-hankL2(Pr)
  
  Pr<-CIRLS_mod(X,crnk,'loess') #CIRLS modification (trend extraction with loess)
  Pr1<-hankL2(Pr)
  
  Pr<-CIRLS_mod(X,crnk,'median') #CIRLS modification (trend extraction with median)
  Pr2<-hankL2(Pr)
  
  Pr<-CIRLS_mod(X,crnk,'lowess') #CIRLS modification (trend extraction with lowess)
  Pr3<-hankL2(Pr)
  
  
  s <- ssa(ser, L = L, kind = "cssa")
  r <- reconstruct(s, groups = list(1:crnk))[[1]]
  
  X<-hankel(Re(ser), L=L)
  
  Pr<-IRLS_mod(X,rnk,'loess') #IRLS modification (trend extraction with loess)
  Pr1.Re<-hankL2(Pr)
  
  Pr<-IRLS_mod(X,rnk,'median') #IRLS modification (trend extraction with median)
  Pr2.Re<-hankL2(Pr)
  
  Pr<-IRLS_mod(X,rnk,'lowess') #IRLS modification (trend extraction with lowess)
  Pr3.Re<-hankL2(Pr)
  
  Pr<-IRLS_orig(X,rnk) #IRLS original
  Pr0.Re<-hankL2(Pr)
  
  
  s <- ssa(Re(ser), L = L, kind = "1d-ssa")
  r.Re <- reconstruct(s, groups = list(1:rnk))[[1]]
  
  X<-hankel(Im(ser), L=L)
  
  Pr<-IRLS_mod(X,rnk,'loess') #IRLS modification (trend extraction with loess)
  Pr1.Im<-hankL2(Pr)
  
  Pr<-IRLS_mod(X,rnk,'median') #IRLS modification (trend extraction with median)
  Pr2.Im<-hankL2(Pr)
  
  Pr<-IRLS_mod(X,rnk,'lowess') #IRLS modification (trend extraction with lowess)
  Pr3.Im<-hankL2(Pr)
  
  Pr<-IRLS_orig(X,rnk) #IRLS original
  Pr0.Im<-hankL2(Pr)
  
  s <- ssa(Im(ser), L = L, kind = "1d-ssa")
  r.Im <- reconstruct(s, groups = list(1:rnk))[[1]]
  
  a1 <- c(a1, mean(abs(r - sig)^2))
  a2 <- c(a2, mean(abs(Pr0 - sig)^2))
  a3 <- c(a3, mean(abs(Pr1 - sig)^2))
  a4 <- c(a4, mean(abs(Pr2 - sig)^2))
  a5 <- c(a5, mean(abs(Pr3 - sig)^2))
  
  a1.Re <- c(a1.Re, mean(abs(r.Re - Re(sig))^2))
  a2.Re <- c(a2.Re, mean(abs(Pr0.Re - Re(sig))^2))
  a3.Re <- c(a3.Re, mean(abs(Pr1.Re - Re(sig))^2))
  a4.Re <- c(a4.Re, mean(abs(Pr2.Re - Re(sig))^2))
  a5.Re <- c(a5.Re, mean(abs(Pr3.Re - Re(sig))^2))
  
  a1.Im <- c(a1.Im, mean(abs(r.Im - Im(sig))^2))
  a2.Im <- c(a2.Im, mean(abs(Pr0.Im - Im(sig))^2))
  a3.Im <- c(a3.Im, mean(abs(Pr1.Im - Im(sig))^2))
  a4.Im <- c(a4.Im, mean(abs(Pr2.Im - Im(sig))^2))
  a5.Im <- c(a5.Im, mean(abs(Pr3.Im - Im(sig))^2))
}

mse <- c(mean(a1), mean(a2), mean(a3), mean(a4), mean(a5))
mse.Re <- c(mean(a1.Re), mean(a2.Re), mean(a3.Re), mean(a4.Re), mean(a5.Re))
mse.Im <- c(mean(a1.Im), mean(a2.Im), mean(a3.Im), mean(a4.Im), mean(a5.Im))

mse
mse.Re
mse.Im



N<-50
L <- 20
repeats <- 100
sig <- (cos(2*pi*(1:N)/10) + 1i*cos(2*pi*(1:N)/10 + pi/2))
crnk <- 1
rnk<-2

set.seed(7)

a <- vector()
for(i in 1:repeats) {
  sig.outl<-sig
  outlier.seq<-sample(1:(N),0.05 * N)
  sig.outl[outlier.seq]<-sig.outl[outlier.seq] + 5*sig.outl[outlier.seq] 
  noise <- 0.1 * rcnorm(N)
  
  ser<-sig.outl
  s <- ssa(ser, L = L, kind = "cssa")
  r <- reconstruct(s, groups = list(1:crnk))[[1]]
  a <- c(a, mean(abs(r - sig)^2))
}

mean(a)