library(QZ, quiet = TRUE)

H.1 <- function(E, U, V) {
  Pcol <- U[[1]] %*% H(U[[1]])
  for(i in 2:length(U)) {
    Pcol <- Pcol + U[[i]] %*% H(U[[i]])
  }
  Prow <- V[[1]] %*% H(V[[1]]) 
  for(i in 2:length(V)) {
    Prow <- Prow + V[[i]] %*% H(V[[i]])
  }
  Pcol.orth <- diag(rep(1,length(U[[1]]))) - Pcol
  Pcol %*% E + Pcol.orth %*% E %*% Prow
}
