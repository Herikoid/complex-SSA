# Define the functions D1, D2, and D3
D1 <- function(alpha, lambda) {
  term1 <- (lambda^2 * (alpha + 1))
  #term2 <- (-2 * lambda * alpha * (1 + alpha)^2)
  term2 <- (-2 * lambda * (1 + alpha)^2)
  term3 <- (4 * alpha * (-3 * alpha + 3 + 2 * alpha^2))
  
  return((1 / (12 * alpha^2 * (1 - alpha)^2)) * (term1 + term2 + term3))
  #return((1 / (2 * lambda * alpha^2 * (1 - alpha)^2)) * (term1 + term2 + term3))
  #lam <- 2*alpha-lambda
  #return((1 / (1.5 * (1 - alpha)^2)) * (6*lam - 6*alpha*lam + 4*alpha^2*lam+2*alpha^2+3-3*alpha))
  lam <- lambda
#  return(1 / (3*(1 - alpha)^2 * alpha^3)*(alpha*(lam/2)^2 + alpha^2*(lam/2)^2 -     alpha*lam/2 - alpha^3*lam/2 - 2*alpha^2*lam/2 + 3*alpha^2 - 3*alpha^3+2*alpha^4))
  #return(1 / (12*(1 - alpha)^2 * alpha^2)*((lam)^2*(1+ alpha) - 2*lam*(1+alpha)^2 + 4*alpha*(3 - 3*alpha+2*alpha^2)))
}

D2 <- function(alpha, lambda) {
  term1 <- (lambda^4)
  term2 <- (2 * lambda^3 * (3 * alpha - 2 - 3 * alpha^2))
  term3 <- (2 * lambda^2 * (3 - 9 * alpha + 12 * alpha^2 - 4 * alpha^3))
  term4 <- (4 * lambda * (4 * alpha^4 - 4 * alpha^3 - 3 * alpha^2 + 4 * alpha - 1))
  constant_term <- (8 * alpha - 56 * alpha^2 + 144 * alpha^3 - 160 * alpha^4 + 64 * alpha^5)
  
term <- 16*alpha*lambda-8*lambda^2*alpha^3-12*alpha^2*lambda+16*alpha^4*lambda+6*lambda^3*alpha+8*alpha-4*lambda-56*alpha^2+6*lambda^2+lambda^4-16*alpha^3*lambda-4*lambda^3+24*lambda^2*alpha^2-18*lambda^2*alpha-6*lambda^3*alpha^2+144*alpha^3+64*alpha^5-160*alpha^4

  return((1 / (6 * alpha^2 * lambda^2 * (alpha - 1)^2)) * 
           (term1 + term2 + term3 + term4 + constant_term))
#return((1 / (6 * alpha^2 * lambda^2 * (alpha - 1)^2)) * term)
}

D3 <- function(alpha, lambda) {
  return(2/(3*alpha))
}

# Main function to calculate Df based on lambda
calculate_Df <- function(alpha, lambda, sigma1, sigma2) {
  if (lambda >= 0 && lambda <= min(2 * (1 - 2 * alpha), 2 * alpha)) {
    D_value <- D1(alpha, lambda)
  } else if (lambda > 2 * (1 - 2 * alpha) && lambda < 2 * alpha) {
    D_value <- D2(alpha, lambda) 
    #print(D_value)
  } else if (lambda >= 2 * alpha && lambda <= 1) {
    D_value <- D3(alpha, lambda)
    #print(D_value)
  } else {
    stop("Lambda is out of bounds.")
  }
  
  # Calculate the final result
  result <- (sigma1^2 + sigma2^2) / N * D_value
  return(result)
}

# Example usage
alpha <- 0.3
sigma1 <- 1
sigma2 <- 1
N <- 3999

i <- 1
lam.set <- seq(0,2,0.01)
Df_result <- numeric(length(lam.set))
for(lambda in lam.set){
#  print(c(lambda*N/2, N*alpha, N*(1-alpha) - N*alpha))
  if(lambda <= 1) Df_result[i] <- calculate_Df(alpha, lambda, sigma1, sigma2)
  else Df_result[i] <- calculate_Df(alpha, 2-lambda, sigma1, sigma2)
  i <- i+1
}

plot(Df_result~lam.set, type = "l")

lambda <- 2*alpha
print(D1(alpha,lambda)*(sigma1^2 + sigma2^2) / N)
print(D3(alpha,1)*(sigma1^2 + sigma2^2) / N)



