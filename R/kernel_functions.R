linear_kernel <- function(x1, x2){
  return(x1 %*% x2)
}

rbf_kernel <- function(x1, x2, gamma= 1){
  temp <- as.matrix(x1 - x2)
  K <-  exp(- gamma * norm(temp, type = '2')^2)
  return(K)
}

poly_kernel <- function(x1, x2, gamma, degree = 3, coef0 = 0){
  K <- (gamma * x1 %*% t(x2) + coef0)^(degree)
  return(K)
}

r_rbf_kernel <- function(x1, x2, gamma= 1){
  m <- nrow(x1)
  n <- nrow(x2)
  K <- matrix(0, nrow = m, ncol = n)
  for(i in 1:m){
    for(j in 1:n){
      temp <- as.matrix(x1[i, ] - x2[j, ])
      K[i, j] <- exp(- gamma * norm(temp, type = '2')^2)
    }
  }
  return(K)
}
