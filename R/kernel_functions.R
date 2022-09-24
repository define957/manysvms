r_linear_kernel <- function(x1, x2){
  return(x1 %*% x2)
}

r_poly_kernel <- function(x1, x2, gamma, degree = 3, coef0 = 0){
  K <- (gamma * x1 %*% t(x2) + coef0)^(degree)
  return(K)
}

r_rbf_kernel <- function(x1, x2, gamma= 1/ncol(x2)){
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

rbf_kernel <- function(x1, x2, gamma = 1/ncol(x1), rcpp = TRUE){
  if(rcpp == TRUE){
    K <- cpp_rbf_kernel(x1, x2, gamma)
  }else if(rcpp == FALSE){
    K <- r_rbf_kernel(x1, x2, gamma)
  }
  return(K)
}


kernel_function <- function(x1, x2,
                            kernel.type = c('linear', 'rbf', 'poly'),
                            gamma = 1/ncol(x1), degree = 3, coef0 = 0,
                            rcpp = TRUE){
  kernel.type <- match.arg(kernel.type)
  if(kernel.type == 'linear'){
    K <- r_linear_kernel(x1, x2)
  }else if(kernel.type == 'rbf'){
    K <- rbf_kernel(x1, x2, gamma, rcpp)
  }else if(kernel.type == 'poly'){
    K <- r_poly_kernel(x1, x2, gamma, degree = 3, coef0 = 0)
  }
  return(K)
}
