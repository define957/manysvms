linear_kernel <- function(x1, x2){
  return(x1 %*% x2)
}

rbf_kernel <- function(x1, x2, gamma= 1){
  temp <- as.matrix(x1 - x2)
  K <-  exp(- gamma * norm(temp, type = '2')^2)
  return(K)
}
