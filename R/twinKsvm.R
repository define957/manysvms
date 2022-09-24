#' Twin K SVM for Multi-classification by Using Ones versus Rest Strategy
#'
#' @author Zhang Jiaqi
#' @param X,y dataset and label
#' @param Ck plenty term vector
#' @param kernel kernel function
#' @param gamma rbf kernel parameter
#' @param  reg regularization tern
#' @param kernel_rect set kernel size. \code{0<= kernel_rect <= 1}
#' @param tol the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param rcpp speed up your code with Rcpp, default \code{rcpp = TRUE}.
#' @return return twinKsvm object
#' @export
twinKsvm <- function(X, y,
                     Ck = rep(1, length(unique(y))),
                     kernel = c('linear', 'rbf', 'poly'),
                     gamma = 1 / ncol(X), reg = 1, kernel_rect = 1,
                     tol = 1e-5, max.steps = 300, rcpp = TRUE){
  kernel <- match.arg(kernel)

  X <- as.matrix(X)
  y <- as.matrix(y)

  m <- nrow(X)
  n <- ncol(X)

  class_set <- unique(as.matrix(y))
  class_num <- length(class_set)

  if(kernel == 'linear'){
    coef_dim <- n
  }else if(kernel == 'rbf'){
    coef_dim <- m*kernel_rect
  }

  coef_list <- matrix(0, nrow = coef_dim, ncol = class_num)
  intercept_list <- matrix(0, ncol = class_num)

  # solve K * K-1 models
  for(i in 1:class_num){
    for(j in i+1:class_num){

    }
  }
}
