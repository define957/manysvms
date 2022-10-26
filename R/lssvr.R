#' Least Squares  Support Vector Regression
#'
#' \code{lssvr} is an R implementation of LSSVR
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and response variable.
#' @param C plenty term.
#' @param kernel kernel function.
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param rcpp speed up your code with Rcpp, default \code{rcpp = TRUE}.
#' @return return lssvr object.
#' @export
#'
lssvr <- function(X, y, C = 1.0,
                  kernel = c('linear', 'rbf', 'poly'),
                  gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                  rcpp = TRUE){
  kernel <- match.arg(kernel)

  X <- as.matrix(X)
  y <- as.matrix(y)

  n <- nrow(X)
  m <- ncol(X)

  kernel_X <- kernel_function(X, X,
                       kernel.type = kernel,
                       gamma = gamma, degree = degree, coef0 = coef0,
                       rcpp = rcpp)
  kernel_X <- kernel_X + diag(1, n)/C
  e <- matrix(1, nrow = n)
  mat1 <- cbind(kernel_X, e)
  mat2 <- cbind(t(e), 0)
  mat <- rbind(mat1, mat2)
  coef <- solve(mat) %*% rbind(y, 0)
  fitted.values <- kernel_X %*% coef[1:n, 1] + coef[n, 1]
  if(kernel == "linear"){
    coef = rbind(t(X) %*% coef[1:n ], coef[n+1, 1])
  }
  lssvr_model <- list("X" = X, "y" = y,
                      "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                      "coef" = coef,
                      "fitted.values" = fitted.values)
  return(lssvr_model)
}
