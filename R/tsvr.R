#' Twin Support Vector Regression
#'
#' \code{tsvr} is an R implementation of TSVR.
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and response variable.
#' @param C1,C2 plenty term.
#' @param kernel kernel function.
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param reg regularization term to take care of problems due to ill-conditioning in dual problem.
#' @param kernel_rect set kernel size. \code{0<= kernel_rect <= 1}.
#' @param tol the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param rcpp speed up your code with Rcpp, default \code{rcpp = TRUE}.
#' @references Xinjun Peng. TSVR: An efficient Twin Support Vector Machine for
#' regression, Neural Networks, Volume 23, Issue 3, 2010, Pages 365-372.
#' @return return tsvr object.
#' @export

tsvr <- function(X, y, C1 = 1.0, C2 = 1.0,
                    kernel = c('linear', 'rbf', 'poly'),
                    reg = 1e-5, kernel_rect = 1,
                    gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                    tol = 1e-5, max.steps  = 300,
                    rcpp = TRUE){

  kernel <- match.arg(kernel)

  X <- as.matrix(X)
  y <- as.matrix(y)

  n <- nrow(X)
  m <- ncol(X)

  e <- matrix(1, nrow = n)
  if(kernel != 'linear'){
    kernel_m <- round(n*kernel_rect, 0)
    X <- kernel_function(X, as.matrix(X[1:kernel_m, ]),
                         kernel.type = kernel,
                         gamma = gamma, degree = degree, coef0 = coef0,
                         rcpp = rcpp)
  }
  G <- cbind(X, e)
  f <- y - e*C1
  h <- y + e*C2
  H <- G %*% solve(t(G) %*% G + diag(rep(reg, ncol(G)))) %*% t(G)
  q1 <- t(t(f) - t(f) %*% H)
  q2 <- t(t(h) %*% H - t(h))
  lb <- matrix(0, nrow = n)
  ub1 <- matrix(C1, nrow = n)
  ub2 <- matrix(C2, nrow = n)
  alphas <- clip_dcd_optimizer(H, q1, lb, ub1, tol, max.steps, rcpp = rcpp)$x
  gammas <- clip_dcd_optimizer(H, q2, lb, ub2, tol, max.steps, rcpp = rcpp)$x
  u1 <- solve(t(G) %*% G + diag(rep(reg, ncol(G)))) %*% t(G) %*% (f - alphas)
  u2 <- solve(t(G) %*% G + diag(rep(reg, ncol(G)))) %*% t(G) %*% (h + gammas)
  coef <- (u1 + u2) / 2
  fitted <- G %*% coef
  tsvr_model <- list("coef" = coef,
                        "coef_1" = u1,
                        "coef_2" = u2,
                        "fitted" = fitted,
                        "C1" = C1,
                        "C2" = C2,
                        "kernel_rect" = kernel_rect,
                        "kernel" = kernel,
                        "gamma" = gamma,
                        "degree" = degree,
                        "coef0" = coef0,
                        "Rcpp" = rcpp,
                        "X" = X,
                        "y" = y
                        )
  class(tsvr_model) <- "tsvr"
  return(tsvr_model)
}


#' Predict Method for Multiple Birth Support Vector Machines
#'
#' @author Zhang Jiaqi
#' @param object a fitted object of class inheriting from \code{twinsvr}.
#' @param X a new data frame for predicting.
#' @param y response variable data frame corresponding to X.
#' @param ... unused parameter.
#' @importFrom stats predict
#' @export

predict.tsvr <- function(object, X, y, ...){
  km <- nrow(object$X)
  if(object$kernel != 'linear'){
    kernel_m <- round(km*object$kernel_rect, 0)
    X <- kernel_function(X, as.matrix(X[1:kernel_m, ]),
                         kernel.type = object$kernel,
                         gamma = object$gamma, degree = object$degree,
                         coef0 = object$coef0,
                         rcpp = object$Rcpp)
  }

  kernelX <- as.matrix(cbind(X, 1))
  y_hat <- kernelX %*% object$coef
  mse <- mean_squared_error(y, y_hat)
  return(y_hat)
}
