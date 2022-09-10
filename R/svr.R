#' epsilon - Support Vector Regression
#'
#' \code{eps.svr} is an R implementation of epsilon - support vector regression
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and explained variable.
#' @param eps epsilon in the insensitive-loss function (default \code{eps = 0.1}).
#' @param kernel kernel function.
#' @param C plenty term (default \code{C = 1}).
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param rcpp speed up your code with Rcpp, default \code{rcpp = TRUE}.
#' @return return eps.svr object.
#' @useDynLib manysvms, .registration = TRUE
#' @import Rcpp
#' @export
#' @examples
#' library(manysvms)
#' library(MASS)
#' library(ggplot2)
#'
#' set.seed(1235)
#' x1 <- mvrnorm(200, mu = 0, Sigma = 3)
#' x2 <- 2*x1  + mvrnorm(200, mu = 0, Sigma = 1)
#'
#' X <- as.matrix(x1)
#' y <- as.matrix(x2)
#'
#' m <- eps.svr(X,y, eps=2, C = 1)
#' dataXy <- as.data.frame(cbind(X, y))
#' ggplot(data = dataXy, aes(x = X, y = y))+
#'   geom_point()+
#'   geom_abline(slope = m$coef, intercept = m$intercept)+
#'   geom_abline(slope = m$coef, intercept = m$intercept + m$epsilon, linetype=2, color = 'red')+
#'   geom_abline(slope = m$coef, intercept = m$intercept - m$epsilon, linetype=2, color = 'red')+
#'   theme_classic()

eps.svr <- function(X, y, eps = 0.1,
                    kernel = c('linear', 'rbf', 'poly'),
                    C = 1, gamma = 1 / ncol(X),
                    max.steps = 1000, rcpp = TRUE){
  X <- as.matrix(X)
  y <- as.matrix(y)

  kernel <- match.arg(kernel)
  m <-  nrow(X)
  if(kernel == 'linear'){
    Q <- X%*%t(X)
  }else if(kernel == 'rbf'){
    if(rcpp == FALSE){
      Q <- matrix(0, nrow = m, ncol = m)
      for(i in 1:m){
        for(j in 1:m){
          Q[i, j] <-  rbf_kernel(X[i, ], X[j, ], gamma = gamma)
        }
      }
    }else{
      Q <- cpp_rbf_kernel(X, X, gamma = gamma)
    }

  }

  Q1 <- cbind(Q, -Q)
  Q2 <- cbind(-Q, Q)
  H <- rbind(Q1, Q2)

  e <- rep(1, nrow(y))

  q1 <- eps * e - y
  q2 <- eps * e + y

  q <- t(rbind(c(q1, q2)))

  lb <- matrix(0, nrow = nrow(q))
  ub <- matrix(C, nrow = nrow(q))
  if(rcpp == TRUE){
    beta <- cpp_clip_dcd_optimizer(H, -q, lb, ub, 1e-12, max.steps)
  }else{
    beta <- clip_dcd_optimizer(H, -q, lb, ub, max.steps = max.steps)$x
  }
  w <- (beta[1:nrow(y)] - beta[-c(1:nrow(y))]) %*% X
  b <- mean(y - X%*%w *  + eps)
  fitted <- (beta[1:nrow(y)] - beta[-c(1:nrow(y))]) %*% Q
  svr <- list('coef' = w, 'intercept' = b, 'epsilon' = eps, 'fitted' = fitted)
  class(svr) <- 'eps.svr'
  return(svr)
}
