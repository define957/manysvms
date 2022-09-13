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
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param tol tolerance of termination criterion, default: \code{1e-5}.
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
#' obs <- 200
#' x1 <- seq(-5, 5, length.out = obs)
#' x2 <- sin(x1)  + mvrnorm(obs, mu = 0, Sigma = 0.03)
#'
#' X <- as.matrix(x1)
#' y <- as.matrix(x2)
#'
#' s <- Sys.time()
#' m <- eps.svr(X,y, eps=0.3, kernel = 'rbf', C = 1,
#'              gamma = 1, max.steps = 500, rcpp = TRUE)
#' e <- Sys.time()
#' print(e - s)
#' dataXy <- as.data.frame(cbind(X, y))
#'
#' ggplot(data = dataXy, aes(x = X, y = y))+
#'   geom_point()+
#'   geom_line(aes(x=X, y=m$fitted))+
#'   geom_line(aes(x=X, y=m$fitted + m$epsilon, color = 'red'), show.legend = FALSE)+
#'   geom_line(aes(x=X, y=m$fitted - m$epsilon, color = 'red'), show.legend = FALSE)+
#'   theme_bw()

eps.svr <- function(X, y, eps = 0.1,
                    kernel = c('linear', 'rbf', 'poly'),
                    C = 1, gamma = 1 / ncol(X), degree = 3,
                    coef0 = 0, max.steps = 1000, tol = 1e-5, rcpp = TRUE){
  X <- as.matrix(X)
  y <- as.matrix(y)

  kernel <- match.arg(kernel)
  m <-  nrow(X)
  if(kernel == 'linear'){
    Q <- X%*%t(X)
  }
  if(rcpp == FALSE){
    if(kernel == 'rbf'){
      Q <- r_rbf_kernel(X, X, gamma = gamma)
    }else if(kernel == 'poly'){
      Q <-  poly_kernel(X, X,
                        gamma = gamma, degree = degree,
                        coef0 = coef0)
    }
  }else if(rcpp == TRUE){
    if(kernel == 'rbf'){
      Q <- cpp_rbf_kernel(X, X, gamma = gamma)
    }else if(kernel == 'poly'){
      Q <- cpp_poly_kernel(X, X, gamma = gamma,
                           degree = degree, coef0 = coef0)
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

  beta <- clip_dcd_optimizer(H, -q, lb, ub, eps = tol, max.steps, rcpp = rcpp)$x
  coef <- (beta[1:nrow(y)] - beta[-c(1:nrow(y))])
  fitted <- coef %*% Q
  svr <- list('coef' = coef, 'epsilon' = eps, 'fitted' = fitted)
  class(svr) <- 'eps.svr'

  return(svr)
}
