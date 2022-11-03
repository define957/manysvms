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
#' x2 <- sin(x1) + mvrnorm(obs, mu = 0, Sigma = 0.03)
#'
#' X <- as.matrix(x1)
#' y <- as.matrix(x2)
#'
#' s <- Sys.time()
#' m <- eps.svr(X, y, eps=0.3, kernel = 'rbf', C = 1,
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

  Q <- kernel_function(X, X,
                       kernel.type = kernel,
                       gamma = gamma, degree = degree, coef0 = coef0,
                       rcpp = rcpp)


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
  fitted <- matrix(coef %*% Q, nrow = m)
  svr <- list('X' = X, 'y' = y,
              'coef' = coef,
              'epsilon' = eps,
              'fitted' = fitted,
              'kernel' = kernel,
              'gamma' = gamma,
              'call' = match.call(),
              "Rcpp" = rcpp)
  class(svr) <- 'eps.svr'

  return(svr)
}


#' Print Method for epsilon - Support Vector Regression
#'
#' This function prints information about \code{eps.svr} model.
#' @param x object of class \code{twinsvm}.
#' @param ... unsed argument.
#' @export

print.eps.svr <- function(x, ...){
  cat("\nCall:", deparse(x$call, 0.8 * getOption("width")), "\n", sep="\n")
  cat("SVM type : ", class(x), "\n")
  cat("SVM kernel : ", x$kernel, "\n")
  if(x$kernel == 'rbf'){
    cat("gamma : ", x$gamma, "\n")
  }
  cat("number of observations : ", nrow(x$X), "\n")
}

cv.svr <- function(X, y , K = 5,
                     eps = 0.1,
                     kernel = c('linear', 'rbf', 'poly'),
                     shuffer = TRUE, seed = NULL){

  m <- nrow(X)
  if(shuffer == TRUE){
    if(is.null(seed) == FALSE){
      set.seed(seed)
    }
    new_idx <- sample(m)

  }else{
    new_idx <- 1:m
  }
  v_size <- m %/% K
  indx_cv <- 1
  mse_list <- c()
  for(i in 1:K){
    new_idx_k <- new_idx[indx_cv:(indx_cv+v_size - 1)] #get test dataset
    indx_cv <- indx_cv + v_size
    test_X <- X[new_idx_k, ]
    train_X <- X[-new_idx_k, ]
    test_y <- y[new_idx_k]
    train_y <- y[-new_idx_k]
    jssvr_model <- svr(train_X, train_y, eps = eps,
                         max.steps = 800, kernel = 'rbf')
    pred <- predict(jssvr_model, test_X, test_y)
    mse <- mean_squared_error(test_y, pred)
    mse_list <- append(mse_list, mse)
    cat('MSE in ',K, 'fold cross validation :', mse, '\n')
  }
  cat('average MSE in ',K, 'fold cross validation :', mean(mse_list), '\n')
  cat('Sd of MSE in ',K, 'fold cross validation :', sd(mse_list), '\n')
}


predict.eps.svr <- function(object, X, y = NULL, ...){
  m <-  nrow(X)
  X <- X %*% object$Projection
  Q <- kernel_function(X, object$X%*%object$Projection,
                       kernel.type = object$kernel,
                       gamma = object$gamma, degree = object$degree, coef0 = object$coef0,
                       rcpp = object$Rcpp)
  pred <- matrix(Q %*% object$coef , nrow = m)
  return(pred)
}

#' Predict Method for epsilon - Support Vector Regression
#'
#' @author Zhang Jiaqi
#' @param object A fitted object of class inheriting from \code{eps.svr}.
#' @param X A new data frame for predicting.
#' @param y A label data frame corresponding to X.
#' @param ... unused parameter.
#' @importFrom stats predict
#' @export

predict.eps.svr <- function(object, X, y, ...){
  X <- as.matrix(X)
  y <- as.matrix(y)
  X <- kernel_function(X, object$X,
                         kernel.type = object$kernel,
                         gamma = object$gamma, degree = object$degree,
                         coef0 = object$coef0,
                         rcpp = object$Rcpp)
  y_hat <- X %*% object$coef
  return(y_hat)
}
