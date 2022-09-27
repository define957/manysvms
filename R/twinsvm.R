#' Twin Support Vector Machines
#'
#' \code{twinsvm} is an R implementation of Twin SVM.
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param C1 plenty term 1.
#' @param C2 plenty term 2.
#' @param kernel kernel function. The definitions of various kernel functions are as follows:
#' \describe{
#'     \item{linear:}{\eqn{u'v}{u'*v}}
#'     \item{poly:}{\eqn{(\gamma u'v + coef0)^{degree}}{(gamma*u'*v + coef0)^degree}}
#'     \item{rbf:}{\eqn{e^{(-\gamma |u-v|^2)}}{exp(-gamma*|u-v|^2)}}
#' }
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param reg regularization term to take care of problems due to ill-conditioning in dual problem.
#' @param kernel_rect set kernel size. \code{0<= kernel_rect <= 1}.
#' @param tol the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param rcpp speed up your code with Rcpp, default \code{rcpp = TRUE}.
#' @return return twinsvm object.
#' @export
#' @examples
#' library(manysvms)
#' data('iris')
#'
#' m <- 100
#' X <- iris[0:m, 1:2]
#' y <- iris[0:m, 5]
#' model <- twinsvm(X, y, kernel = 'linear')
#'
#' coef(model)
#' pred <-predict(model, X, y)

twinsvm <- function(X, y,
                    C1 = 1.0, C2 = 1.0,
                    kernel = c('linear', 'rbf', 'poly'),
                    gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                    reg = 1, kernel_rect = 1,
                    tol = 1e-5, max.steps = 300, rcpp = TRUE){

  kernel <- match.arg(kernel)

  X <- as.matrix(X)
  y <- as.matrix(y)

  m <- nrow(X)
  n <- ncol(X)
  class_set <- unique(y)
  num_class <- length(class_set)
  if(num_class > 2){
    return(0)
  }

  idxA <- which(y == class_set[1])
  idxB <- which(y == class_set[2])

  A <- X[idxA, ]
  B <- X[idxB, ]

  mA <- nrow(A)
  mB <- nrow(B)

  e1 <- matrix(1, nrow = mA)
  e2 <- matrix(1, nrow = mB)

  if(kernel == 'linear'){
    S <- A
    R <- B

  }else{
    kernel_m <- round(m*kernel_rect, 0)

    S <- kernel_function(A, X[1:kernel_m, ],
                         kernel.type = kernel,
                         gamma = gamma, degree = degree, coef0 = coef0,
                         rcpp = rcpp)
    R <- kernel_function(B, X[1:kernel_m, ],
                         kernel.type = kernel,
                         gamma = gamma, degree = degree, coef0 = coef0,
                         rcpp = rcpp)
  }
  S <- cbind(S, e1)
  R <- cbind(R, e2)

  # Solve QP 1
  STS_reg_inv <- solve(t(S) %*% S + diag(rep(reg, ncol(S))))
  H <- R %*% STS_reg_inv %*% t(R)
  lbA <- matrix(0, nrow = mB)
  ubA <- matrix(C1, nrow = mB)
  AA <- diag(rep(1, nrow(H)))
  qp1_solver <- clip_dcd_optimizer(H, e2, lbA, ubA, tol, max.steps, rcpp)
  alphas <- as.matrix(qp1_solver$x)
  Z1 <- - STS_reg_inv %*% t(R) %*% alphas

  # Solve QP 2
  RTR_reg_inv <- solve(t(R) %*% R + diag(rep(reg, ncol(R))))
  H <- S %*% RTR_reg_inv %*% t(S)
  lbB <- matrix(0, nrow = mA)
  ubB <- matrix(C2, nrow = mA)
  AB <- diag(rep(1, nrow(H)))
  qp2_solver <- clip_dcd_optimizer(H, e1, lbB, ubB, tol, max.steps, rcpp)
  gammas <- as.matrix(qp2_solver$x)
  Z2 <- RTR_reg_inv %*% t(S) %*% gammas
  u_idx <- length(Z2) - 1
  u1 <- Z1[1:u_idx]
  u2 <- Z2[1:u_idx]
  b1 <- Z1[u_idx+1]
  b2 <- Z2[u_idx+1]

  twinsvm_model <- list('X' = X,
                        'y' = y,
                        'u1' = u1,
                        'u2' = u2,
                        'b1' = b1,
                        'b2' = b2,
                        'class_set' = class_set,
                        'kernel' = kernel,
                        'gamma' = gamma,
                        'coef0' = coef0,
                        'kernel_rect' = kernel_rect,
                        'Rcpp' = rcpp,
                        'call' = match.call())
  class(twinsvm_model)<-"twinsvm"

  return(twinsvm_model)
}


#' Coef Method for Twin Support Vector Machines
#'
#' This function solve coefficients based upon a model trained by \code{twinsvm}.
#' @author Zhang Jiaqi
#' @param object object of class `twinsvm`.
#' @param ... unused parameter.
#' @importFrom stats predict
#' @export
#' @examples
#' library(manysvms)
#' data('iris')
#'
#' m <- 100
#' X <- iris[0:m, 1:2]
#' y <- iris[0:m, 5]
#' model <- twinsvm(X, y, kernel = 'linear')
#' pred <-predict(model, X, y)

coef.twinsvm <- function(object, ...){
  coefs <- list('u1' = object$u1, 'u2' = object$u2,
                'b1' = object$b1, 'b2'= object$b2)
  return(coefs)
}


#' Predict Method for Twin Support Vector Machines
#'
#' This function predicts values based upon a model trained by \code{twinsvm}.
#'
#' @author Zhang Jiaqi
#' @param object object of class `twinsvm`.
#' @param X a new data frame for predicting.
#' @param y a label data frame corresponding to X.
#' @param ... unused parameter.
#' @importFrom stats predict
#' @export
#' @examples
#' library(manysvms)
#' data('iris')
#'
#' m <- 100
#' X <- iris[0:m, 1:2]
#' y <- iris[0:m, 5]
#' model <- twinsvm(X, y, kernel = 'linear')
#' pred <-predict(model, X, y)

predict.twinsvm <- function(object, X, y, ...){

  X <- as.matrix(X)
  y <- as.matrix(y)
  if(object$kernel == 'linear'){
    kernelX <- X
  }else{
    kernel_m <- round(m*object$kernel_rect, 0)
    kernelX <- kernel_function(X, X[1:kernel_m, ],
                               kernel.type = object$kernel,
                               gamma = object$gamma,
                               degree = object$degree,
                               coef0 = object$coef0,
                               rcpp = object$Rcpp)
  }

  disA <- abs(kernelX %*% object$u1 + object$b1) / norm(object$u1, type = '2')
  disB <- abs(kernelX %*% object$u2 + object$b2) / norm(object$u2, type = '2')
  dis <- cbind(disA, disB)
  predict_idx <- as.vector(apply(dis, 1, which.min))
  conf <- as.vector(apply(dis, 1, which.min))

  m <- length(predict_idx)
  predict_value <- rep(0, m)
  for(i in 1:m){
    predict_value[i] <- object$class_set[predict_idx[i]]
  }
  accuracy <- sum(predict_value == y) * 100 / m
  predlist <- list('predict.value' = predict_value, "accuracy"= accuracy,
                   'A.distance' = disA, 'B.distance' = disB)
  cat('use kernel : ', object$kernel, '\n')
  cat('total accuracy :', accuracy, '%\n')
  return(predlist)
}

#' @rdname twinsvm
#' @format NULL
#' @usage NULL
#' @export

print.twinsvm <- function(x, ...){
  cat("\nCall:", deparse(x$call, 0.8 * getOption("width")), "\n", sep="\n")
  cat("SVM type : ", class(x), "\n")
  cat("SVM kernel : ", x$kernel, "\n")
  if(x$kernel == 'rbf'){
    cat("gamma : ", x$gamma, "\n")
  }
  cat("number of observations : ", nrow(x$X), "\n")
  cat("number of class : ", length(x$class_set), "\n")
}
