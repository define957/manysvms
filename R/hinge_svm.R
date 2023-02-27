hinge_svm_dual_solver <- function (KernelX, y, C = 1,
                                   eps = 1e-5, max.steps = 80, rcpp = TRUE) {
  D <- diag(as.vector(y))
  n <- nrow(KernelX)
  H <- D %*% KernelX %*% D
  e <- matrix(1, nrow = n)
  lb <- matrix(0, nrow = n)
  ub <- matrix(C, nrow = n)
  
  alphas <- clip_dcd_optimizer(H, e, lb, ub, eps, max.steps, rcpp)$x
  coef <- D %*% alphas
  BaseDualHingeSVMClassifier <- list(coef = as.matrix(coef))
  class(BaseDualHingeSVMClassifier) <- "BaseDualHingeSVMClassifier"
  return(BaseDualHingeSVMClassifier)
}


hinge_svm_primal_solver <- function (KernelX, y, C = 1, eps = 1e-5,
                                     max.steps = 80, batch_size = nrow(KernelX) / 10,
                                     seed = NULL, sample_seed = NULL,
                                     optimizer = pegasos, ...) {
  sgHinge <- function(KernelX, y, v, ...) { # sub-gradient of hinge loss function
    C <- list(...)$C
    xn <- nrow(X)
    xp <- ncol(X)
    sg <- matrix(0, nrow = xp, ncol = 1)
    u <- 1 - y * (X %*% v)
    u[u <= 0] <- 0
    u[u > 0] <- 1
    sg <- v - (C/xn) * t(X)%*%(u*y)
    return(sg)
  }
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  w0 <- matrix(0, nrow = xp, ncol = 1)
  wt <- optimizer(KernelX, y, w0, batch_size, max.steps, sgHinge, sample_seed, C = C, ...)
  wnorm <- norm(wt[1:xp], type = "2")
  BasePrimalHingeSVMClassifier <- list(coef = as.matrix(wt[1:xp]))
  class(BasePrimalHingeSVMClassifier) <- "BasePrimalHingeSVMClassifier"
  return(BasePrimalHingeSVMClassifier)
}


#' Hinge Support Vector Machine
#'
#' \code{hinge_svm} is an R implementation of Hinge-SVM
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param C plenty term.
#' @param kernel kernel function. The definitions of various kernel functions are as follows:
#' \describe{
#'     \item{linear:}{\eqn{u'v}{u'*v}}
#'     \item{poly:}{\eqn{(\gamma u'v + coef0)^{degree}}{(gamma*u'*v + coef0)^degree}}
#'     \item{rbf:}{\eqn{e^{(-\gamma |u-v|^2)}}{exp(-gamma*|u-v|^2)}}
#' }
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param eps the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param batch_size mini-batch size for primal solver.
#' @param solver \code{"dual"} and \code{"primal"} are available.
#' @param rcpp speed up your code with Rcpp, default \code{rcpp = TRUE}.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param optimizer default primal optimizer pegasos.
#' @param randx parameter for reduce SVM, default \code{randx = 0.1}.
#' @param ... unused parameters.
#' @return return \code{HingeSVMClassifier} object.
#' @export
hinge_svm <- function (X, y, C = 1, kernel = c("linear", "rbf", "poly"),
                       gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                       eps = 1e-5, max.steps = 80, batch_size = nrow(X) / 10,
                       solver = c("dual", "primal"), rcpp = TRUE,
                       fit_intercept = TRUE, optimizer = pegasos, randx = 0.1, ...) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  class_set <- unique(y)
  idx <- which(y == class_set[1])
  y[idx] <- 1
  y[-idx] <- -1
  if (length(class_set) > 2) {
    cat(errorCondition("Error: The number of class should less 2!"))
  }
  kernel <- match.arg(kernel)
  solver <- match.arg(solver)
  if (fit_intercept == TRUE) {
    X <- cbind(X, 1)
  }
  if (kernel == "linear" & solver == "primal") {
    KernelX <- X
  } else if (kernel != "linear" & solver == "primal"){
    if(randx > 0 ){
      randX = x[sample(nrow(X), floor(randx*nrow(X))),] 
    }
    KernelX <- kernel_function(X, randX,
                               kernel.type = kernel,
                               gamma = gamma, degree = degree, coef0 = coef0,
                               rcpp = rcpp)
    X <- randx
  } else if (solver == "dual"){
    KernelX <- kernel_function(X, X,
                               kernel.type = kernel,
                               gamma = gamma, degree = degree, coef0 = coef0,
                               rcpp = rcpp)
  }
  if (solver == "primal") {
    solver.res <- hinge_svm_primal_solver(KernelX, y, C, eps,
                                          max.steps, batch_size,
                                          optimizer, ...)
  } else if(solver == "dual") {
    solver.res <- hinge_svm_dual_solver(KernelX, y, C, eps,
                                        max.steps, rcpp)
  }
  SVMClassifier <- list("X" = X, "y" = y, "class_set" = class_set,
                        "C" = C, "kernel" = kernel,
                        "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                        "solver" = solver, "coef" = solver.res$coef,
                        "fit_intercept" = fit_intercept,
                        "rcpp" = rcpp)
  class(SVMClassifier) <- "SVMClassifier"
  return(SVMClassifier)
}


#' Predict Method for Hinge Support Vector Machine
#'
#' @author Zhang Jiaqi
#' @param object a fitted object of class inheriting from \code{HingeSVMClassifier}.
#' @param X new data for predicting.
#' @param ... unused parameter.
#' @importFrom stats predict
#' @export
predict.SVMClassifier <- function(object, X, ...) {
  if (object$fit_intercept == TRUE) {
    X <- cbind(X, 1)
  }
  if (object$kernel == "linear" & object$solver == "primal") {
    KernelX <- X
  } else {
    KernelX <- kernel_function(X, object$X,
                               kernel.type = object$kernel,
                               gamma = object$gamma,
                               degree = object$degree,
                               coef0 = object$coef0,
                               rcpp = object$rcpp)
  }
  fx <- KernelX %*% object$coef
  decf <- sign(fx)
  decf[decf > 0] <- object$class_set[1]
  decf[decf < 0] <- object$class_set[2]
  return(decf)
}




