ls_svm_dual_solver <- function(KernelX, y, C = 1, rcpp = TRUE) {
  D <- diag(as.vector(y))
  H <- D %*% KernelX %*% D
  m <- nrow(KernelX)
  if (rcpp == TRUE) {
    alphas <- cpp_chol_solve(H + diag(1/C, m), matrix(1, nrow = m))
  } else {
    alphas <- chol2inv(chol((H + diag(1/C, m)))) %*% matrix(1, nrow = m)
  }
  coef <- D %*% alphas
  BaseDualLeastSquaresSVMClassifier <- list(coef = as.matrix(coef))
  class(BaseDualLeastSquaresSVMClassifier) <- "BaseDualLeastSquaresSVMClassifier"
  return(BaseDualLeastSquaresSVMClassifier)
}


ls_svm_primal_solver <- function(KernelX, y, C = 1, eps = 1e-5,
                                 max.steps = 80, batch_size = nrow(KernelX) / 10,
                                 seed = NULL, sample_seed = NULL,
                                 optimizer = pegasos, ...) {
   gLeastSquares <- function(KernelX, y, v, ...) { # gradient of Least Squares loss function
    C <- list(...)$C
    xn <- nrow(KernelX)
    xp <- ncol(KernelX)
    sg <- matrix(0, nrow = xp, ncol = 1)
    u <- 1 - y * (KernelX %*% v)
    sg <- v - (C/xn) * t(KernelX) %*% (u*y)
    return(sg)
  }
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  w0 <- matrix(0, nrow = xp, ncol = 1)
  wt <- optimizer(KernelX, y, w0, batch_size, max.steps, gLeastSquares, sample_seed, C = C, ...)
  BasePrimalLeastSquaresSVMClassifier <- list(coef = as.matrix(wt[1:xp]))
  class(BasePrimalLeastSquaresSVMClassifier) <- "BasePrimalLeastSquaresSVMClassifier"
  return(BasePrimalLeastSquaresSVMClassifier)
}


#' Least Squares Support Vector Machine
#'
#' \code{ls_svm} is an R implementation of Hinge-SVM
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
#' @return return \code{SVMClassifier} object.
#' @export
ls_svm <- function(X, y, C = 1, kernel = c("linear", "rbf", "poly"),
                   gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                   eps = 1e-5, max.steps = 80, batch_size = nrow(X) / 10,
                   solver = c("dual", "primal"), rcpp = TRUE,
                   fit_intercept = TRUE, optimizer = pegasos, randx = 0.1, ...) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  class_set <- unique(y)
  idx <- which(y == (sort(class_set)[1]))
  y[idx] <- 1
  y[-idx] <- -1
  y <- as.matrix(as.numeric(y))
  if (length(class_set) > 2) {
    stop("The number of class should less 2!")
  }
  kernel <- match.arg(kernel)
  solver <- match.arg(solver)
  if (fit_intercept == TRUE) {
    X <- cbind(X, 1)
  }
  if (kernel == "linear" & solver == "primal") {
    KernelX <- X
  } else if (kernel != "linear" & solver == "primal") {
    if (randx > 0) {
      randX = X[sample(nrow(X), floor(randx*nrow(X))),]
    } else {
      randX <- X
    }
    KernelX <- kernel_function(X, randX,
                               kernel.type = kernel,
                               gamma = gamma, degree = degree, coef0 = coef0,
                               rcpp = rcpp)
    X <- randX
  } else if (solver == "dual") {
    KernelX <- kernel_function(X, X,
                               kernel.type = kernel,
                               gamma = gamma, degree = degree, coef0 = coef0,
                               rcpp = rcpp)
  }
  if (solver == "primal") {
    solver.res <- ls_svm_primal_solver(KernelX, y, C, eps,
                                       max.steps, batch_size,
                                       optimizer, ...)
  } else if (solver == "dual") {
    solver.res <- ls_svm_dual_solver(KernelX, y, C, rcpp)
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
