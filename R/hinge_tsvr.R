hinge_tsvr_dual_solver <- function(KernelX, y, C1, C2, epsilon1, epsilon2,
                                   eps, max.steps) {
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  H <- KernelX
  HTH_inv_H <- cholsolve(t(H) %*% H + diag(1e-7, xp), t(H))
  dualH <- H %*% HTH_inv_H

  f <- as.matrix(y - epsilon1)
  h <- as.matrix(y + epsilon2)

  q1 <- dualH %*% f - f
  q2 <- h - dualH %*% h

  lb <- matrix(0, xn, 1)
  ub1 <- matrix(C1, xn, 1)
  ub2 <- matrix(C2, xn, 1)
  x0 <- lb
  alphas <- clip_dcd_optimizer(dualH, q1, lb, ub1, eps, max.steps, x0)$x
  gammas <- clip_dcd_optimizer(dualH, q2, lb, ub2, eps, max.steps, x0)$x

  u1 <- HTH_inv_H %*% (f - alphas)
  u2 <- HTH_inv_H %*% (h + gammas)
  BaseDualHingeTSVRRegressor <- list("coef1" = as.matrix(u1),
                                     "coef2" = as.matrix(u2),
                                     "lag1" = alphas,
                                     "lag2" = gammas)
}

#' Hinge Twin Support Vector Regression
#'
#' \code{hinge_tsvr} is an R implementation of Hinge-TSVR
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param C1,C2 plenty term.
#' @param epsilon1,epsilon2 parameter for epsilon tabe.
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
#' @param solver \code{"dual"} is available.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param randx parameter for reduce SVM, default \code{randx = 0.1}.
#' @param ... unused parameters.
#' @return return \code{TSVRClassifier} object.
#' @export
hinge_tsvr <- function(X, y, C1 = 1, C2 = C1,
                       epsilon1 = 0.1, epsilon2 = epsilon1,
                       kernel = c("linear", "rbf", "poly"),
                       gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                       eps = 1e-7, max.steps = 4000,
                       solver = c("dual"), fit_intercept = TRUE,
                       reduce_set = NULL, ...) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  kernel <- match.arg(kernel)
  solver <- match.arg(solver)
  kso <- kernel_select_option_(X, kernel, reduce_set, gamma, degree, coef0)
  KernelX <- kso$KernelX
  if (fit_intercept == TRUE) {
    KernelX <- cbind(KernelX, 1)
  }
  solver.res <- hinge_tsvr_dual_solver(KernelX, y, C1, C2,
                                         epsilon1, epsilon2,
                                         eps, max.steps)
  TSVRegressor <- list("X" = X, "y" = y,
                       "C1" = C1, "C2" = C2,
                       "epsilon1" = epsilon1, "epsilon2" = epsilon2,
                       "kernel" = kernel,
                       "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                       "solver" = solver, "coef1" = solver.res$coef1,
                       "coef2" = solver.res$coef2,
                       "fit_intercept" = fit_intercept,
                       "Kw" = Kw, solver.res = solver.res)
  class(TSVRegressor) <- "TSVRegressor"
  return(TSVRegressor)
}


#' Predict Method for Twin Support Vector Regression
#'
#' @author Zhang Jiaqi
#' @param object a fitted object of class inheriting from \code{SVMClassifier}.
#' @param X new data for predicting.
#' @param ... unused parameter.
#' @importFrom stats predict
#' @export
predict.TSVRegressor <- function(object, X, ...) {
  X <- as.matrix(X)
  if (object$kernel == "linear") {
    KernelX <- X
  } else {
    KernelX <- kernel_function(X, object$X,
                               kernel.type = object$kernel,
                               gamma = object$gamma,
                               degree = object$degree,
                               coef0 = object$coef0)
  }
  if (object$fit_intercept == TRUE) {
    KernelXe <- cbind(KernelX, 1)
  }
  f <- KernelXe %*% (object$coef1 + object$coef2) / 2
  return(f)
}
