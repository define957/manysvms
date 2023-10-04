bls_svr_primal_solver <- function(KernelX, y, C, b, max.steps,  batch_size,
                                  optimizer, ...) {
  gBLeastSquares <- function(KernelX, y, w, pars,...) {
    # gradient of Bounded Least Squares loss function
    C <- pars$C
    b <- pars$b
    xn <- pars$xn
    xmn <- nrow(KernelX)
    xmp <- ncol(KernelX)
    resi <- y - KernelX%*%w
    delta1 <- 1 + b*(resi^2)
    gterm1 <- resi / (delta1^2)
    g <- w - 2*(C*xn/xmn)*b* t(apply(KernelX, 2, "*", gterm1)) %*% matrix(1, xmn)
    return(g)
  }
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  w0 <- matrix(0, xp, 1)
  pars <- list("C" = C, "xn" = xn, "b" = b)
  wt <- optimizer(KernelX, y, w0, batch_size, max.steps, gBLeastSquares, pars, ...)
  BasePrimalBLSSVMRegressor <- list(coef = as.matrix(wt))
  class(BasePrimalBLSSVMRegressor) <- "BasePrimalBLSSVMRegressor"
  return(BasePrimalBLSSVMRegressor)
}

#' Bounded Least Squares Loss Support Vector Regression
#'
#' \code{bls_svr} is an R implementation of Bounded Least Squares SVR.
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param C plenty term.
#' @param b paramter of epsilon insensitive zone.
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
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param optimizer default primal optimizer pegasos.
#' @param randx parameter for reduce SVM, default \code{randx = 0.1}.
#' @param ... unused parameters.
#' @return return \code{SVMRegressor} object.
#' @export
bls_svr <- function(X, y, C = 1, b = 0.2, kernel = c("linear", "rbf", "poly"),
                    gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                    eps = 1e-5, max.steps = 80, batch_size = nrow(X) / 10,
                    solver = c("primal"),
                    fit_intercept = TRUE, optimizer = nesterov, randx = 0.1, ...) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  kernel <- match.arg(kernel)
  solver <- match.arg(solver)
  if (fit_intercept == TRUE) {
    X <- cbind(X, 1)
  }
  kso <- kernel_select_option(X, kernel, solver, randx,
                              gamma, degree, coef0)

  KernelX <- kso$KernelX
  X <- kso$X
  if (solver == "primal") {
    solver.res <- bls_svr_primal_solver(KernelX, y, C, b, max.steps, batch_size,
                                        optimizer, ...)
  }
  SVMRegressor <- list("X" = X, "y" = y,
                       "C" = C, "kernel" = kernel,
                       "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                       "solver" = solver, "coef" = solver.res$coef,
                       "fit_intercept" = fit_intercept)
  class(SVMRegressor) <- "SVMRegressor"
  return(SVMRegressor)
}
