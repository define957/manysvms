sigmoid_svm_dual_solver <- function(KernelX, y, C = 1, update_deltak,
                                    epsilon = 0, lambda = 1,
                                    eps = 1e-5, eps.cccp = 1e-2, max.steps = 80, cccp.steps = 10) {
  n <- nrow(KernelX)
  H <- calculate_svm_H(KernelX, y)
  e <- matrix(1, nrow = n, ncol = 1)
  u0 <- matrix(0, nrow = n, ncol = 1)
  for (i in 1:cccp.steps) {
    f <- 1 - H %*% u0 - epsilon
    delta_k <- update_deltak(f, u0, epsilon, lambda)
    lb <- -C*delta_k
    ub <- -C*delta_k + lambda*C
    u <- clip_dcd_optimizer(H, (1 - epsilon)*e, lb, ub, eps, max.steps, u0)$x
    if (norm(u - u0, type = "2") < eps.cccp) {
      break
    } else {
      u0 <- u
    }
  }
  coef <- y*u
  BaseDualSigmoidSVMClassifier <- list(coef = as.matrix(coef))
  class(BaseDualSigmoidSVMClassifier) <- "BaseDualSigmoidSVMClassifier"
  return(BaseDualSigmoidSVMClassifier)
}


sigmoid_svm_primal_solver <- function(KernelX, y, C = 1, update_deltak,
                                      epsilon = 0, lambda = 1,
                                      eps = 1e-5, eps.cccp = 1e-2,
                                      max.steps = 80, cccp.steps = 10, batch_size = nrow(KernelX) / 10,
                                      optimizer = pegasos, ...) {
  sgSigmoid <- function(KernelX, y, v, epsilon, lambda, deltak, At, C, ...) { # sub-gradient of Sigmoid loss function
    xn <- nrow(KernelX)
    xp <- ncol(KernelX)
    sg <- matrix(0, nrow = xp, ncol = 1)
    u <- 1 - y*(KernelX %*% v) - epsilon
    u[u < 0] <- 0
    u[u >= 0] <- 1
    sg <- v - lambda*(C/xn) * t(KernelX) %*% (u*y) +
      C/xn*t(KernelX) %*% (y*deltak[At, ])
    return(sg)
  }
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  wt <- matrix(0, nrow = xp, ncol = 1)
  for (i in 1:cccp.steps) {
    f <- 1 - y*(KernelX %*% wt)
    deltak <- update_deltak(f, wt, epsilon, lambda)
    wt <- optimizer(KernelX, y, wt, batch_size, max.steps,
                    sgSigmoid, eps, C = C,
                    epsilon = epsilon, lambda = lambda, deltak = deltak, ...)
  }
  BasePrimalSigmoidSVMClassifier <- list(coef = as.matrix(wt[1:xp]))
  class(BasePrimalSigmoidSVMClassifier) <- "BasePrimalSigmoidSVMClassifier"
  return(BasePrimalSigmoidSVMClassifier)
}


#' Sigmoid Support Vector Machine
#'
#' \code{sigmoid_svm} is an R implementation of Sigmoid-SVM
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
#' @param epsilon parameter for sigmoid loss (epsilon-insensitive zone).
#' @param lambda parameter for sigmoid loss (loss increase speed).
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param eps the precision of the optimization algorithm.
#' @param eps.cccp the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param cccp.steps the number of iterations of CCCP.
#' @param batch_size mini-batch size for primal solver.
#' @param solver \code{"dual"} and \code{"primal"} are available.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param optimizer default primal optimizer pegasos.
#' @param randx parameter for reduce SVM, default \code{randx = 0.1}.
#' @param ... unused parameters.
#' @return return \code{SVMClassifier} object.
#' @export
sigmoid_svm <- function(X, y, C = 1, kernel = c("linear", "rbf", "poly"),
                        gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                        epsilon = 0, lambda = 1,
                        eps = 1e-5, eps.cccp = 1e-2, max.steps = 80, cccp.steps = 10,
                        batch_size = nrow(X) / 10,
                        solver = c("dual", "primal"),
                        fit_intercept = TRUE, optimizer = pegasos, randx = 0.1, ...) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  class_set <- unique(y)
  idx <- which(y == class_set[1])
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
  kso <- kernel_select_option(X, kernel, solver, randx,
                              gamma, degree, coef0)
  KernelX <- kso$KernelX
  X <- kso$X
  update_deltak <- function(f, u, epsilon, lambda) {
    idx <- which(f >= 0)
    ef <- exp(-lambda*f)
    delta_k <- lambda*(1 - 2*ef/((1 + ef)^2))
    delta_k[-idx] <- 0
    return(delta_k)
  }
  if (solver == "primal") {
    solver.res <- sigmoid_svm_primal_solver(KernelX, y, C, update_deltak,
                                            epsilon, lambda,
                                            eps, eps.cccp, max.steps, cccp.steps, batch_size,
                                            optimizer, ...)
  } else if (solver == "dual") {
    solver.res <- sigmoid_svm_dual_solver(KernelX, y, C, update_deltak,
                                          epsilon, lambda, eps, eps.cccp,
                                          max.steps, cccp.steps)
  }
  SVMClassifier <- list("X" = X, "y" = y, "class_set" = class_set,
                        "C" = C, "kernel" = kernel,
                        "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                        "solver" = solver, "coef" = solver.res$coef,
                        "fit_intercept" = fit_intercept)
  class(SVMClassifier) <- "SVMClassifier"
  return(SVMClassifier)
}









