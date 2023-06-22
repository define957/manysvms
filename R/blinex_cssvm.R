blinex_cssvm_primal_solver <- function(KernelX, y, C = 1,
                                       a, b, eps = 1e-5,
                                       max.steps = 80, batch_size = nrow(KernelX) / 10,
                                       optimizer = pegasos, ...) {
  sgBlinex <- function(KernelX, y, v, ...) {
    add_param <- list(...)
    a <- add_param$a
    b <- add_param$b
    C <- add_param$C
    n <- add_param$n
    m <- nrow(KernelX)
    xp <- ncol(KernelX)
    sg <- matrix(0, nrow = xp, ncol = 1)
    xi <- 1 - y * (KernelX %*% v)
    na_idx <- which(is.na(xi)== TRUE)
    print(v)
    xi <- ifelse(xi > 0, xi, 0)
    sgterm <- a*y*xi
    exp_sgterm <- exp(sgterm)
    idx <- which(exp_sgterm!=Inf)
    exp_sgterm <- exp_sgterm[idx]
    sgweight1 <- (1 - exp_sgterm*sign(xi[idx]))
    sgweight2 <- (1 + b*(exp_sgterm - sgterm[idx] - 1))^2
    sg <- v/n + (a*b*C/length(idx))*t((t(sgweight1/sgweight2)%*%KernelX[idx, ]))
    return(sg)
  }
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  w0 <- matrix(0, nrow = xp, ncol = 1)
  wt <- optimizer(KernelX, y, w0, batch_size, max.steps = 80, sgBlinex, C = C,
                  a = a, b = b, n = xn, ...)
  BasePrimalBlinexCSSVMClassifier <- list(coef = as.matrix(wt[1:xp]))
  class(BasePrimalBlinexCSSVMClassifier) <- "BasePrimalBlinexCSSVMClassifier"
  return(BasePrimalBlinexCSSVMClassifier)
}

#' Cost-Sensitive Support Vector Machine with Blinex loss
#'
#' \code{blinex_cssvm} is an R implementation of CSKB
#'
#' @author Li Feihong
#' @param X,y dataset and label.
#' @param C plenty term.
#' @param a,b,c parameters for blinex loss.
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
#' @param solver \code{"primal"} are available.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param optimizer default primal optimizer pegasos.
#' @param randx parameter for reduce SVM, default \code{randx = 0.1}.
#' @param ... unused parameters.
#' @return return \code{HingeSVMClassifier} object.
#' @export
blinex_cssvm <- function(X, y, C = 1, kernel = c("linear", "rbf", "poly"),
                         gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                         a = 1, b = 1, c = 1,
                         eps = 1e-2, max.steps = 80, batch_size = nrow(X) / 10,
                         solver = c("primal"),
                         fit_intercept = TRUE, optimizer = pegasos, randx = 0.1, ...) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  class_set <- sort(unique(y))
  idx <- which(y == class_set[1])
  y[idx] <- 1
  y[-idx] <- -1
  y <- as.matrix(as.numeric(y))
  if (length(class_set) > 2) {
    stop("The number of class should less 2!")
  }
  kernel <- match.arg(kernel)
  if (fit_intercept == TRUE) {
    X <- cbind(X, 1)
  }
  kernel <- match.arg(kernel)
  solver <- match.arg(solver)
  kso <- kernel_select_option(X, kernel, solver, randx,
                              gamma, degree, coef0)
  KernelX <- kso$KernelX
  X <- kso$X
  if (solver == "primal") {
    solver.res <- blinex_cssvm_primal_solver(KernelX, y, C, a, b,
                                             eps, max.steps, batch_size,
                                             optimizer, ...)
  }
  SVMClassifier <- list("X" = X, "y" = y, "class_set" = class_set,
                        "C" = C, "kernel" = kernel,
                        "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                        "solver" = solver, "coef" = solver.res$coef,
                        "fit_intercept" = fit_intercept)
  class(SVMClassifier) <- "SVMClassifier"
  return(SVMClassifier)
}