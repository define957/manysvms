roboss_svm_primal_solver <- function(KernelX, X, y, C, a, lambda,
                                     max.steps, batch_size,
                                     optimizer, kernel, reduce_flag,
                                     reduce_set,
                                     ...) {
  groboss <- function(batch_KernelX, y, w, pars, ...) {
    C <- pars$C
    a <- pars$a
    lambda <- pars$lambda
    xn <- pars$xn
    KernelX <- pars$KernelX
    xmn <- nrow(batch_KernelX)
    xmp <- ncol(batch_KernelX)
    u <- 1 - y * (batch_KernelX %*% w)
    idx <- which(u >= 0)
    g <- matrix(0, xmp, 1)
    if (pars$kernel == "linear" || pars$reduce_flag) {
      g <- w / xn -
        (C * lambda * a^2/xmn) * t(batch_KernelX[idx, , drop = FALSE]) %*% (y[idx]*u[idx]*exp(-a*u[idx]))
    } else if (pars$kernel != "linear") {
      g <- KernelX %*% w * xn -
        (C * lambda * a^2/xmn) * t(batch_KernelX[idx, , drop = FALSE]) %*% (y[idx]*u[idx]*exp(-a*u[idx]))
    }
    return(g)
  }
  xn <- nrow(X)
  if (kernel == "linear") { xp <- ncol(X) } else { xp <- ncol(KernelX)}
  w0 <- matrix(0, nrow = xp, ncol = 1)
  pars <- list("C" = C, "a" = a, "lambda" = lambda, "xn" = xn, "kernel" = kernel,
               "KernelX" = KernelX, "reduce_flag" = reduce_flag)
  if (kernel == "linear") {
    wt <- optimizer(X, y, w0, batch_size, max.steps, groboss, pars, ...)
  } else {
    wt <- optimizer(KernelX, y, w0, batch_size, max.steps, groboss, pars, ...)
  }
  BasePrimalrobossSVMClassifier <- list(coef = as.matrix(wt[1:xp]))
  class(BasePrimalrobossSVMClassifier) <- "BasePrimalrobossSVMClassifier"
  return(BasePrimalrobossSVMClassifier)
}


#' Support Vector Machine With RoBoSS Loss Function
#'
#' \code{roboss_svm} is an R implementation of RoBoSS-SVM
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
#' @param a,lambda parameter for roboss loss.
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param batch_size mini-batch size for primal solver.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param optimizer default primal optimizer pegasos.
#' @param reduce_set reduce set for reduce SVM, default \code{reduce_set = NULL}.
#' @param ... unused parameters.
#' @return return \code{SVMClassifier} object.
#' @export
roboss_svm <- function(X, y, C = 1, kernel = c("linear", "rbf", "poly"),
                       gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                       a = 1, lambda = 1,
                       max.steps = 4000, batch_size = nrow(X) / 10,
                       fit_intercept = TRUE, optimizer = nesterov,
                       reduce_set = NULL,
                       ...) {
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
  solver <- "primal"
  if (fit_intercept == TRUE) {
    X <- cbind(X, 1)
  }
  reduce_flag <- is.null(reduce_set) == FALSE
  if (solver == "dual" && reduce_flag == TRUE) {
    reduce_flag <- FALSE
    reduce_set <- NULL
    cat("The dual solver does not support the reduce set; it has been set to NULL.\n")
  }
  kso <- kernel_select_option_(X, kernel, reduce_set, gamma, degree, coef0)
  KernelX <- kso$KernelX
  if (solver == "primal") {
    solver.res <- roboss_svm_primal_solver(KernelX, X, y, C, a, lambda,
                                           max.steps, batch_size,
                                           optimizer, kernel, reduce_flag,
                                           reduce_set,
                                           ...)
  }
  SVMClassifier <- list("X" = X, "y" = y,
                        "reduce_flag" = reduce_flag,
                        "reduce_set" = reduce_set,
                        "class_set" = class_set,
                        "C" = C, "kernel" = kernel,
                        "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                        "solver" = solver, "coef" = solver.res$coef,
                        "fit_intercept" = fit_intercept)
  class(SVMClassifier) <- "SVMClassifier"
  return(SVMClassifier)
}

