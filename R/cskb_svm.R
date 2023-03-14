#' Hinge Support Vector Machine
#'
#' \code{hinge_svm} is an R implementation of Hinge-SVM
#'
#' @author Li Feihong
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
cskb_svm <- function(X, y, a = 1, b = 1, C = 1, kernel = c("linear", "rbf", "poly"),
                       gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                       eps = 1e-2, max.steps = 80, batch_size = nrow(X) / 10,
                       fit_intercept = TRUE, randx = 0.1, ...) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  class_set <- unique(y)
  idx <- which(y == class_set[1])
  y[idx] <- 1
  y[-idx] <- -1
  if (length(class_set) > 2) {
    stop("The number of class should less 2!")
  }
  kernel <- match.arg(kernel)
  if (fit_intercept == TRUE) {
    X <- cbind(X, 1)
  }
  # define a function to compute the stochastic gradient for linear case
  sgFun <- function(w, sqnInd){#
    xiTerm <- 1 - y[sqnInd] * c( X[sqnInd,] %*% w)
    xi <-  ifelse(xiTerm > 0, xiTerm, 0)
    sgTerm <- a*y[sqnInd]*xi
    sgWeight1 <- (1 - exp(sgTerm))*sign(xi)
    sgWeight2 <- (1 + b*(exp(sgTerm) - sgTerm - 1))^2
    sg <- w/n + (a * b * C / m)*t((sgWeight1/sgWeight2) %*% X[sqnInd,])
    return(sg) # sg is a column vector
  }
  n <- nrow(X)
  p <- ncol(X)
  if (kernel == "linear") {
    KernelX <- X
  }
  else if (kernel != "linear") {
    if (randx > 0) {
      randX = X[sample(nrow(X), floor(randx*nrow(X))),]
    }
    KernelX <- kernel_function(X, randX,
                               kernel.type = kernel,
                               gamma = gamma, degree = degree, coef0 = coef0)
    X <- randX
  }

  if (kernel == "linear") {
    w0 <- matrix(0, nrow = p, ncol = 1) # initial normal vector
    v0 <- matrix(0, nrow = p, ncol = 1) # initial velocity
    eta0 <- 1 # initial learning rate
    gam <- 0.5 # the momentum coefficient
    k <- 0.5 # the learning rate decay factor
    m <- 10 # the number of the randomly selected samples

    # the first update

    v <- gam*v0 - eta0*sgFun(w0 + gam*v0, sample(n, m))
    w <- w0 + v
    eta <- eta0*exp(-k)

    t <- 2
    while (t <= max.steps && sqrt(sum(v^2)) >=  eps) {
      eta0 <- eta
      v0 <- v
      w0 <- w

      # update velocity
      v <- gam*v0 - eta0*sgFun(w0 + gam*v0, sample(n, m))
      v[is.na(v)] <- 0
      # update model parameter
      w <- w0 + v
      # update the learning rate
      eta <- eta0*exp(-k*t)
      t <- t + 1
    }
    coef <- w
  }
  if (kernel == "gaussian") {
    #给定一系列的初始值
    a0 <- c(rep(0,n))
    it <- 1
    buchang <- 1
    while (it < max.steps && abs(buchang >= eps)) {
      index <- sample(1:n,size = 1,replace = F)
      if (it == 1) {
        kesai <- 1
      }
      if ( it > 1) {
        kesai = 1 - (y[index] * a * b * C/(it - 1))*a0 %*% KernelX[index,]
      }
      if (!is.na(kesai)) {
        if (kesai > 0) {
          buchang = (exp( a * y[index]*kesai) - 1)/((1 + b*(exp(a * y[index] * kesai) - a * y[index]*kesai - 1))^(2))
          a0[index] = a0[index] + buchang
        }else{
          buchang = 0
        }
      }
      it = it + 1
    }
    coef = a0
  }
  SVMClassifier <- list("X" = X, "y" = y, "class_set" = class_set,
                        "C" = C, "kernel" = kernel,
                        "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                        "coef" = coef, "solver" = "primal",
                        "fit_intercept" = fit_intercept)
  class(SVMClassifier) <- "SVMClassifier"
  return(SVMClassifier)
}
