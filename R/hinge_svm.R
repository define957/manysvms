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
  return(BaseDualHingeSVMClassifier)
}


hinge_svm_primal_solver <- function (X, y, C = 1, eps = 1e-5,
                                     max.steps = 80, batch_size = nrow(X) / 10, ...) {
  sgHinge <- function(X, y, v, ...) { # sub-gradient of hinge loss function
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
  xn <- nrow(X)
  xp <- ncol(X)
  xx <- cbind(X, 1)
  w0 <- matrix(0, nrow = xp+1, ncol = 1)
  wt <- pegasos(xx, y, w0, batch_size, max.steps, sgHinge, C = C)
  wnorm <- norm(wt[1:xp], type = "2")
  if (fit_intercept == TRUE) {
    BasePrimalHingeSVMClassifier <- list(coef = as.matrix(wt[1:xp]),
                                         intercept = wt[xp+1])
  } else {
    BasePrimalHingeSVMClassifier <- list(coef = as.matrix(wt[1:xp]),
                                         intercept = 0)
  }
  class(BasePrimalHingeSVMClassifier) <- "BasePrimalHingeSVMClassifier"
  return(BasePrimalHingeSVMClassifier)
}

hinge_svm <- function (X, y, C = 1, kernel = c("linear", "rbf", "poly"),
                       gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                       eps = 1e-5, max.steps = 80, batch_size = nrow(X) / 10,
                       solver = c("dual", "primal"), rcpp = TRUE) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  kernel <- match.arg(kernel)
  solver <- match.arg(solver)
  if (kernel == "linear" & solver == "primal") {
    KernelX <- X
  } else {
    KernelX <- kernel_function(X, X,
                               kernel.type = kernel,
                               gamma = gamma, degree = degree, coef0 = coef0,
                               rcpp = rcpp)
  }
  if (solver == "primal") {
    solver.res <- hinge_svm_primal_solver(KernelX, y, C, eps,
                                          max.steps, batch_size)
  } else if(solver == "dual") {
    solver.res <- hinge_svm_dual_solver(KernelX, y, C, eps,
                                        max.steps, rcpp)
  }
  HingeSVMClassifier <- list("C" = C, "kernel" = kernel,
                             "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                             "solver" = solver,
                             "coef" = solver.res$coef)
  class(HingeSVMClassifier) <- "HingeSVMClassifier"
  return(HingeSVMClassifier)
}









