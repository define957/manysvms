#'  Primal Estimated sub-Gradient SOlver for SVM (Pegasos)
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param w initial point.
#' @param m mini-batch size for pegasos solver.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param fx sub-gradient of objective function.
#' @param eps the precision of the optimization algorithm.
#' @param ... additional settings for the sub-gradient.
#' @return return optimal solution.
#' @references ${1:Pegasos: Primal Estimated sub-GrAdient SOlver for SVM}
#' @export
pegasos <- function(X, y, w, m, max.steps, fx, eps = 1e-5, ...) {
  C <- list(...)$C
  v <- w
  nx = nrow(X)
  px = ncol(X)
  for (t in 1:max.steps) {
    At <- sample(nx, m)
    xm <- X[At, ]
    dim(xm) <- c(m, px)
    ym <- as.matrix(y[At])
    # update parameter
    dF <- fx(xm, ym, v, At = At, ...)
    v <- v - (C/t)*dF
    v <- min(1, sqrt(C)/norm(v, type = "2"))*v
    if (norm(v - w, type = "2") < eps) {
      break
    } else {
      w <- v
    }
  }
  return(v)
}

#'  Conjugate Gradient Method for Solving Linear Equation Ax = b
#'
#' @author Zhang Jiaqi.
#' @param A,b matrix of linear equation Ax = b (Since A is a PSD matrix).
#' @param x initial point.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param eps the precision of the optimization algorithm.
#' @param ... additional settings for the sub-gradient.
#' @return return optimal solution.
#' @export
conjugate_gradient_method <- function(A, b, x, max.steps, eps = 1e-5, ...) {
  rk <- b - A%*%x
  pk <- rk
  for (t in 1:max.steps) {
    rk2 <- (t(rk) %*% rk)
    alphak <-  as.numeric(rk2 / (t(pk) %*% A %*% pk))
    x <- x + alphak*pk
    rk_1 <- rk - alphak*(A%*%pk)
    if (norm(rk_1, type = "2") < eps) {
      cat("converge after", t, "steps", "\n")
      break
    } else {
      betak <- as.numeric((t(rk_1) %*% rk_1) / rk2)
      pk <- rk_1 + betak*pk
    }
    rk <- rk_1
  }
  return(x)
}
