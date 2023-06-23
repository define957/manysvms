#'  Primal Estimated sub-Gradient SOlver for SVM (Pegasos)
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param w initial point.
#' @param m mini-batch size for pegasos solver.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param fx sub-gradient of objective function.
#' @param eps the precision of the optimization algorithm.
#' @param C penalty term of SVMs.
#' @param ... additional settings for the sub-gradient.
#' @return return optimal solution.
#' @references ${1:Pegasos: Primal Estimated sub-GrAdient SOlver for SVM}
#' @export
pegasos <- function(X, y, w, m, max.steps, fx, eps = 1e-5, C = 1, ...) {
  v <- w
  nx = nrow(X)
  px = ncol(X)
  for (t in 1:max.steps) {
    At <- sample(nx, m)
    xm <- X[At, ]
    dim(xm) <- c(m, px)
    ym <- as.matrix(y[At])
    # update parameter
    dF <- fx(xm, ym, v, At = At, C = C, ...)
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

#' Nesterov's Accelerated Gradient Descent
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param w initial point.
#' @param m mini-batch size for pegasos solver.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param fx sub-gradient of objective function.
#' @param eps the precision of the optimization algorithm.
#' @param v initial velocity.
#' @param eta initial learning rate.
#' @param gam momentum parameter.
#' @param k decay rate of learning rate.
#' @param ... additional settings for the sub-gradient.
#' @return return optimal solution.
#' @export
nesterov <- function(X, y, w, m, max.steps, fx, eps = 1e-5,
                     v = matrix(0, nrow(w)), eta = 1, gam = 0.5, k = 0.5, ...) {
  nx <- nrow(X)
  px <- ncol(X)
  for (t in 1:max.steps) {
    At <- sample(nx, m)
    xm <- X[At, ]
    dim(xm) <- c(m, px)
    ym <- as.matrix(y[At])
    v <- gam*v - eta*fx(xm, ym, w, At = At, ...)
    wk <- w + gam*v
    eta <- eta*exp(-k)
    if (norm(wk - w, type = "2") < eps) {
      break
    } else {
      w <- wk
    }
  }
  return(w)
}

#' Rmsprop Optimizer
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param w initial point.
#' @param m mini-batch size for pegasos solver.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param fx sub-gradient of objective function.
#' @param eps the precision of the optimization algorithm.
#' @param epsilon initial stepsize.
#' @param rho momentum parameter.
#' @param delta avoid division by 0.
#' @param ... additional settings for the sub-gradient.
#' @return return optimal solution.
#' @export
rmsprop <- function(X, y, w, m, max.steps, fx, eps = 1e-5,
                    epsilon = 0.001, rho = 0.9, delta = 1e-5,...) {
  v <- w
  xn <- nrow(X)
  xp <- ncol(X)
  r <- matrix(0.1,nrow = ncol(X),ncol = 1)
  g <- matrix(0.1,nrow = ncol(X),ncol = 1)
  for (t in 1:max.steps) {
    At <- sample(xn, m)
    xm <- X[At, ]
    dim(xm) <- c(m, xp)
    ym <- as.matrix(y[At])
    # update parameter
    dF <- fx(xm, ym, v, ...)
    rk <- rho*r + (1 - rho)*g*g
    vk <- v - (epsilon/sqrt(delta + rk)) * dF
    g <- dF
    r <- rk
    if (norm(vk - v, type = "2") < eps) {
      break
    } else {
      v <- vk
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
