#  Primal Estimated sub-Gradient SOlver for SVM (Pegasos)
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param w initial point.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param fx sub-gradient of objective function.
#' @return return optimal solution.
#' @export
#' @examples
pegasos <- function(X, y, w, m, max.steps, fx, seed, ...) {
  if (is.null(seed)==FALSE) {
    set.seed(seed)
  }
  C <- list(...)$C
  v <- w
  nx = nrow(X)
  px = ncol(X)
  for (t in 1:max.steps) {
    At <- sample(nx, m)
    xm <- X[At, ]
    dim(xm) <- c(m, px)
    ym <- y[At]
    # update parameter
    dF <- fx(xm, ym, v, ...)
    v <- v - (C/t)*dF
  }
  return(v)
}
