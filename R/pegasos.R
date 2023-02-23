pegasos <- function(X, y, w, m, max.steps, fx, ...) {
  C <- list(...)$C
  v <- w
  nx = nrow(X)
  px = ncol(X)
  for (t in 1:max.steps) {
    At <- sample(nx, m)
    xm <- X[At, ]
    ym <- y[At]
    # update parameter
    dF <- fx(xm, ym, v, ...)
    v <- v - (C/t)*dF
  }
  return(v)
}
