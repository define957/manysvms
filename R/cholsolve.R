cholsolve <- function(a, b) {
  Lt <- chol(a, pivot = T)
  pivot <- attr(Lt, "pivot")
  oo <- order(pivot)
  x <- backsolve(Lt, forwardsolve(Lt, b[pivot, ], upper.tri = T, transpose = T))
  dim(x) <- c(ncol(a), ncol(b))
  x <- x[oo, ]
  return(x)
}
 