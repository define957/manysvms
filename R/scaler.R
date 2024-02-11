#' Standard Scaler
#'
#' @author Zhang Jiaqi
#' @param X data matrix.
#' @export
standar_scaler <- function(X) {
  X <- scaler_check(X)
  X_mean <- apply(X, 2, mean)
  X_sd  <- apply(X, 2, sd)
  StandardScaler <- list("X_mean" = X_mean, "X_sd" = X_sd)
  class(StandardScaler) <- "StandardScaler"
  return(StandardScaler)
}

#' Scaler Transform Function
#'
#' @author Zhang Jiaqi
#' @param scaler data matrix.
#' @param X data matrix.
#' @export
trans <- function(scaler, X) {
  UseMethod("trans")
}

#' @exportS3Method manysvms::trans
trans.StandardScaler <- function(scaler, X) {
  X_mean <- scaler$X_mean
  X_sd   <- scaler$X_sd
  X_scaled <- t((t(X) - X_mean) / X_sd)
  return(X_scaled)
}

scaler_check <- function(X) {
  if (is.vector(X) == TRUE) {
    X <- as.matrix(X)
  }
  return(X)
}
