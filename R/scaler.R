scaler_check <- function(X) {
  if (is.vector(X) == TRUE) {
    X <- as.matrix(X)
  }
  ScalerChecker <- list("X" = X, "num_features" = ncol(X))
  class(ScalerChecker) <- "ScalerChecker"
  return(ScalerChecker)
}

trans_check <- function(scaler, X) {
  if (is.vector(X) == TRUE) {
    X <- as.matrix(X)
  }
  X_features <- ncol(X)
  if (X_features != scaler$num_features) {
    stop("X should have sample features with scaler!")
  }
  TransChecker <- list("X" = X,
                       "num_samples" = nrow(X), "num_features" = ncol(X))
  class(TransChecker) <- "TransChecker"
  return(TransChecker)
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

#' Standard Scaler
#'
#' @author Zhang Jiaqi
#' @param X data matrix.
#' @export
standar_scaler <- function(X) {
  scaler_checker <- scaler_check(X)
  X <- scaler_checker$X
  X_mean <- apply(X, 2, mean)
  X_sd  <- apply(X, 2, sd)
  StandardScaler <- list("X_mean" = X_mean, "X_sd" = X_sd,
                         "num_features" = scaler_checker$num_features)
  class(StandardScaler) <- "StandardScaler"
  return(StandardScaler)
}

#' @exportS3Method manysvms::trans
trans.StandardScaler <- function(scaler, X) {
  TransChecker <- trans_check(scaler, X)
  X            <- TransChecker$X
  X_scaled     <- matrix(0, TransChecker$num_samples, TransChecker$num_features)
  X_mean       <- scaler$X_mean
  X_sd         <- scaler$X_sd
  idx          <- which(X_sd != 0)
  # Calculate centered X
  X_centered_T     <- (t(X) - X_mean)
  # Calculate Standarded X (sd != 0)
  X_scaled[, idx]  <- t(X_centered_T[idx, ] / X_sd[idx])
  X_scaled[, -idx] <- X_centered_T[-idx, ]
  return(X_scaled)
}

