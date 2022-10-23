#' Mean Squared Error
#' return mean squared error for two vector
#' @author Zhang Jiaqi
#' @param y,y_hat real values and fitted values.
#' @export
mean_squared_error <- function(y, y_hat){
  y <- as.matrix(y)
  y_hat <- as.matrix(y_hat)
  n1 <- nrow(y)
  n2 <- nrow(y_hat)
  if (n1 != n2) {
    stop("y and y_hat should have same rows.")
  }
  n <- n1
  mse <- sum((y - y_hat)^2) / n
  return(mse)
}

