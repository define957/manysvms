#' Mean Squared Error
#' return mean squared error for two vector
#' @author Zhang Jiaqi
#' @param y,y_hat real values and fitted values.
#' @param dof degree of freedom.
#' @export
mean_squared_error <- function(y, y_hat, dof = 0) {
  y <- as.matrix(y)
  y_hat <- as.matrix(y_hat)
  n1 <- nrow(y)
  n2 <- nrow(y_hat)
  if (n1 != n2) {
    stop("y and y_hat should have same rows.")
  }
  n <- n1
  mse <- sum((y - y_hat)^2) / (n - dof)
  return(mse)
}

#' Compute Accuracy
#'
#' @author Zhang Jiaqi
#' @param y,y_hat real values and fitted values.
#' @export
accuracy <- function(y, y_hat){
  y <- as.matrix(y)
  y_hat <- as.matrix(y_hat)
  n1 <- nrow(y)
  n2 <- nrow(y_hat)
  if (n1 != n2) {
    stop("y and y_hat should have same rows.")
  }
  acc <- sum((y == y_hat)) / n1
  return(acc)
}

#' Mean Absolute Error
#'
#' @author Zhang Jiaqi
#' @param y,y_hat real values and fitted values.
#' @export
mean_absolute_error <- function(y, y_hat) {
  y <- as.matrix(y)
  y_hat <- as.matrix(y_hat)
  n1 <- nrow(y)
  n2 <- nrow(y_hat)
  if (n1 != n2) {
    stop("y and y_hat should have same rows.")
  }
  n <- n1
  mae <- sum(abs((y - y_hat))) / (n)
  return(mae)
}

#' Root Mean Squared Error
#' return Root Mean Squared Error for two vector
#' @author Zhang Jiaqi
#' @param y,y_hat real values and fitted values.
#' @param dof degree of freedom.
#' @export
root_mean_squared_error <- function(y, y_hat, dof = 0) {
  mse <- mean_squared_error(y, y_hat, dof)
  rmse <- sqrt(mse)
  return(rmse)
}

#' F1 Score
#' return F1 Score for two vector
#' @author Zhang Jiaqi
#' @param y,y_hat real values and fitted values.
#' @export
f1score <- function(y, y_hat) {
  y <- as.matrix(y)
  y_hat <- as.matrix(y_hat)
  n1 <- nrow(y)
  n2 <- nrow(y_hat)
  if (n1 != n2) {
    stop("y and y_hat should have same rows.")
  }
  y_uni <- unique(y)
  idx <- which(y == y_uni[1])
  y[idx] <- 1
  y[-idx] <- -1
  idx <- which(y_hat == y_uni[1])
  y_hat[idx] <- 1
  y_hat[-idx] <- -1
  y <- factor(y, levels = c("1", "-1"))
  y_hat <- factor(y_hat, levels = c("1", "-1"))
  conf_mat <- table(y_hat, y)
  tp <- conf_mat[1, 1]
  tn <- conf_mat[2, 2]
  fp <- conf_mat[2, 1]
  fn <- conf_mat[1, 2]
  f1s <- 2*tp / (2*tp+fp+fn)
  return(f1s)
}

#' Recall
#' return recall for two vector
#' @author Zhang Jiaqi
#' @param y,y_hat real values and fitted values.
#' @export
recall <- function(y, y_hat) {
  y <- as.matrix(y)
  y_hat <- as.matrix(y_hat)
  n1 <- nrow(y)
  n2 <- nrow(y_hat)
  if (n1 != n2) {
    stop("y and y_hat should have same rows.")
  }
  y_uni <- unique(y)
  idx <- which(y == y_uni[1])
  y[idx] <- 1
  y[-idx] <- -1
  idx <- which(y_hat == y_uni[1])
  y_hat[idx] <- 1
  y_hat[-idx] <- -1
  y <- factor(y, levels = c("1", "-1"))
  y_hat <- factor(y_hat, levels = c("1", "-1"))
  conf_mat <- table(y_hat, y)
  tp <- conf_mat[1, 1]
  fn <- conf_mat[1, 2]
  rec <- tp / (tp+fn)
  return(rec)
}

#' Precision
#' return precision for two vector
#' @author Zhang Jiaqi
#' @param y,y_hat real values and fitted values.
#' @export
precision <- function(y, y_hat) {
  y <- as.matrix(y)
  y_hat <- as.matrix(y_hat)
  n1 <- nrow(y)
  n2 <- nrow(y_hat)
  if (n1 != n2) {
    stop("y and y_hat should have same rows.")
  }
  y_uni <- unique(y)
  idx <- which(y == y_uni[1])
  y[idx] <- 1
  y[-idx] <- -1
  idx <- which(y_hat == y_uni[1])
  y_hat[idx] <- 1
  y_hat[-idx] <- -1
  y <- factor(y, levels = c("1", "-1"))
  y_hat <- factor(y_hat, levels = c("1", "-1"))
  conf_mat <- table(y_hat, y)
  tp <- conf_mat[1, 1]
  fp <- conf_mat[2, 1]
  pre <- tp / (tp+fp)
  return(pre)
}
