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

#' Binary F1 Score
#' return F1 Score for two vector
#' @author Zhang Jiaqi
#' @param y,y_hat real values and fitted values.
#' @param positive positive class, default "Class-1".
#' @export
binaryf1score <- function(y, y_hat, positive = "Class-1") {
  y_hat <- as.character(y_hat)
  n1 <- length(y)
  n2 <- length(y_hat)
  if (n1 != n2) {
    stop("y and y_hat should have same rows.")
  }
  class_set <- sort(unique(y))
  if (length(class_set) > 2) {
    stop("The number of class should be 2.")
  }
  if (positive %in% class_set == FALSE) {
    return(0)
  } else if (any(unique(y_hat) %in% y) == FALSE) {
    return(0)
  }
  conf_DF <- as.data.frame(confusion_mat(y, y_hat))
  conf_DF$Freq <- as.integer(conf_DF$Freq)
  TP <- as.integer(subset(conf_DF, y == positive &
                            y_hat == positive)["Freq"])
  FN <- as.integer(sum(subset(conf_DF, y == positive &
                                y_hat != positive)["Freq"]))
  FP <- as.integer(sum(subset(conf_DF, y != positive &
                                y_hat == positive)["Freq"]))
  if (is.na(TP)) {
    TP <- 0
  }
  if (is.na(FP)) {
    FP <- 0
  }
  if (is.na(FN)) {
    FN <- 0
  }
  F1 <- 2 * TP/(2 * TP + FP + FN)
  return(F1)
}

confusion_mat <- function(y, y_hat) {
  confusion_mat <- table(y, y_hat)
  return(confusion_mat)
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
