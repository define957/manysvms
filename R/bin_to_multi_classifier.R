#' Convert Binary Classifier to multiclass classifier
#' by Using Ones versus Rest Strategy
#'
#' \code{pegasos} is an R implementation of Pegasos algorithm
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param bin_model binary classifier
#' @param ... parameter of model.
#' @return return \code{Classifiers_OVR} object.
#' @export
OVR_Classifier <- function(X, y, bin_model, ...) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  model <- as.character(match.call(bin_model)[-1L])[3]
  class_set <- unique(y)
  class_num <- length(class_set)
  classifier_num <- class_num
  Classifiers_OVR <- list()
  classifier_list <- list()
  y_temp <- rep(0, nrow(X))
  for (i in 1:classifier_num) {
    idx_pos <- which(y == class_set[i])
    y_temp[idx_pos] <- -1
    y_temp[-idx_pos] <- 1
    classifier_list[[i]] <- do.call("bin_model", list(X, y_temp, ...))
  }
  Classifiers_OVR$classifier_list <- classifier_list
  Classifiers_OVR$class_set <- class_set
  Classifiers_OVR$class_num <- class_num
  class(Classifiers_OVR) <- "Classifiers_OVR"
  return(Classifiers_OVR)
}

#' Predict Method of OVR Classifier
#'
#' @author Zhang Jiaqi.
#' @param object a fitted object of class inheriting from \code{Classifiers_OVR}.
#' @param X new data for predicting.
#' @param predict_func predict function of your model.
#' @param ... parameters of your model.
#' @return return predict results.
#' @export
predict.Classifiers_OVR <- function(object, X, predict_func = predict, ...) {
  n <-  nrow(X)
  classifiers <- object$classifier_list
  vote_mat <- matrix(0, nrow = n, ncol = object$class_num)
  for (i in 1:length(classifiers)) {
    pred <- predict_func(classifiers[[i]], X, ...)
    idx <- which(pred == -1)
    vote_mat[idx, i] <- vote_mat[idx, i] + 1
  }
  res_vote <- apply(vote_mat, 1, which.max)
  res <- matrix(0, nrow = n)
  for (i in 1:object$class_num) {
    idx <- which(res_vote == i)
    res[idx] <- object$class_set[i]
  }
  return(res)
}
