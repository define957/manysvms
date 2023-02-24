OVR_Classifier <- function(X, y, model, ...) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  model <- model
  class_set <- unique(y)
  class_num <- length(class_set)
  classifier_num <- class_num
  Classifiers_OVR <- list()
  classifier_list <- list()
  y_temp <- rep(0, nrow(X))
  for (i in 1:classifier_num) {
    idx_pos <- which(y == class_set[i])

    y_temp[idx_pos] <- 1
    y_temp[-idx_pos] <- -1
    classifier_list[[i]] <- do.call("model", list(X, y_temp, ...))
  }
  Classifiers_OVR$classifier_list <- classifier_list
  Classifiers_OVR$class_set <- class_set
  Classifiers_OVR$class_num <- class_num
  class(Classifiers_OVR) <- "Classifiers_OVR"
  return(Classifiers_OVR)
}

predic.OVR_Classifier <- function(object, X, y, predict_func = predict, ...) {
  n <-  nrow(X)
  classifiers <- object$classifier_list
  vote_mat <- matrix(0, nrow = n, ncol = object$class_num)
  for (i in 1:length(classifiers)) {
    pred <- predict.hinge_svm(classifiers[[i]], X, y, ...)
    idx <- which(pred+1 == y)
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
