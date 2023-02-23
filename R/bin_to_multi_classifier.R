OVR_Classifier <- function(X, y, model, ...) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  model <- model
  class_num <- length(unique(y))
  Classifiers_OVR <- list()
  if (class_num == 2) {
    class_num <- class_num - 1
  }
  for (i in 1:class_num) {
    Classifiers_OVR[[i]] <- do.call("model", list(X, y, ...))
  }
  class(Classifiers_OVR) <- "Classifiers_OVR"
  return(Classifiers_OVR)
}

