#' Split Training and Testting data
#'
#' Return training and testting dataset.
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param test_size the size of test dataset in [0, 1].
#' @param shuffle shuffle the dataset.
#' @param seed random seed option.
#' @export

train_test_split <- function(X, y, test_size = 0.3,
                             shuffle = TRUE, seed = NULL){
  X <- as.matrix(X)
  y <- as.matrix(y)
  m <- nrow(X)
  idx <- 1:round(m*(1 - test_size))
  if (shuffle == TRUE) {
    if (is.null(seed) == FALSE) {
      set.seed(seed)
    }
    new_idx <- sample(idx)
  }else{
    new_idx <- idx
  }
  train_X <- X[new_idx, ]
  train_y <- y[new_idx]
  test_X <- X[-new_idx, ]
  test_y <- y[-new_idx]
  train_test_data <- list("train_X" = train_X, "train_y" = train_y,
                          "test_X" = test_X, "test_y" = test_y)
  return(train_test_data)
}


#' Generate noisy label
#'
#' @author Zhang Jiaqi
#' @param y clean label.
#' @param p ratio of noisy.
#' @param seed random seed.
#' @export

noisy_label_generator <- function(y, p, seed = NULL){
  if (is.null(seed) == FALSE) {
    set.seed(seed)
  }
  y <- as.matrix(y)
  class_set <- unique(y)
  class_num <- length(class_set)
  class_idx <- list()
  for (i in 1:class_num) {
    idx <- which(y == class_set[i])
    class_idx[[i]] <- idx
  }
  for (i in 1:class_num) {
    m <- length(class_idx[[i]])
    n <- round(m*p, 0)
    idx_temp <- sample(m, n)
    noisy_y <- sample(class_set[class_set != class_set[i]], n, replace = T)
    y[class_idx[[i]][idx_temp]] <- noisy_y
  }
  return(y)
}
