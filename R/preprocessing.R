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
  m <- nrow(m)
  if(shuffle == TRUE){
    if(is.null(seed)==FALSE){
      set.seed(seed)
    }
    new_idx <- sample(1:m)
  }else{
    new_idx <- 1:m
  }
  idx <- round(m*0.3)
  train_X <- X[1:idx, ]
  train_y <- X[1:idx]
  test_X <- X[(idx+1):m, ]
  test_y <- y[(idx+1):m]
  train_test_data <- list("train_X" = train_X, "train_y" = train_y,
                          "test_X" = test_X, "test_y" = test_y)
  return(train_test_data)
}
