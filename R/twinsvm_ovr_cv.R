#' Twin SVM for Multi-classification by Using Ones versus Rest Strategy
#' An R implementation of \code{twinsvm} (Cross Validation)
#' @author Zhang Jiaqi
#' @param X,y dataset and label
#' @param K `K = 10` means K-fold cross validation
#' @param Ck plenty term vector
#' @param kernel kernel function
#' @param reg regularization tern
#' @param gamma rbf kernel parameter
#' @param shuffer `shuffer = TRUE` means shuffer the dataset
#' @param seed is shuffer random seed
#' @return average accuracy in K-fold cross validation
#' @export
#' @examples
#' library(manysvms)
#' data("iris")
#'
#' X <- iris[, 1:4]
#' y <- iris[, 5]
#' cv.twinsvm_ovr(X, y, K = 10, kernel = 'rbf', gamma = 1/8, seed = 1234)

cv.twinsvm_ovr <- function(X, y , K = 5, Ck = rep(1, length(unique(y))),
                           kernel = c('linear', 'rbf', 'poly'),
                           reg = 1e-3, gamma = 1/ncol(X),
                           shuffer = TRUE, seed = NULL){

  m <- nrow(X)
  if(shuffer == TRUE){
    if(is.null(seed) == FALSE){
      set.seed(seed)
    }
    new_idx <- sample(m)

  }else{
    new_idx <- 1:m
  }
  v_size <- m %/% K
  indx_cv <- 1
  accuracy_list <- c()
  for(i in 1:K){
    new_idx_k <- new_idx[indx_cv:(indx_cv+v_size - 1)] #get test dataset
    indx_cv <- indx_cv + v_size
    test_X <- X[new_idx_k, ]
    train_X <- X[-new_idx_k, ]
    test_y <- y[new_idx_k]
    train_y <- y[-new_idx_k]
    twinsvm_ovr <- twinsvm_ovr(train_X, train_y, Ck = Ck, kernel = kernel, reg = reg, gamma = gamma)
    pred <- predict(twinsvm_ovr, test_X, test_y)
    accuracy_list <- append(accuracy_list, pred$accuracy)
  }
  avg_acc <- mean(accuracy_list)
  cat('average accuracy in ',K, 'fold cross validation :', 100*avg_acc, '%\n')
}
