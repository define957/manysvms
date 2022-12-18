#' Twin SVM for Multi-classification by Using Ones versus Rest Strategy
#' An R implementation of \code{twinsvm} (Cross Validation)
#' @author Zhang Jiaqi
#' @param X,y dataset and label
#' @param K `K = 10` means K-fold cross validation
#' @param C plenty term.
#' @param kernel kernel function
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param reg regularization term to take care of problems due to ill-conditioning in dual problem.
#' @param kernel_rect set kernel size. \code{0<= kernel_rect <= 1}
#' @param shuffle if set \code{shuffer==TRUE}, This function will shuffle the dataset.
#' @param seed random seed for \code{shuffle} option.
#' @param tol the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param rcpp speed up your code with Rcpp, default \code{rcpp = TRUE}.
#' @param threads.num The number of threads used for parallel execution.
#' @import foreach
#' @import doParallel
#' @import doSNOW
#' @import stats
#' @export
#' @examples
#' library(manysvms)
#' data("iris")
#'
#' X <- iris[, 1:4]
#' y <- iris[, 5]
#' cv.twinsvm_ovr(X, y, K = 10, kernel = 'rbf',
#'                gamma = 1/8, seed = 1234, threads.num = 2)

cv.twinsvm_ovr <- function(X, y , K = 5, C = 1,
                           kernel = c('linear', 'rbf', 'poly'),
                           gamma = 1/ncol(X), degree = 3, coef0 = 0,
                           reg = 1e-7, kernel_rect = 1,
                           shuffle = TRUE, seed = NULL,
                           tol = 1e-5, max.steps = 200, rcpp = TRUE,
                           threads.num = parallel::detectCores() - 1){
  X <- as.matrix(X)
  y <- as.matrix(y)

  param <- expand.grid(C, gamma, degree, coef0)
  m <- nrow(X)
  if (shuffle == TRUE) {
    if (is.null(seed) == FALSE) {
      set.seed(seed)
    }
    new_idx <- sample(m)

  }else{
    new_idx <- 1:m
  }
  v_size <- m %/% K
  j = 1

  cl <- parallel::makeCluster(threads.num)
  pb <- utils::txtProgressBar(max = nrow(param), style = 3)
  progress <- function(n){utils::setTxtProgressBar(pb, n)}
  opts <- list(progress = progress)
  doSNOW::registerDoSNOW(cl)
  res <- foreach::foreach(j = 1:nrow(param), .combine = rbind,
                          .packages = c('manysvms', 'Rcpp'),
                          .options.snow = opts) %dopar% {
    indx_cv <- 1
    accuracy_list <- c()
    for (i in 1:K) {
      new_idx_k <- new_idx[indx_cv:(indx_cv + v_size - 1)] #get test dataset
      indx_cv <- indx_cv + v_size
      test_X <- X[new_idx_k, ]
      train_X <- X[-new_idx_k, ]
      test_y <- y[new_idx_k]
      train_y <- y[-new_idx_k]
      twinsvm_ovr <- twinsvm_ovr(train_X, train_y,
                                 Ck = param[j, 1]*rep(1, length(unique(train_y))),
                                 kernel = kernel, reg = reg,
                                 gamma = param[j, 2],
                                 degree = param[j, 3],
                                 coef0 = param[j, 4],
                                 max.steps = max.steps,
                                 tol = tol,
                                 rcpp = rcpp)
      pred <- predict(twinsvm_ovr, test_X, test_y)
      accuracy_list <- append(accuracy_list, pred$accuracy)
    }
    avg_acc <- mean(accuracy_list)
    sd_acc <- sd(accuracy_list)

    cv_list <- list("accuracy" = avg_acc,
                    "sd_acc" = sd_acc)
    cv_list
  }
  close(pb)
  parallel::stopCluster(cl)
  res <- matrix(res, ncol = 2)

  max_idx <- which.max(res[ ,1])
  call <- match.call()
  cat("\nCall:", deparse(call, 0.8 * getOption("width")), "\n", sep = "\n")
  cat("Total Parameters:", nrow(param), "\n")
  cat("Best Parameters :",
      "C = ", param[max_idx, 1],
      "\n",
      "gamma = ", param[max_idx, 2],
      "degree = ",param[max_idx, 3],
      "coef0 =", param[max_idx, 4],
      "\n")
  cat("Accuracy :", as.numeric(res[max_idx, 1]),
      "Sd :", as.numeric(res[max_idx, 2]), "\n")
  return(res)
}
