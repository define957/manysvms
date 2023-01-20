#' SVM for Multi-classification by Using Ones versus Rest Strategy
#'
#' \code{svm_ovr} is an R implementation of support vector machine
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param C plenty term.
#' @param kernel kernel function. The definitions of various kernel functions are as follows:
#' \describe{
#'     \item{linear:}{\eqn{u'v}{u'*v}}
#'     \item{poly:}{\eqn{(\gamma u'v + coef0)^{degree}}{(gamma*u'*v + coef0)^degree}}
#'     \item{rbf:}{\eqn{e^{(-\gamma |u-v|^2)}}{exp(-gamma*|u-v|^2)}}
#' }
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param tol the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param rcpp speed up your code with Rcpp, default \code{rcpp = TRUE}.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @return return \code{svm_ovr} object.
#' @export

svm_ovr <- function(X, y,
                    C = 1.0,
                    kernel = c('linear', 'rbf', 'poly'),
                    gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                    tol = 1e-5, max.steps = 200,
                    rcpp = TRUE, fit_intercept = TRUE) {
  kernel <- match.arg(kernel)

  X <- as.matrix(X)
  y <- as.matrix(y)

  m <- nrow(X)
  n <- ncol(X)

  X <- as.matrix(X)
  y <- as.matrix(y)

  class_set <- unique(as.matrix(y))
  class_num <- length(class_set)

  KernelX <- kernel_function(X, X,
                             kernel.type = kernel,
                             gamma = gamma, degree = degree, coef0 = coef0,
                             rcpp = rcpp)
  y_temp <- rep(0, m)
  coef_list <- matrix(0, nrow = m, ncol = class_num)
  intercept_list <- matrix(0, ncol = class_num)

  for (k in 1:class_num) {
    idx_pos <- which(y == class_set[k])

    y_temp[idx_pos] <- 1
    y_temp[-idx_pos] <- -1

    y_temp_mat <- diag(y_temp)

    H <- y_temp_mat %*% KernelX %*% y_temp_mat
    e <- matrix(1, nrow = m)
    lb <- matrix(0, nrow = m)
    ub <- matrix(C, nrow = m)

    x <- clip_dcd_optimizer(H, e, lb, ub, tol, max.steps, rcpp)$x

    alphas <- as.matrix(x)
    coef_list[, k] <- y_temp_mat %*% alphas
    if (fit_intercept == TRUE) {
      temp <- y_temp - KernelX %*% coef_list[, k]
      intercept_list[k] <- mean(temp[ub > alphas && alphas > lb ])
    }
  }

  svm_model <- list('X' = X, 'y' = y,
                    'class_set' = class_set,'class_num' = class_num,
                    'coef' = coef_list,
                    'kernel' = kernel,
                    'gamma' = gamma,
                    'degree' = degree,
                    'coef0' = coef0,
                    'Rcpp' = rcpp,
                    'fit_intercept' = fit_intercept
                    )
  if (fit_intercept == TRUE) {
    svm_model$intercept <- intercept_list
  }
  class(svm_model) <- "svm_ovr"
  return(svm_model)
}


#' Predict Method for Support Vector Machine
#'
#' @author Zhang Jiaqi
#' @param object A fitted object of class inheriting from \code{svm_ovr}.
#' @param X A new data frame for predicting.
#' @param y A label data frame corresponding to X.
#' @param ... unused parameter.
#' @importFrom stats predict
#' @export

predict.svm_ovr <- function(object, X, y, ...) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  m <- nrow(X)
  KernelX <- kernel_function(X, object$X,
                             kernel.type = object$kernel,
                             gamma = object$gamma, degree = object$degree,
                             coef0 = object$coef0,
                             rcpp = object$Rcpp)
  if (object$fit_intercept == TRUE) {
    intercept_mat <- as.matrix(rep(object$intercept, m))
    dim(intercept_mat) <- c(object$class_num, m)
    intercept_mat <- t(intercept_mat)
    pred_value <- KernelX %*% object$coef + intercept_mat
  } else {
    pred_value <- KernelX %*% object$coef
  }
  vote_mat <- matrix(0, nrow = m, ncol = object$class_num)
  for (i in 1:object$class_num) {
    idx <- which(pred_value[, i] > 0)
    vote_mat[idx, i] <- vote_mat[idx, i] + 1
  }
  pred <- rep(0, m)
  for (i in 1:m) {
    idx <- which.max(vote_mat[i, ])
    pred[i] <- object$class_set[idx]
  }
  if (is.null(y) == FALSE) {
    y <- as.matrix(y)
    acc <- sum(pred == y)/length(y)
    cat('kernel type :', object$kernel, '\n')
    cat(paste("total accuracy :", acc*100, '% \n'))
    predlist <- list("accuracy" = acc,
                     'vote_mat' = vote_mat,
                     'predict' = pred,
                     'pred_value' = pred_value)
  }else{
    predlist <- list('vote_mat' = vote_mat,
                     'predict' = pred,
                     'pred_value' = pred_value)
  }
  return(predlist)
}



#' Computes K-fold Cross-Validation Accuracy for svm_ovr
#'
#' @author Zhang Jiaqi
#' @param X,y dataset and label.
#' @param C plenty term.
#' @param K Number of folds.
#' @param kernel kernel type. Default value \code{kernel = 'linear'}.
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param tol the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param rcpp speed up your code with Rcpp, default \code{rcpp = TRUE}.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param shuffle if set \code{shuffle==TRUE}, This function will shuffle the dataset.
#' @param seed random seed for \code{shuffer} option.
#' @param threads.num The number of threads used for parallel execution.
#' @import foreach
#' @import doParallel
#' @import doSNOW
#' @import stats
#' @export

cv.svm_ovr <- function(X, y , K = 5, C = 1,
                       kernel = c('linear', 'rbf', 'poly'),
                       gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                       tol = 1e-5,
                       max.steps = 200, rcpp = TRUE, fit_intercept = TRUE,
                       shuffle = TRUE, seed = NULL,
                       threads.num = parallel::detectCores() - 1) {

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

  cl <- parallel::makeCluster(threads.num)
  pb <- utils::txtProgressBar(max = nrow(param), style = 3)
  progress <- function(n){utils::setTxtProgressBar(pb, n)}
  opts <- list(progress = progress)
  doSNOW::registerDoSNOW(cl)
  j <- 1
  res <- foreach::foreach(j = 1:nrow(param), .combine = rbind,
                          .packages = c('manysvms', 'Rcpp'),
                          .options.snow = opts) %dopar% {
    indx_cv <- 1
    accuracy_list <- rep(0, K)
    for (i in 1:K) {
      new_idx_k <- new_idx[indx_cv:(indx_cv + v_size - 1)] #get test dataset
      indx_cv <- indx_cv + v_size
      test_X <- X[new_idx_k, ]
      train_X <- X[-new_idx_k, ]
      test_y <- y[new_idx_k]
      train_y <- y[-new_idx_k]
      mbsvm_model <- svm_ovr(train_X, train_y,
                             C = param[j, 1],
                             kernel = kernel, tol = tol,
                             gamma = param[j, 2], degree = param[j, 3],
                             coef0 = param[j, 4],
                             max.steps = max.steps,
                             fit_intercept = fit_intercept)
      pred <- predict(mbsvm_model, test_X, test_y)
      accuracy_list[i] <- pred$accuracy
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
      "C = ", param[max_idx, 1], "\n")
  cat("gamma = ", param[max_idx, 2],
      "degree = ",param[max_idx, 3],
      "coef0 =", param[max_idx, 4],
      "\n")
  cat("Accuracy :", as.numeric(res[max_idx, 1]),
      "Sd :", as.numeric(res[max_idx, 2]), "\n")
  return(res)
}
