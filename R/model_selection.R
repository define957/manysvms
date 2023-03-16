#' K-Fold Cross Validation
#'
#' @author Zhang Jiaqi.
#' @param model your model.
#' @param X,y dataset and label.
#' @param K number of folds.
#' @param metric this parameter receive a metric function.
#' @param predict_func this parameter receive a function for predict.
#' @param shuffle if set \code{shuffle==TRUE}, This function will shuffle the dataset.
#' @param seed random seed for \code{shuffle} option.
#' @param ... additional parameters for your model.
#' @return return a metric matrix
#' @export
cross_validation <- function(model, X, y, K = 5, metric, predict_func = predict,
                             ...) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  n <- nrow(X)
  index = c(0:K)*n/K

  metric_mat <- matrix(0, nrow = 1, ncol = K)
  for (i in 1:K) {
    X_test <- X[c(index[i]:index[i + 1]), ]
    y_test <- y[c(index[i]:index[i + 1])]
    X_train <- X[-c(index[i]:index[i + 1]), ]
    y_train <- y[-c(index[i]:index[i + 1])]
    model_res <- do.call("model", list("X" = X_train, "y" = y_train, ...))
    y_test_hat <- predict_func(model_res, X_test)
    metric_mat[i] <- metric(y_test, y_test_hat)
  }
  return(metric_mat)
}

#' Grid Search and Cross Validation
#'
#' @author Zhang Jiaqi.
#' @param model your model.
#' @param X,y dataset and label.
#' @param K number of folds.
#' @param metric this parameter receive a metric function.
#' @param param_list parameter list.
#' @param predict_func this parameter receive a function for predict.
#' @param shuffle if set \code{shuffle==TRUE}, This function will shuffle
#'                the dataset.
#' @param seed random seed for \code{shuffle} option.
#' @param threads.num the number of threads used for parallel execution.
#' @param ... additional parameters for your model.
#' @return return a metric matrix
#' @import foreach
#' @import doParallel
#' @import doSNOW
#' @import stats
#' @export
grid_search_cv <- function(model, X, y, K = 5, metric, param_list,
                           predict_func = predict,
                           shuffle = TRUE, seed = NULL,
                           threads.num = parallel::detectCores() - 1, ...) {
  s <- Sys.time()
  X <- as.matrix(X)
  y <- as.matrix(y)
  if (is.null(seed) == FALSE) {
    set.seed(seed)
  }
  if (shuffle == TRUE) {
    idx <- sample(n)
    X <- X[idx, ]
    y <- y[idx]
  }
  param_grid <- expand.grid(param_list)
  n_param <- nrow(param_grid)
  param_names <- colnames(param_grid)
  cl <- parallel::makeCluster(threads.num)
  pb <- utils::txtProgressBar(max = n_param, style = 3)
  progress <- function(n){utils::setTxtProgressBar(pb, n)}
  opts <- list(progress = progress)
  doSNOW::registerDoSNOW(cl)
  i <- 1
  cv_res <- foreach::foreach(i = 1:n_param, .combine = rbind,
                             .packages = c('manysvms', 'Rcpp'),
                             .options.snow = opts) %dopar% {
    temp <- data.frame(param_grid[i, ])
    colnames(temp) <- param_names
    params_cv <- append(list("model" = model,
                              "X" = X, "y" = y, "K" = K,
                              "metric" = metric,
                              "predict_func" =  predict_func,
                               ...),
                               temp)
    cv_res <- do.call("cross_validation", params_cv)
  }
  close(pb)
  parallel::stopCluster(cl)
  cv_res = cbind(apply(cv_res, 1, mean), apply(cv_res, 1, sd),
                 cv_res,param_grid)
  colnames(cv_res) = c("avg", "sd", as.character(c(1:K)), names(param_list))
  e <- Sys.time()
  idx.best <- which.max(cv_res$avg)
  cv_model <- list("results" = cv_res,
                   "idx.best" = idx.best,
                   "num.parameters" = n_param,
                   "K" = K,
                   "time" = e - s)
  class(cv_model) <- "cv_model"
  return(cv_model)
}

#' Print Method for Grid-Search and Cross Validation Results
#'
#' @param x object of class \code{eps.svr}.
#' @param ... unsed argument.
#' @export
print.cv_model <- function(x, ...) {
  cat("Number of Fold", x$K, "\n")
  cat("Total Parameters:", x$num.parameters, "\n")
  cat("Time Cost:", x$time, "\n")
  cat("Best Avg.:", x$results[x$idx.best, 1], "\n")
  cat("Best Sd:", x$results[x$idx.best, 2], "\n")
  cat("Best Parameter:", "\n")
  print(x$results[x$idx.best, (3 + x$K):ncol(x$results)])
}
