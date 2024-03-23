metrics_check_cv <- function(metrics) {
  if (is.list(metrics) == F) {
    metrics <- list(metrics)
    names(metrics) <- paste("metric", length(metrics), sep = "")
  }
  return(metrics)
}

metrics_params_check_cv <- function(num_metrics, metrics_params) {
  if (is.null(metrics_params)) {
    metrics_params <- rep(list(NULL), num_metrics)
  } else {
    len_params <- length(metrics_params)
    metrics_params <- append(metrics_params,
                             rep(NULL, num_metrics - len_params))
  }
  return(metrics_params)
}

metric_evaluate <- function(metric_func, y, y_hat, metric_params) {
  metric_params <- append(list("y" = y, "y_hat" = y_hat), metric_params)
  evaluate_res <- do.call("metric_func", metric_params)
  return(evaluate_res)
}

predict_model <- function(model_res, X_test, y_test,
                          predict_params, predict_func) {
  start_predict <- Sys.time()
  y_test_hat <- do.call("predict_func", append(list(model_res, X_test),
                                               predict_params))
  end_predict <- Sys.time()
  predict_time <- end_predict - start_predict
  return(list("y_test_hat" = y_test_hat, "predict_time" = predict_time))
}

#' K-Fold Cross Validation
#'
#' @author Zhang Jiaqi.
#' @param model your model.
#' @param X,y dataset and label.
#' @param K number of folds.
#' @param metrics this parameter receive a metric function.
#' @param predict_func this parameter receive a function for predict.
#' @param pipeline preprocessing pipline.
#' @param metrics_params set parameters for each metrics (need a list).
#' @param predict_params set parameters for each predict method (need a list).
#' @param model_settings set parameters for model (need a list).
#' @return return a metric matrix
#' @export
cross_validation <- function(model, X, y, K = 5, metrics, predict_func = predict,
                             pipeline = NULL,
                             metrics_params = NULL, predict_params = NULL,
                             model_settings = NULL) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  n <- nrow(X)
  metrics <- metrics_check_cv(metrics)
  num_metric <- length(metrics)
  metrics_params_check_cv(num_metric, metrics_params)
  metric_mat <- matrix(0, num_metric, K)
  index <- sort(rep(1:K, length.out = n))
  for (i in 1:K) {
    idx <- which(index == i)
    X_test <- X[idx, ]
    y_test <- y[idx]
    if (K == 1) {
      X_train <- X_test
      y_train <- y_test
    }else {
      X_train <- X[-idx, ]
      y_train <- y[-idx]
    }
    if (is.null(pipeline) == F) {
      for (pipi in 1:length(pipeline)) {
        pip_temp <- pipeline[[pipi]](X_train)
        X_train <- trans(pip_temp, X_train)
        X_test <- trans(pip_temp, X_test)
      }
    }
    model_res <- do.call("model", append(list("X" = X_train, "y" = y_train),
                                         model_settings))
    predict_res <- predict_model(model_res, X_test, y_test,
                                 predict_params, predict_func)
    for (j in 1:num_metric) {
      metric_j <- metrics[[j]]
      metric_mat[j, i] <- metric_evaluate(metric_j,
                                          y_test, predict_res$y_test_hat,
                                          metrics_params[[j]])
    }
  }
  rownames(metric_mat) <- names(metrics)
  return(metric_mat)
}


#' Grid Search and Cross Validation
#'
#' @author Zhang Jiaqi.
#' @param model your model.
#' @param X,y dataset and label.
#' @param K number of folds.
#' @param metrics this parameter receive a metric function.
#' @param param_list parameter list.
#' @param predict_func this parameter receive a function for predict.
#' @param pipeline preprocessing pipline.
#' @param metrics_params set parameters for each metrics (need a list).
#' @param predict_params set parameters for each predict method (need a list).
#' @param model_settings set parameters for model (need a list).
#' @param shuffle if set \code{shuffle==TRUE}, This function will shuffle
#'                the dataset.
#' @param seed random seed for \code{shuffle} option.
#' @param threads.num the number of threads used for parallel execution.
#' @return return a metric matrix
#' @import foreach
#' @import doParallel
#' @import doSNOW
#' @import stats
#' @export
grid_search_cv <- function(model, X, y, K = 5, metrics, param_list,
                           predict_func = predict,
                           pipeline = NULL,
                           metrics_params = NULL, predict_params = NULL,
                           model_settings = NULL,
                           shuffle = TRUE, seed = NULL,
                           threads.num = parallel::detectCores() - 1) {
  s <- Sys.time()
  X <- as.matrix(X)
  y <- as.matrix(y)
  if (is.list(metrics) == F) {
    metrics <- list(metrics)
    names(metrics) <- paste("metric", length(metrics), sep = "")
  }
  n <- nrow(X)
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
    params_cv <- list("model" = model,
                      "X" = X, "y" = y, "K" = K,
                      "metrics" = metrics,
                      "predict_func" =  predict_func,
                      "pipeline" = pipeline,
                      "metrics_params" = metrics_params,
                      "model_settings" = append(model_settings, as.list(temp))
                       )
    cv_res <- do.call("cross_validation", params_cv)
    cv_res <- rbind(c(apply(cv_res, 1, mean), apply(cv_res, 1, sd)))
  }
  close(pb)
  parallel::stopCluster(cl)
  cat("\n")
  num_metrics <- length(metrics)
  name_matrics <- names(metrics)
  colnames(cv_res)[(num_metrics + 1):(2*num_metrics)] <- paste(name_matrics, "- sd")
  e <- Sys.time()
  idx_max <- apply(as.matrix(cv_res[,1:num_metrics]), 2, which.max)
  idx_min <- apply(as.matrix(cv_res[,1:num_metrics]), 2, which.min)
  score_mat <- matrix(0, 2, 2*num_metrics)
  rownames(score_mat) <- c("max", "min")
  colnames(score_mat) <- c(name_matrics, paste(name_matrics, "- sd"))
  for (i in 1:num_metrics) {
    score_mat[1, i] <- cv_res[idx_max[i], i]
    score_mat[2, i] <- cv_res[idx_min[i], i]
    score_mat[1, num_metrics + i] <- cv_res[idx_max[i], num_metrics + i]
    score_mat[2, num_metrics + i] <- cv_res[idx_min[i], num_metrics + i]
  }
  cv_res <- cbind(cv_res, param_grid)
  cv_model <- list("results" = cv_res,
                   "idx_max" = idx_max,
                   "idx_min" = idx_min,
                   "num.parameters" = n_param,
                   "K" = K,
                   "time" = e - s,
                   "score_mat" = score_mat,
                   "param_grid" = param_grid
                   )
  class(cv_model) <- "cv_model"
  return(cv_model)
}


#' Print Method for Grid-Search and Cross Validation Results
#'
#' @param x object of class \code{eps.svr}.
#' @param ... unsed argument.
#' @export
print.cv_model <- function(x, ...) {
  cat("Results of Grid Search and Cross Validation\n\n")
  cat("Number of Fold", x$K, "\n")
  cat("Total Parameters:", x$num.parameters, "\n\n")
  cat("Time Cost:\n")
  print(x$time)
  cat("Summary of Metrics\n\n")
  print(x$score_mat)
}


#' Grid Search and Cross Validation with Noisy (Simulation Only)
#'
#' @author Zhang Jiaqi.
#' @param model your model.
#' @param X,y dataset and label.
#' @param y_noisy label (contains label noise)
#' @param K number of folds.
#' @param metrics this parameter receive a metric function.
#' @param param_list parameter list.
#' @param predict_func this parameter receive a function for predict.
#' @param pipeline preprocessing pipline.
#' @param metrics_params set parameter for each metrics (need a list).
#' @param predict_params set parameters for each predict method (need a list).
#' @param model_settings set parameters for model (need a list).
#' @param shuffle if set \code{shuffle==TRUE}, This function will shuffle
#'                the dataset.
#' @param seed random seed for \code{shuffle} option.
#' @param threads.num the number of threads used for parallel execution.
#' @return return a metric matrix
#' @import foreach
#' @import doParallel
#' @import doSNOW
#' @import stats
#' @export
grid_search_cv_noisy <- function(model, X, y, y_noisy, K = 5, metrics, param_list,
                                 predict_func = predict,
                                 pipeline = NULL,
                                 metrics_params = NULL, predict_params = NULL,
                                 model_settings = NULL,
                                 shuffle = TRUE, seed = NULL,
                                 threads.num = parallel::detectCores() - 1) {
  s <- Sys.time()
  X <- as.matrix(X)
  y <- as.matrix(y)
  if (is.list(metrics) == F) {
    metrics <- list(metrics)
    names(metrics) <- paste("metric", length(metrics), sep = "")
  }
  n <- nrow(X)
  if (is.null(seed) == FALSE) {
    set.seed(seed)
  }
  if (shuffle == TRUE) {
    idx <- sample(n)
    X <- X[idx, ]
    y <- y[idx]
    y_noisy <- y_noisy[idx]
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
    params_cv <- list("model" = model,
                      "X" = X, "y" = y, "y_noisy" = y_noisy, "K" = K,
                      "metrics" = metrics,
                      "predict_func" =  predict_func,
                      "pipeline" = pipeline,
                      "metrics_params" = metrics_params,
                      "model_settings" = append(model_settings, as.list(temp))
                      )
    cv_res <- do.call("cross_validation_noisy", params_cv)
    cv_res <- rbind(c(apply(cv_res, 1, mean), apply(cv_res, 1, sd)))
  }
  parallel::stopCluster(cl)
  cat("\n")
  num_metrics <- length(metrics)
  name_matrics <- names(metrics)
  colnames(cv_res)[(num_metrics + 1):(2*num_metrics)] <- paste(name_matrics, "- sd")
  e <- Sys.time()
  idx_max <- apply(as.matrix(cv_res[,1:num_metrics]), 2, which.max)
  idx_min <- apply(as.matrix(cv_res[,1:num_metrics]), 2, which.min)
  score_mat <- matrix(0, 2, 2*num_metrics)
  rownames(score_mat) <- c("max", "min")
  colnames(score_mat) <- c(name_matrics, paste(name_matrics, "- sd"))
  for (i in 1:num_metrics) {
    score_mat[1, i] <- cv_res[idx_max[i], i]
    score_mat[2, i] <- cv_res[idx_min[i], i]
    score_mat[1, num_metrics + i] <- cv_res[idx_max[i], num_metrics + i]
    score_mat[2, num_metrics + i] <- cv_res[idx_min[i], num_metrics + i]
  }
  cv_res <- cbind(cv_res, param_grid)
  cv_model <- list("results" = cv_res,
                   "idx_max" = idx_max,
                   "idx_min" = idx_min,
                   "num.parameters" = n_param,
                   "K" = K,
                   "time" = e - s,
                   "score_mat" = score_mat,
                   "param_grid" = param_grid
                   )
  class(cv_model) <- "cv_model"
  return(cv_model)
}

#' K-Fold Cross Validation with Noisy (Simulation Only)
#'
#' \code{cross_validation_noisy} function use noisy data for training,
#' then calculates the average and standard deviation of your metric
#' using clean samples
#'
#' @author Zhang Jiaqi.
#' @param model your model.
#' @param X,y dataset and label.
#' @param y_noisy label with label noise.
#' @param K number of folds.
#' @param metrics this parameter receive a metric function.
#' @param predict_func this parameter receive a function for predict.
#' @param pipeline preprocessing pipline.
#' @param metrics_params set parameters for each metrics (need a list).
#' @param predict_params set parameters for each predict method (need a list).
#' @param model_settings set parameters for model (need a list).
#' @return return a metric matrix
#' @export
cross_validation_noisy <- function(model, X, y, y_noisy, K = 5, metrics,
                                   predict_func = predict,
                                   pipeline = NULL,
                                   metrics_params = NULL, predict_params = NULL,
                                   model_settings = NULL) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  y_noisy <- as.matrix(y_noisy)
  metrics <- metrics_check_cv(metrics)
  num_metric <- length(metrics)
  metrics_params_check_cv(num_metric, metrics_params)
  n <- nrow(X)
  num_metric <- length(metrics)
  metric_mat <- matrix(0, num_metric, K)
  index <- sort(rep(1:K, length.out = n))
  for (i in 1:K) {
    idx <- which(index == i)
    X_test <- X[idx, ]
    y_test <- y[idx]
    if (K == 1) {
      X_train <- X_test
      y_train <- y_test
    }else{
      X_train <- X[-idx, ]
      y_train <- y_noisy[-idx]
    }
    if (is.null(pipeline) == F) {
      for (pipi in 1:length(pipeline)) {
        pip_temp <- pipeline[[pipi]](X_train)
        X_train <- trans(pip_temp, X_train)
        X_test <- trans(pip_temp, X_test)
      }
    }
    model_res <- do.call("model", append(list("X" = X_train, "y" = y_train),
                                         model_settings))
    predict_res <- predict_model(model_res, X_test, y_test,
                                 predict_params, predict_func)
    for (j in 1:num_metric) {
      metric_j <- metrics[[j]]
      metric_mat[j, i] <- metric_evaluate(metric_j,
                                          y_test, predict_res$y_test_hat,
                                          metrics_params[[j]])
    }
  }
  return(metric_mat)
}
