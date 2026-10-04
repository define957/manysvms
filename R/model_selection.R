metrics_check_cv <- function(metrics) {
  if (is.list(metrics) == F) {
    metrics <- list(metrics)
    names(metrics) <- paste("metric", length(metrics), sep = "")
  }
  return(metrics)
}

metrics_params_check_cv <- function(num_metrics, metrics_params) {
  if (is.null(metrics_params)) {
    metrics_params <- vector("list", num_metrics)
  } else {
    if (!is.list(metrics_params)) {
      stop("'metrics_params' must be a list.")
    }
    len_params <- length(metrics_params)
    if (len_params != num_metrics) {
      stop(sprintf("Length of 'metrics_params' (%d) does not match number of metrics (%d).", len_params, num_metrics))
    }
  }
  return(metrics_params)
}

#' Generate K Folds Indices for Cross Validation
#'
#' Split sample indices into \code{K} folds, optionally with shuffling
#' and/or stratification by a label vector. The result can be passed to
#' \code{cross_validation} and related functions via their \code{folds}
#' argument.
#'
#' @author Zhang Jiaqi.
#' @param num_data total number of samples (a positive integer).
#' @param K number of folds, must be a positive integer no larger than
#'          \code{num_data}.
#' @param shuffle if \code{TRUE}, shuffle samples before splitting.
#' @param stratified if \code{TRUE}, split each class of \code{stratify_by}
#'                   separately so that every fold keeps (approximately)
#'                   the same class proportions as the full dataset.
#' @param stratify_by a vector of length \code{num_data} used for
#'                    stratification (e.g. the label vector). Must be
#'                    provided when \code{stratified = TRUE}; ignored
#'                    otherwise.
#' @param seed random seed for shuffling.
#' @return a list of length \code{K}; each element is an integer vector
#'         containing the test-set row indices of one fold. The folds
#'         partition \code{1:num_data} exactly.
#' @export
make_k_folds <- function(
    num_data,
    K = 5,
    shuffle = FALSE,
    stratified = FALSE,
    stratify_by = NULL,
    seed = NULL
) {
  if (!is.numeric(num_data) || length(num_data) != 1 ||
      is.na(num_data) || num_data %% 1 != 0 || num_data < 1) {
    stop("'num_data' must be a positive integer.")
  }
  if (!is.numeric(K) || length(K) != 1 ||
      is.na(K) || K %% 1 != 0 || K < 1) {
    stop("'K' must be a positive integer.")
  }
  if (K > num_data) {
    stop(sprintf("Cannot have K = %d with num_data = %d.", K, num_data))
  }
  if (stratified == TRUE && is.null(stratify_by)) {
    stop("`stratify_by` must be provided when `stratified = TRUE`.")
  }
  if (!is.null(stratify_by) && any(is.na(stratify_by))) {
    stop("`stratify_by` contains NA.")
  }
  if (!is.null(stratify_by)) {
    stratify_by <- as.vector(stratify_by)
    if (num_data != length(stratify_by)) {
      stop("length of 'stratify_by' does not match 'num_data'.")
    }
  }
  if (is.null(seed) == FALSE) {
    set.seed(seed)
  }
  perm <- if (shuffle) {
    sample.int(num_data)
  } else {
    seq_len(num_data)
  }

  index <- integer(num_data)
  if (stratified) {
    for (cls in unique(stratify_by)) {
      idx_cls <- perm[stratify_by[perm] == cls]
      index[idx_cls] <- sort(rep(1:K, length.out = length(idx_cls)))
    }
  } else {
    index[perm] <- sort(rep(1:K, length.out = num_data))
  }
  folds <- lapply(1:K, function(i) { which(index == i) })
  return(folds)
}

folds_check_cv <- function(folds, n) {
  if (is.numeric(folds) && length(folds) == 1) {
    index <- sort(rep(1:folds, length.out = n))
    folds <- lapply(1:folds, function(i) which(index == i))
  }
  if (is.list(folds) == FALSE) {
    stop("'folds' must be a positive integer or a list of index vectors.")
  }
  K <- length(folds)
  if (K < 1) {
    stop("'folds' must contain at least one fold.")
  }
  for (i in 1:K) {
    if (is.numeric(folds[[i]]) == FALSE || any(folds[[i]] %% 1 != 0)) {
      stop("Each element of 'folds' must be an integer vector of row indices.")
    }
    if (any(folds[[i]] < 1 | folds[[i]] > n)) {
      stop(sprintf("Indices in folds[[%d]] out of range [1, %d].", i, n))
    }
  }
  if (length(unique(unlist(folds))) != n) {
    stop("The folds must cover each sample exactly once.")
  }
  return(folds)
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
#' @param folds a positive integer indicating the number of folds (sequential
#'              split, compatible with the old \code{K} argument) or a list of
#'              index vectors, where each element contains the test-set row
#'              indices of one fold.
#' @param metrics this parameter receive a metric function.
#' @param predict_func this parameter receive a function for predict.
#' @param pipeline preprocessing pipline.
#' @param metrics_params set parameters for each metrics (need a list).
#' @param predict_params set parameters for each predict method (need a list).
#' @param model_settings set parameters for model (need a list).
#' @param transy apply transforms defined in `pipeline` on y, default FALSE.
#' @param model_seed random_seed for model.
#' @return return a metric matrix
#' @export
cross_validation <- function(model, X, y, folds = 5, metrics, predict_func = predict,
                             pipeline = NULL,
                             metrics_params = NULL, predict_params = NULL,
                             model_settings = NULL, transy = F,
                             model_seed = NULL) {
  if (is.null(model_seed) == FALSE) {
    set.seed(model_seed)
  }
  X <- as.matrix(X)
  y <- as.matrix(y)
  n <- nrow(X)
  folds <- folds_check_cv(folds, n)
  K <- length(folds)
  metrics <- metrics_check_cv(metrics)
  num_metric <- length(metrics)
  metrics_params <- metrics_params_check_cv(num_metric, metrics_params)
  metric_mat <- matrix(0, num_metric, K)
  for (i in 1:K) {
    idx <- folds[[i]]
    X_test <- X[idx, , drop = FALSE]
    y_test <- y[idx]
    if (K == 1) {
      X_train <- X_test
      y_train <- y_test
    }else {
      X_train <- X[-idx, , drop = FALSE]
      y_train <- y[-idx]
    }
    if (is.null(pipeline) == F) {
      for (pipi in 1:length(pipeline)) {
        pip_temp <- pipeline[[pipi]](X_train)
        X_train <- trans(pip_temp, X_train)
        X_test <- trans(pip_temp, X_test)
        if (transy == T) {
          pip_temp <- pipeline[[pipi]](y_train)
          y_train <- trans(pip_temp, y_train)
          y_test <- trans(pip_temp, y_test)
        }
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
#' @param folds a positive integer indicating the number of folds (sequential
#'              split, compatible with the old \code{K} argument) or a list of
#'              index vectors, where each element contains the test-set row
#'              indices of one fold.
#' @param metrics this parameter receive a metric function.
#' @param param_list parameter list.
#' @param predict_func this parameter receive a function for predict.
#' @param pipeline preprocessing pipline.
#' @param metrics_params set parameters for each metrics (need a list).
#' @param predict_params set parameters for each predict method (need a list).
#' @param model_settings set parameters for model (need a list).
#' @param transy apply transforms defined in `pipeline` on y, default FALSE.
#' @param shuffle if set \code{shuffle==TRUE}, This function will shuffle
#'                the dataset (only used when \code{folds} is a number).
#' @param seed random seed for \code{shuffle} option (only used when
#'             \code{folds} is a number).
#' @param model_seed random_seed for model.
#' @param threads.num the number of threads used for parallel execution.
#' @return return a metric matrix
#' @import foreach
#' @import doParallel
#' @import doSNOW
#' @import stats
#' @export
grid_search_cv <- function(model, X, y, folds = 5, metrics, param_list,
                           predict_func = predict,
                           pipeline = NULL,
                           metrics_params = NULL, predict_params = NULL,
                           model_settings = NULL, transy = FALSE,
                           shuffle = TRUE, seed = NULL, model_seed = NULL,
                           threads.num = parallel::detectCores() - 1) {
  s <- Sys.time()
  X <- as.matrix(X)
  y <- as.matrix(y)
  if (is.list(metrics) == F) {
    metrics <- list(metrics)
    names(metrics) <- paste("metric", length(metrics), sep = "")
  }
  n <- nrow(X)
  if (is.numeric(folds) && length(folds) == 1) {
    K <- folds
    if (is.null(seed) == FALSE) {
      set.seed(seed)
    }
    if (shuffle == TRUE) {
      idx <- sample(n)
      X <- X[idx, , drop = FALSE]
      y <- y[idx]
    }
  } else {
    K <- length(folds)
  }
  param_grid <- expand.grid(param_list, stringsAsFactors = FALSE)
  n_param <- nrow(param_grid)
  param_names <- colnames(param_grid)
  cl <- parallel::makeCluster(threads.num)
  op <- options(cli.progress_show_after = 0)
  on.exit(options(op), add = TRUE)
  pb <- cli::cli_progress_bar("Grid search", total = n_param, clear = FALSE)
  # pb <- utils::txtProgressBar(max = n_param, style = 3)
  progress <- function(n){cli::cli_progress_update(id = pb, set = n)}
  # progress <- function(n){utils::setTxtProgressBar(pb, n)}
  opts <- list(progress = progress)
  doSNOW::registerDoSNOW(cl)
  i <- 1
  cv_res <- foreach::foreach(i = 1:n_param, .combine = rbind,
                             .packages = c('manysvms', 'Rcpp'),
                             .options.snow = opts) %dopar% {
    params_to_add <- setNames(
     lapply(param_names, function(nm) param_grid[[nm]][[i]]),
     param_names
    )
    params_cv <- list("model" = model,
                      "X" = X, "y" = y, "folds" = folds,
                      "metrics" = metrics,
                      "predict_func" =  predict_func,
                      "pipeline" = pipeline,
                      "metrics_params" = metrics_params,
                      "model_settings" = append(model_settings, params_to_add),
                      "transy" = transy,
                      "model_seed" = model_seed
                       )
    cv_res <- do.call("cross_validation", params_cv)
    cv_res <- rbind(c(apply(cv_res, 1, mean), apply(cv_res, 1, sd)))
  }
  cli::cli_progress_done(id = pb)
  # close(pb)
  parallel::stopCluster(cl)
  # cat("\n")
  num_metrics <- length(metrics)
  name_matrics <- names(metrics)
  colnames(cv_res)[(num_metrics + 1):(2*num_metrics)] <- paste(name_matrics, "- sd")
  e <- Sys.time()
  # idx_max <- apply(as.matrix(cv_res[,1:num_metrics]), 2, which.max)
  # idx_min <- apply(as.matrix(cv_res[,1:num_metrics]), 2, which.min)

  idx_max <- sapply(seq_len(num_metrics),
                    function(j) select_best_idx(cv_res[, j],
                                                cv_res[, num_metrics + j], "max"))
  idx_min <- sapply(seq_len(num_metrics),
                    function(j) select_best_idx(cv_res[, j],
                                                cv_res[, num_metrics + j], "min"))

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
#' @param folds a positive integer indicating the number of folds (sequential
#'              split, compatible with the old \code{K} argument) or a list of
#'              index vectors, where each element contains the test-set row
#'              indices of one fold.
#' @param metrics this parameter receive a metric function.
#' @param param_list parameter list.
#' @param predict_func this parameter receive a function for predict.
#' @param pipeline preprocessing pipline.
#' @param metrics_params set parameter for each metrics (need a list).
#' @param predict_params set parameters for each predict method (need a list).
#' @param model_settings set parameters for model (need a list).
#' @param transy apply transforms defined in `pipeline` on y, default FALSE.
#' @param shuffle if set \code{shuffle==TRUE}, This function will shuffle
#'                the dataset (only used when \code{folds} is a number).
#' @param seed random seed for \code{shuffle} option (only used when
#'             \code{folds} is a number).
#' @param model_seed random_seed for model.
#' @param threads.num the number of threads used for parallel execution.
#' @return return a metric matrix
#' @import foreach
#' @import doParallel
#' @import doSNOW
#' @import stats
#' @export
grid_search_cv_noisy <- function(model, X, y, y_noisy, folds = 5, metrics, param_list,
                                 predict_func = predict,
                                 pipeline = NULL,
                                 metrics_params = NULL, predict_params = NULL,
                                 model_settings = NULL, transy = FALSE,
                                 shuffle = TRUE, seed = NULL, model_seed = NULL,
                                 threads.num = parallel::detectCores() - 1) {
  s <- Sys.time()
  X <- as.matrix(X)
  y <- as.matrix(y)
  if (is.list(metrics) == F) {
    metrics <- list(metrics)
    names(metrics) <- paste("metric", length(metrics), sep = "")
  }
  n <- nrow(X)
  if (is.numeric(folds) && length(folds) == 1) {
    K <- folds
    if (is.null(seed) == FALSE) {
      set.seed(seed)
    }
    if (shuffle == TRUE) {
      idx <- sample(n)
      X <- X[idx, , drop = FALSE]
      y <- y[idx]
      y_noisy <- y_noisy[idx]
    }
  } else {
    K <- length(folds)
  }
  param_grid <- expand.grid(param_list, stringsAsFactors = FALSE)
  n_param <- nrow(param_grid)
  param_names <- colnames(param_grid)
  cl <- parallel::makeCluster(threads.num)
  op <- options(cli.progress_show_after = 0)
  on.exit(options(op), add = TRUE)
  pb <- cli::cli_progress_bar("Grid search", total = n_param, clear = FALSE)
  # pb <- utils::txtProgressBar(max = n_param, style = 3)
  progress <- function(n){cli::cli_progress_update(id = pb, set = n)}
  # progress <- function(n){utils::setTxtProgressBar(pb, n)}
  opts <- list(progress = progress)
  doSNOW::registerDoSNOW(cl)
  i <- 1
  cv_res <- foreach::foreach(i = 1:n_param, .combine = rbind,
                             .packages = c('manysvms', 'Rcpp'),
                             .options.snow = opts) %dopar% {
    params_to_add <- setNames(
      lapply(param_names, function(nm) param_grid[[nm]][[i]]),
      param_names
    )
    params_cv <- list("model" = model,
                      "X" = X, "y" = y, "y_noisy" = y_noisy, "folds" = folds,
                      "metrics" = metrics,
                      "predict_func" =  predict_func,
                      "pipeline" = pipeline,
                      "metrics_params" = metrics_params,
                      "model_settings" = append(model_settings, params_to_add),
                      "transy" = transy,
                      "model_seed" = model_seed
                      )
    cv_res <- do.call("cross_validation_noisy", params_cv)
    cv_res <- rbind(c(apply(cv_res, 1, mean), apply(cv_res, 1, sd)))
  }
  cli::cli_progress_done(id = pb)
  # close(pb)
  parallel::stopCluster(cl)
  # cat("\n")
  num_metrics <- length(metrics)
  name_matrics <- names(metrics)
  colnames(cv_res)[(num_metrics + 1):(2*num_metrics)] <- paste(name_matrics, "- sd")
  e <- Sys.time()
  # idx_max <- apply(as.matrix(cv_res[,1:num_metrics]), 2, which.max)
  # idx_min <- apply(as.matrix(cv_res[,1:num_metrics]), 2, which.min)

  idx_max <- sapply(seq_len(num_metrics),
                    function(j) select_best_idx(cv_res[, j],
                                                cv_res[, num_metrics + j], "max"))
  idx_min <- sapply(seq_len(num_metrics),
                    function(j) select_best_idx(cv_res[, j],
                                                cv_res[, num_metrics + j], "min"))

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
#' @param folds a positive integer indicating the number of folds (sequential
#'              split, compatible with the old \code{K} argument) or a list of
#'              index vectors, where each element contains the test-set row
#'              indices of one fold.
#' @param metrics this parameter receive a metric function.
#' @param predict_func this parameter receive a function for predict.
#' @param pipeline preprocessing pipline.
#' @param metrics_params set parameters for each metrics (need a list).
#' @param predict_params set parameters for each predict method (need a list).
#' @param model_settings set parameters for model (need a list).
#' @param transy apply transforms defined in `pipeline` on y, default FALSE.
#' @param model_seed random_seed for model.
#' @return return a metric matrix
#' @export
cross_validation_noisy <- function(model, X, y, y_noisy, folds = 5, metrics,
                                   predict_func = predict,
                                   pipeline = NULL,
                                   metrics_params = NULL, predict_params = NULL,
                                   model_settings = NULL, transy = FALSE,
                                   model_seed = NULL
                                   ) {
  if (is.null(model_seed) == FALSE) {
    set.seed(model_seed)
  }
  X <- as.matrix(X)
  y <- as.matrix(y)
  y_noisy <- as.matrix(y_noisy)
  n <- nrow(X)
  folds <- folds_check_cv(folds, n)
  K <- length(folds)
  metrics <- metrics_check_cv(metrics)
  num_metric <- length(metrics)
  metrics_params <- metrics_params_check_cv(num_metric, metrics_params)
  metric_mat <- matrix(0, num_metric, K)
  for (i in 1:K) {
    idx <- folds[[i]]
    X_test <- X[idx, , drop = FALSE]
    y_test <- y[idx]
    if (K == 1) {
      X_train <- X_test
      y_train <- y_test
    }else{
      X_train <- X[-idx, , drop = FALSE]
      y_train <- y_noisy[-idx]
      y_train_clean <- y[-idx]
    }
    if (is.null(pipeline) == F) {
      for (pipi in 1:length(pipeline)) {
        pip_temp <- pipeline[[pipi]](X_train)
        X_train <- trans(pip_temp, X_train)
        X_test <- trans(pip_temp, X_test)
        if (transy == T) {
          pip_temp <- pipeline[[pipi]](y_train_clean)
          y_train <- trans(pip_temp, y_train)
          y_test <- trans(pip_temp, y_test)
        }
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

#' Grid Search and Cross Validation with Noisy (Simulation Only)
#'
#' @author Zhang Jiaqi.
#' @param model your model.
#' @param X,y dataset and label.
#' @param X_noisy dataset with noise.
#' @param y_noisy label (contains label noise)
#' @param folds a positive integer indicating the number of folds (sequential
#'              split, compatible with the old \code{K} argument) or a list of
#'              index vectors, where each element contains the test-set row
#'              indices of one fold.
#' @param metrics this parameter receive a metric function.
#' @param param_list parameter list.
#' @param predict_func this parameter receive a function for predict.
#' @param pipeline preprocessing pipline.
#' @param metrics_params set parameter for each metrics (need a list).
#' @param predict_params set parameters for each predict method (need a list).
#' @param model_settings set parameters for model (need a list).
#' @param transy apply transforms defined in `pipeline` on y, default FALSE.
#' @param shuffle if set \code{shuffle==TRUE}, This function will shuffle
#'                the dataset (only used when \code{folds} is a number).
#' @param seed random seed for \code{shuffle} option (only used when
#'             \code{folds} is a number).
#' @param model_seed random_seed for model.
#' @param threads.num the number of threads used for parallel execution.
#' @return return a metric matrix
#' @import foreach
#' @import doParallel
#' @import doSNOW
#' @import stats
#' @export
grid_search_cv_Xynoisy <- function(model, X, y, X_noisy, y_noisy, folds = 5, metrics, param_list,
                                   predict_func = predict,
                                   pipeline = NULL,
                                   metrics_params = NULL, predict_params = NULL,
                                   model_settings = NULL, transy = FALSE,
                                   shuffle = TRUE, seed = NULL, model_seed = NULL,
                                   threads.num = parallel::detectCores() - 1) {
  s <- Sys.time()
  X <- as.matrix(X)
  y <- as.matrix(y)
  if (is.list(metrics) == F) {
    metrics <- list(metrics)
    names(metrics) <- paste("metric", length(metrics), sep = "")
  }
  n <- nrow(X)
  if (is.numeric(folds) && length(folds) == 1) {
    K <- folds
    if (is.null(seed) == FALSE) {
      set.seed(seed)
    }
    if (shuffle == TRUE) {
      idx <- sample(n)
      X <- X[idx, ]
      X_noisy <- X_noisy[idx, , drop = FALSE]
      y <- y[idx]
      y_noisy <- y_noisy[idx]
    }
  } else {
    K <- length(folds)
  }
  param_grid <- expand.grid(param_list, stringsAsFactors = FALSE)
  n_param <- nrow(param_grid)
  param_names <- colnames(param_grid)
  cl <- parallel::makeCluster(threads.num)
  op <- options(cli.progress_show_after = 0)
  on.exit(options(op), add = TRUE)
  pb <- cli::cli_progress_bar("Grid search", total = n_param, clear = FALSE)
  # pb <- utils::txtProgressBar(max = n_param, style = 3)

  progress <- function(n){cli::cli_progress_update(id = pb, set = n)}
  # progress <- function(n){utils::setTxtProgressBar(pb, n)}
  opts <- list(progress = progress)
  doSNOW::registerDoSNOW(cl)
  i <- 1
  cv_res <- foreach::foreach(i = 1:n_param, .combine = rbind,
                             .packages = c('manysvms', 'Rcpp'),
                             .options.snow = opts) %dopar% {
    params_to_add <- setNames(
      lapply(param_names, function(nm) param_grid[[nm]][[i]]),
      param_names
    )
    params_cv <- list("model" = model,
                     "X" = X, "y" = y, "X_noisy" = X_noisy, "y_noisy" = y_noisy, "folds" = folds,
                     "metrics" = metrics,
                     "predict_func" =  predict_func,
                     "pipeline" = pipeline,
                     "metrics_params" = metrics_params,
                     "model_settings" = append(model_settings, params_to_add),
                     "transy" = transy,
                     "model_seed" = model_seed
    )
    cv_res <- do.call("cross_validation_Xynoisy", params_cv)
    cv_res <- rbind(c(apply(cv_res, 1, mean), apply(cv_res, 1, sd)))
  }
  cli::cli_progress_done(id = pb)
  # close(pb)
  parallel::stopCluster(cl)
  # cat("\n")
  num_metrics <- length(metrics)
  name_matrics <- names(metrics)
  colnames(cv_res)[(num_metrics + 1):(2*num_metrics)] <- paste(name_matrics, "- sd")
  e <- Sys.time()
  # idx_max <- apply(as.matrix(cv_res[,1:num_metrics]), 2, which.max)
  # idx_min <- apply(as.matrix(cv_res[,1:num_metrics]), 2, which.min)

  idx_max <- sapply(seq_len(num_metrics),
                           function(j) select_best_idx(cv_res[, j],
                                                       cv_res[, num_metrics + j], "max"))
  idx_min <- sapply(seq_len(num_metrics),
                           function(j) select_best_idx(cv_res[, j],
                                                       cv_res[, num_metrics + j], "min"))

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
#' \code{cross_validation_Xynoisy} function use noisy data for training,
#' then calculates the average and standard deviation of your metric
#' using clean samples
#'
#' @author Zhang Jiaqi.
#' @param model your model.
#' @param X,y dataset and label.
#' @param X_noisy dataset with noise.
#' @param y_noisy label with label noise.
#' @param folds a positive integer indicating the number of folds (sequential
#'              split, compatible with the old \code{K} argument) or a list of
#'              index vectors, where each element contains the test-set row
#'              indices of one fold.
#' @param metrics this parameter receive a metric function.
#' @param predict_func this parameter receive a function for predict.
#' @param pipeline preprocessing pipline.
#' @param metrics_params set parameters for each metrics (need a list).
#' @param predict_params set parameters for each predict method (need a list).
#' @param model_settings set parameters for model (need a list).
#' @param transy apply transforms defined in `pipeline` on y, default FALSE.
#' @param model_seed random_seed for model.
#' @return return a metric matrix
#' @export
cross_validation_Xynoisy <- function(model, X, y, X_noisy, y_noisy, folds = 5, metrics,
                                     predict_func = predict,
                                     pipeline = NULL,
                                     metrics_params = NULL, predict_params = NULL,
                                     model_settings = NULL, transy = FALSE,
                                     model_seed = NULL) {
  if (is.null(model_seed) == FALSE) {
    set.seed(model_seed)
  }
  X <- as.matrix(X)
  y <- as.matrix(y)
  X_noisy <- as.matrix(X_noisy)
  y_noisy <- as.matrix(y_noisy)
  n <- nrow(X)
  folds <- folds_check_cv(folds, n)
  K <- length(folds)
  metrics <- metrics_check_cv(metrics)
  num_metric <- length(metrics)
  metrics_params <- metrics_params_check_cv(num_metric, metrics_params)
  metric_mat <- matrix(0, num_metric, K)
  for (i in 1:K) {
    idx <- folds[[i]]
    X_test <- X[idx, , drop = FALSE]
    y_test <- y[idx]
    if (K == 1) {
      X_train <- X_test
      y_train <- y_test
    }else{
      X_train_clean <- X[-idx, , drop = FALSE]
      X_train <- X_noisy[-idx, , drop = FALSE]
      y_train <- y_noisy[-idx]
      y_train_clean <- y[-idx]
    }
    if (is.null(pipeline) == F) {
      for (pipi in 1:length(pipeline)) {
        pip_temp <- pipeline[[pipi]](X_train_clean)
        X_train <- trans(pip_temp, X_train)
        X_test <- trans(pip_temp, X_test)
        if (transy == T) {
          pip_temp <- pipeline[[pipi]](y_train_clean)
          y_train <- trans(pip_temp, y_train)
          y_test <- trans(pip_temp, y_test)
        }
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

select_best_idx <- function(scores, sds, direction = c("max", "min")) {
  direction <- match.arg(direction)
  extreme <- if (direction == "max") max(scores) else min(scores)
  ties <- which(scores == extreme)
  if (length(ties) > 1) {
    ties <- ties[which.min(sds[ties])]
  }
  return(ties)
}
