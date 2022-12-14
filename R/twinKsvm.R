#' Twin K Support Vector Machine for Multi-classification
#'
#' @author Zhang Jiaqi
#' @param X,y dataset and label.
#' @param Ck plenty term vector.
#' @param kernel kernel function.
#' @param gamma rbf kernel parameter.
#' @param reg regularization term.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param kernel_rect set kernel size. \code{0<= kernel_rect <= 1}
#' @param eps parameter for rest class.
#' @param tol the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param rcpp speed up your code with Rcpp, default \code{rcpp = TRUE}.
#' @return return twinKsvm object
#' @export
twinKsvm <- function(X, y,
                     Ck = rep(1, 4),
                     kernel = c('linear', 'rbf', 'poly'),
                     gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                     reg = 1e-7, kernel_rect = 1,
                     eps = 0.1,
                     tol = 1e-5, max.steps = 200, rcpp = TRUE){
  kernel <- match.arg(kernel)

  X <- as.matrix(X)
  y <- as.matrix(y)

  m <- nrow(X)
  n <- ncol(X)

  class_set <- unique(as.matrix(y))
  class_num <- length(class_set)

  if (class_num <= 2) {
    return(0)
  }

  if (kernel == 'linear') {
    coef_dim <- n
  }else{
    coef_dim <- round(m * kernel_rect, 0)
  }

  coef_list_pos <- matrix(0, nrow = coef_dim, ncol = class_num * (class_num - 1)/2)
  coef_list_neg <- matrix(0, nrow = coef_dim, ncol = class_num * (class_num - 1)/2)

  intercept_list_pos <- matrix(0, ncol = class_num * (class_num - 1)/2)
  intercept_list_neg <- matrix(0, ncol = class_num * (class_num - 1)/2)

  if (kernel != 'linear') {
    kernel_m <- round(m*kernel_rect, 0)
    KernelX <- kernel_function(X, X[1:kernel_m, ],
                               kernel.type = kernel,
                               gamma = gamma, degree = degree, coef0 = coef0,
                               rcpp = rcpp)
  }else{
    KernelX <- X
  }

  # solve K * K-1 models
  idx <- 0
  for (i in 1:class_num) {
    for (j in i:class_num) {
      if (i == j) {
        next
      }
      idx <- idx + 1

      idxA <- which(y == class_set[i])
      idxB <- which(y == class_set[j])
      idxC <- which(y != class_set[i] & y != class_set[j])

      mA <- length(idxA)
      mB <- length(idxB)
      mC <- m - mA - mB

      S <- KernelX[idxA, ]
      dim(S) <- c(mA, coef_dim)
      R <- KernelX[idxB, ]
      dim(R) <- c(mB, coef_dim)
      W <- KernelX[idxC, ]
      dim(W) <- c(mC, coef_dim)

      e1 <- matrix(1, nrow = mA)
      e2 <- matrix(1, nrow = mB)
      e3 <- matrix(1, nrow = mC)

      S <- cbind(S, e1)
      R <- cbind(R, e2)
      W <- cbind(W, e3)

      # slove plane 1
      STS_reg_inv <- solve(t(S) %*% S + diag(rep(reg, ncol(S))))
      V <- rbind(R, W)
      H <- V %*% STS_reg_inv %*% t(V)
      e4 <- rbind(e2, e3 * (1 - eps))
      lbA <- matrix(0, nrow = m - mA)
      ubA1 <- matrix(Ck[1], nrow = mB)
      ubA2 <- matrix(Ck[2], nrow = mC)
      ubA <- rbind(ubA1, ubA2)
      qp1_solver <- clip_dcd_optimizer(H, e4, lbA, ubA, tol, max.steps, rcpp)
      gammas <- as.matrix(qp1_solver$x)
      Z1 <- -STS_reg_inv %*% (t(R) %*% gammas[0:mB] +
                              t(W) %*% gammas[(mB + 1):length(gammas)])
      coef_list_pos[, idx] <- Z1[1:coef_dim, ]
      intercept_list_pos[idx] <- Z1[coef_dim + 1]

      RTR_reg_inv <- solve(t(R) %*% R + diag(rep(reg, ncol(R))))
      J <- rbind(S, W)
      H <- J %*% RTR_reg_inv %*% t(J)
      e5 <- rbind(e1, e3 * (1 - eps))
      lbB <- matrix(0, nrow = m - mB)
      ubB1 <- matrix(Ck[3], nrow = mA)
      ubB2 <- matrix(Ck[4], nrow = mC)
      ubB <- rbind(ubB1, ubB2)
      qp2_solver <- clip_dcd_optimizer(H, e5, lbB, ubB, tol, max.steps, rcpp)
      rhos <- as.matrix(qp2_solver$x)
      Z2 <- RTR_reg_inv %*% t(J) %*% rhos
      coef_list_neg[, idx] <- Z2[1:coef_dim, ]
      intercept_list_neg[idx] <- Z2[coef_dim + 1]
    }
  }
  twinKsvm <- list('X' = X, 'y' = y,
                   'class_set' = class_set,'class_num' = class_num,
                   'coef_pos' = coef_list_pos, 'coef_neg' = coef_list_neg,
                   'intercept_pos' = intercept_list_pos,
                   'intercept_neg' = intercept_list_neg,
                   'eps' = eps,
                   'kernel' = kernel,
                   'gamma' = gamma,
                   'coef0' = coef0,
                   'kernel_rect' = kernel_rect,
                   'Rcpp' = rcpp,
                   'call' = match.call())
  class(twinKsvm) <- "twinKsvm"
  return(twinKsvm)
}

#' Predict Method for Twin K Support Vector Machine
#'
#' @author Zhang Jiaqi
#' @param object Object of class `twinKsvm`.
#' @param X A new data frame for predicting.
#' @param y A label data frame corresponding to X.
#' @param ... unused parameter.
#' @importFrom stats predict
#' @export

predict.twinKsvm <- function(object, X, y, ...){
  X <- as.matrix(X)
  y <- as.matrix(y)
  m <- nrow(X)
  km <- nrow(object$X)
  vote_mat <- matrix(0, nrow = nrow(X), ncol = object$class_num)
  idx <- 0
  class_num <- object$class_num

  if (object$kernel == 'linear') {
    kernelX <- X
  }else{
    kernel_m <- round(km*object$kernel_rect, 0)
    kernelX <- kernel_function(X, object$X[1:kernel_m, ],
                         kernel.type = object$kernel,
                         gamma = object$gamma,
                         degree = object$degree,
                         coef0 = object$coef0,
                         rcpp = object$Rcpp)
  }

  for (i in 1:class_num) {
    for (j in i:class_num) {
      if (i == j) {
        next
      }
      idx <- idx + 1
      A <- kernelX %*% object$coef_pos[, idx] + object$intercept_pos[idx]
      B <- kernelX %*% object$coef_neg[, idx] + object$intercept_neg[idx]

      idxA <- as.matrix(which(A > object$eps - 1))
      idxB <- as.matrix(which(B < 1 - object$eps))

      idx_uni <- unique(rbind(idxA, idxB))
      if (length(idxA) != 0) {
      vote_mat[idxA, i] <- vote_mat[idxA, i] + 1
      }
      if (length(idxB) != 0) {
        vote_mat[idxB, j] <- vote_mat[idxB, j] + 1
      }
      if (length(idx_uni) != 0) {
        vote_mat[-idx_uni, i] <- vote_mat[-idx_uni, i] - 1
        vote_mat[-idx_uni, j] <- vote_mat[-idx_uni, j] - 1
      }
    }
  }

  idx <- apply(vote_mat, 1, which.max)

  pred <- rep(0, m)
  for (i in 1:m) {
    pred[i] <- object$class_set[idx[i]]
  }
  acc <- sum(pred == y) / m
  cat('kernel type :', object$kernel, '\n')
  cat(paste("total accuracy :", acc*100, '% \n'))
  predlist <- list("accuracy" = acc,
                   'vote_mat' = vote_mat, 'predict' = pred)
  return(predlist)
}

#' @rdname twinKsvm
#' @format NULL
#' @usage NULL
#' @export

print.twinKsvm <- function(x, ...){
  cat("\nCall:", deparse(x$call, 0.8 * getOption("width")), "\n", sep = "\n")
  cat("SVM type : ", class(x), "\n")
  cat("SVM kernel : ", x$kernel, "\n")
  if (x$kernel == 'rbf') {
    cat("gamma : ", x$gamma, "\n")
  }
  cat("number of observations : ", nrow(x$X), "\n")
  cat("number of class : ", length(x$class_set), "\n")
}


#' Computes K-fold Cross-Validation Accuracy for Twin K svm
#'
#' @author Zhang Jiaqi
#' @param X a new data frame for predicting.
#' @param y a label data frame corresponding to X.
#' @param K number of folds.
#' @param C1,C2 plenty term.
#' @param kernel kernel function.
#' @param gamma rbf kernel parameter.
#' @param reg regularization term.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param kernel_rect set kernel size. \code{0<= kernel_rect <= 1}
#' @param eps parameter for rest class.
#' @param tol the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param rcpp speed up your code with Rcpp, default \code{rcpp = TRUE}.
#' @param shuffle if set \code{shuffer==TRUE}, This function will shuffle the dataset.
#' @param seed random seed for \code{shuffle} option.
#' @param threads.num The number of threads used for parallel execution.
#' @export

cv.twinKsvm <- function(X, y, K = 5,
                        C1 = 1, C2 = 1,
                        kernel = c('linear', 'rbf', 'poly'),
                        gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                        reg = 1e-7, kernel_rect = 1,
                        eps = 0.1,
                        tol = 1e-5, max.steps = 200, rcpp = TRUE,
                        shuffle = TRUE, seed = NULL,
                        threads.num = parallel::detectCores() - 1){
  X <- as.matrix(X)
  y <- as.matrix(y)

  param <- expand.grid(C1, C2, gamma, degree, coef0, eps)
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
  j <- 1

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
      twinKsvm_model <- twinKsvm(train_X, train_y,
                                 Ck = c(param[j, 1], param[j, 2],
                                        param[j, 1], param[j, 2]),
                                 kernel = kernel,
                                 gamma = param[j, 3],
                                 degree = param[j, 4], coef0 = param[j, 5],
                                 reg = reg, kernel_rect = kernel_rect,
                                 eps = param[j, 6],
                                 tol = tol, max.steps = max.steps, rcpp = TRUE)
      pred <- predict(twinKsvm_model, test_X, test_y)
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
      "C1 = C3 = ", param[max_idx, 1],
      "C2 = C4 = ", param[max_idx, 2],
      "eps =", param[max_idx, 6],
      "\n",
      "gamma = ", param[max_idx, 3],
      "degree = ",param[max_idx, 4],
      "coef0 =", param[max_idx, 5],
      "\n")
  cat("Accuracy :", as.numeric(res[max_idx, 1]),
      "Sd :", as.numeric(res[max_idx, 2]), "\n")
  return(res)
}
