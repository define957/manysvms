#' Ramp loss Twin K-Class Support Vector Machine for Multi-classification
#'
#' @author Zhang Jiaqi
#' @param X,y dataset and label.
#' @param Ck plenty term vector.
#' @param sk parameter for ramp loss.
#' @param kernel kernel function.
#' @param gamma rbf kernel parameter.
#' @param reg regularization term.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param kernel_rect set kernel size. \code{0<= kernel_rect <= 1}
#' @param eps parameter for rest class.
#' @param tol the precision of the optimization algorithm.
#' @param sig the precision of the CCCP algorithm.
#' @param step_cccp the number of iterations of Concave–Convex Procedure (CCCP).
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param rcpp speed up your code with Rcpp, default \code{rcpp = TRUE}.
#' @return return ramptwinKsvm object
#' @export

ramptwinKsvm <- function(X, y,
                         Ck = rep(1, 4),
                         sk = rep(0.5, 4),
                         kernel = c('linear', 'rbf', 'poly'),
                         gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                         reg = 1e-7, kernel_rect = 1,
                         eps = 0.1,
                         tol = 1e-5, step_cccp = 10, max.steps = 200,
                         rcpp = TRUE, sig = 1e-1){
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
    coef_dim <- m * kernel_rect
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

      X1 <- KernelX[idxA, ]
      dim(X1) <- c(mA, coef_dim)
      X2 <- KernelX[idxB, ]
      dim(X2) <- c(mB, coef_dim)
      X3 <- KernelX[idxC, ]
      dim(X3) <- c(mC, coef_dim)

      delta_pos <- matrix(0, nrow = mB)
      delta_neg <- matrix(0, nrow = mA)

      theta_pos <- matrix(0, nrow = mC)
      theta_neg <- matrix(0, nrow = mC)

      e1 <- matrix(1, nrow = mA)
      e2 <- matrix(1, nrow = mB)
      e3 <- matrix(1, nrow = mC)

      X1 <- cbind(X1, e1)
      X2 <- cbind(X2, e2)
      X3 <- cbind(X3, e3)

      N <- rbind(X2, X3)
      P <- rbind(X1, X3)

      e4 <- rbind(e2, e3 * (1 - eps))
      e5 <- rbind(e1, e3 * (1 - eps))

      inv_mat <- t(X1) %*% X1 + diag(rep(reg, ncol(X1)))
      X1TX1_reg_inv <- chol2inv(chol(inv_mat))
      X1TX1_reg_inv_NT <- X1TX1_reg_inv %*% t(N)
      H1 <- N %*% X1TX1_reg_inv_NT

      inv_mat <- t(X2) %*% X2 + diag(rep(reg, ncol(X2)))
      X2TX2_reg_inv <- chol2inv(chol(inv_mat))
      X2TX2_reg_inv_PT <- X2TX2_reg_inv %*% t(P)
      H2 <- P %*% X2TX2_reg_inv_PT

      lb_pos <- matrix(0, nrow = m - mA)
      ub_pos1 <- matrix(Ck[1], nrow = mB)
      ub_pos2 <- matrix(Ck[2], nrow = mC)
      ub_pos <- rbind(ub_pos1, ub_pos2)

      lb_neg <- matrix(0, nrow = m - mB)
      ub_neg1 <- matrix(Ck[3], nrow = mA)
      ub_neg2 <- matrix(Ck[4], nrow = mC)
      ub_neg <- rbind(ub_neg1, ub_neg2)

      u0_pos <- (lb_pos + ub_pos) / 2
      u0_neg <- (lb_neg + ub_neg) / 2
      for (step in 1:step_cccp) {

        tau_pos <- rbind(delta_pos, theta_pos)
        tau_neg <- rbind(delta_neg, theta_neg)

        q_pos <- t(t(tau_pos) %*% H1) + e4
        q_neg <- t(t(tau_neg) %*% H2) + e5

        qp1_solver <- clip_dcd_optimizer(H1, q_pos, lb_pos, ub_pos,
                                         tol, max.steps, rcpp)
        gammas_pos <- as.matrix(qp1_solver$x)
        u_pos <- -X1TX1_reg_inv_NT %*% (gammas_pos - tau_pos)

        qp2_solver <- clip_dcd_optimizer(H2, q_neg, lb_neg, ub_neg,
                                         tol, max.steps, rcpp)
        gammas_neg <- as.matrix(qp2_solver$x)
        u_neg <- X2TX2_reg_inv_PT %*% (gammas_neg - tau_neg)

        idx_delta_pos <- which(-X2 %*% u_pos < sk[1])
        idx_delta_neg <- which(X1 %*% u_neg < sk[3])

        idx_theta_pos <- which(-X3 %*% u_pos < sk[2])
        idx_theta_neg <- which(X3 %*% u_neg < sk[4])

        delta_pos_new <- matrix(0, nrow = mB)
        delta_neg_new <- matrix(0, nrow = mA)

        theta_pos_new <- matrix(0, nrow = mC)
        theta_neg_new <- matrix(0, nrow = mC)

        delta_pos_new[idx_delta_pos] <- Ck[1]
        delta_pos_new[-idx_delta_pos] <- 0

        delta_neg_new[idx_delta_neg] <- Ck[3]
        delta_neg_new[-idx_delta_neg] <- 0

        theta_pos_new[idx_theta_pos] <- Ck[2]
        theta_pos_new[-idx_theta_pos] <- 0

        theta_neg_new[idx_theta_neg] <- Ck[4]
        theta_neg_new[-idx_theta_neg] <- 0

        if ((norm(u0_pos - gammas_pos, type = "2") < sig) &&
            (norm(u0_neg - gammas_neg, type = "2") < sig)) {
          break
        }else{
          u0_pos <- gammas_pos
          u0_neg <- gammas_neg
        }

        # cnt <- 0
        # if (sum(delta_pos_new == delta_pos) != mB) {
        #   delta_pos <- delta_pos_new
        #   cnt <- cnt + 1
        # }
        # if (sum(delta_neg_new == delta_neg) != mA) {
        #   delta_neg <- delta_neg_new
        #   cnt <- cnt + 1
        # }
        # if (sum(theta_pos_new == theta_pos) != mC) {
        #   theta_pos <- theta_pos_new
        #   cnt <- cnt + 1
        # }
        # if (sum(theta_neg_new == theta_neg) != mC) {
        #   theta_neg <- theta_neg_new
        #   cnt <- cnt + 1
        # }
        # if (cnt == 0) {
        #   break
        # }
      }

      coef_list_pos[, idx] <- u_pos[1:coef_dim, ]
      intercept_list_pos[idx] <- u_pos[coef_dim + 1]

      coef_list_neg[, idx] <- u_neg[1:coef_dim, ]
      intercept_list_neg[idx] <- u_neg[coef_dim + 1]
    }
  }
  ramptwinKsvm <- list('X' = X, 'y' = y,
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
  class(ramptwinKsvm) <- "ramptwinKsvm"
  return(ramptwinKsvm)
}


#' Predict Method for Ramp loss Twin K-Class Support Vector Machine
#'
#' @author Zhang Jiaqi
#' @param object Object of class `ramptwinKsvm`.
#' @param X A new data frame for predicting.
#' @param y A label data frame corresponding to X.
#' @param ... unused parameter.
#' @importFrom stats predict
#' @export
#
predict.ramptwinKsvm <- function(object, X, y, ...){
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


#' Computes K-fold Cross-Validation Accuracy for ramp Twin K svm
#'
#' @author Zhang Jiaqi
#' @param X,y dataset and label.
#' @param K number of folds.
#' @param C1,C2 plenty term.
#' @param s parameter for ramp loss.
#' @param kernel kernel function.
#' @param gamma rbf kernel parameter.
#' @param reg regularization term.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param kernel_rect set kernel size. \code{0<= kernel_rect <= 1}
#' @param eps parameter for rest class.
#' @param tol the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param step_cccp the number of iterations of Concave–Convex Procedure (CCCP).
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param rcpp speed up your code with Rcpp, default \code{rcpp = TRUE}.
#' @param shuffle if set \code{shuffer==TRUE}, This function will shuffle the dataset.
#' @param seed random seed for \code{shuffer} option.
#' @param threads.num The number of threads used for parallel execution.
#' @export

cv.ramptwinKsvm <- function(X, y, K = 5,
                            C1 = 1, C2 = 1, s = 0.5,
                            kernel = c('linear', 'rbf', 'poly'),
                            gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                            reg = 1e-7, kernel_rect = 1,
                            eps = 0.1,
                            tol = 1e-5, max.steps = 200, step_cccp = 30, rcpp = TRUE,
                            shuffle = TRUE, seed = NULL,
                            threads.num = parallel::detectCores() - 1){
  X <- as.matrix(X)
  y <- as.matrix(y)

  param <- expand.grid(C1, C2, s, gamma, degree, coef0, eps)
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
    accuracy_list <- rep(0, K)
    for (i in 1:K) {
      new_idx_k <- new_idx[indx_cv:(indx_cv + v_size - 1)] #get test dataset
      indx_cv <- indx_cv + v_size
      test_X <- X[new_idx_k, ]
      train_X <- X[-new_idx_k, ]
      test_y <- y[new_idx_k]
      train_y <- y[-new_idx_k]
      ramptwinKsvm_model <- ramptwinKsvm(train_X, train_y,
                                         Ck = c(param[j, 1], param[j, 2],
                                                param[j, 1], param[j, 2]),
                                         sk = param[j, 3] * rep(1, 4),
                                         kernel = kernel,
                                         gamma = param[j, 4],
                                         degree = param[j, 5],
                                         coef0 = param[j, 6],
                                         reg = reg, kernel_rect = kernel_rect,
                                         eps = param[j, 7],
                                         tol = tol, max.steps = max.steps,
                                         step_cccp = step_cccp,rcpp = rcpp)
      pred <- predict(ramptwinKsvm_model, test_X, test_y)
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
      "C1 = C3 = ", param[max_idx, 1],
      "C2 = C4 = ", param[max_idx, 2],
      "s = ", param[max_idx, 3],
      "eps =", param[max_idx, 7],
      "\n",
      "gamma = ", param[max_idx, 4],
      "degree = ",param[max_idx, 5],
      "coef0 =", param[max_idx, 6],
      "\n")
  cat("Accuracy :", as.numeric(res[max_idx, 1]),
      "Sd :", as.numeric(res[max_idx, 2]), "\n")
  return(res)
}
