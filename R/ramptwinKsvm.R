#' Ramp loss Twin K Support Vector Machine for Multi-classification
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
#' @param cccp.steps the number of iterations of Concave–Convex Procedure (CCCP).
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param rcpp speed up your code with Rcpp, default \code{rcpp = TRUE}.
#' @return return ramptwinKsvm object
#' @export

ramptwinKsvm <- function(X, y,
                         Ck = rep(1, 4),
                         sk = rep(0.5, 4),
                         kernel = c('linear', 'rbf', 'poly'),
                         gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                         reg = 1, kernel_rect = 1,
                         eps = 0.1,
                         tol = 1e-5, cccp.steps = 30,max.steps = 300,
                         rcpp = TRUE){
  kernel <- match.arg(kernel)

  X <- as.matrix(X)
  y <- as.matrix(y)

  m <- nrow(X)
  n <- ncol(X)

  class_set <- unique(as.matrix(y))
  class_num <- length(class_set)

  if(class_num <=2){
    return(0)
  }

  if(kernel == 'linear'){
    coef_dim <- n
  }else{
    coef_dim <- m * kernel_rect
  }

  coef_list_pos <- matrix(0, nrow = coef_dim, ncol = class_num * (class_num - 1)/2)
  coef_list_neg <- matrix(0, nrow = coef_dim, ncol = class_num * (class_num - 1)/2)

  intercept_list_pos <- matrix(0, ncol = class_num * (class_num - 1)/2)
  intercept_list_neg <- matrix(0, ncol = class_num * (class_num - 1)/2)

  idx <- 0
  for(i in 1:class_num){
    for(j in i:class_num){
      if(i == j){
        next
      }
      idx <- idx + 1

      idxA <- which(y == class_set[i])
      idxB <- which(y == class_set[j])
      idxC <- which(y != class_set[i] & y != class_set[j])

      A <- X[idxA, ]
      B <- X[idxB, ]
      C <- X[idxC, ]

      mA <- nrow(A)
      mB <- nrow(B)
      mC <- m - mA -mB

      delta_pos <- matrix(0, nrow = mB)
      delta_neg <- matrix(0, nrow = mA)

      theta_pos <- matrix(0, nrow = mC)
      theta_neg <- matrix(0, nrow = mC)

      e1 <- matrix(1, nrow = mA)
      e2 <- matrix(1, nrow = mB)
      e3 <- matrix(1, nrow = mC)

      if(kernel == 'linear'){
        X1 <- A
        X2 <- B
        X3 <- C
      }else{
        kernel_m <- round(m*kernel_rect, 0)

        X1 <- kernel_function(A, X[1:kernel_m, ],
                             kernel.type = kernel,
                             gamma = gamma, degree = degree, coef0 = coef0,
                             rcpp = rcpp)
        X2 <- kernel_function(B, X[1:kernel_m, ],
                             kernel.type = kernel,
                             gamma = gamma, degree = degree, coef0 = coef0,
                             rcpp = rcpp)
        X3 <- kernel_function(C, X[1:kernel_m, ],
                             kernel.type = kernel,
                             gamma = gamma, degree = degree, coef0 = coef0,
                             rcpp = rcpp)
      }

      X1 <- cbind(X1, e1)
      X2 <- cbind(X2, e2)
      X3 <- cbind(X3, e3)

      N <- rbind(X2, X3)
      P <- rbind(X1, X3)

      e4 <- rbind(e2, e3 * (1 - eps))
      e5 <- rbind(e1, e3 * (1 - eps))

      for (step in 1:cccp.steps){

        X1TX1_reg_inv <- solve(t(X1) %*% X1 + diag(rep(reg, ncol(X1))))
        H <- N %*% X1TX1_reg_inv %*% t(N)

        lb_pos <- matrix(0, nrow = m - mA)
        ub_pos1 <- matrix(Ck[1], nrow = mB)
        ub_pos2 <- matrix(Ck[2], nrow = mC)
        ub_pos <- rbind(ub_pos1, ub_pos2)

        qp1_solver <- clip_dcd_optimizer(H, e4, lb_pos, ub_pos,
                                         tol, max.steps, rcpp)
        gammas_pos <- as.matrix(qp1_solver$x)
        u_pos <- - X1TX1_reg_inv %*% (
                   t(X2) %*% (gammas_pos[0:mB] - delta_pos) +
                   t(X3) %*% (gammas_pos[(mB+1):length(gammas_pos)] - theta_pos))

        X2TX2_reg_inv <- solve(t(X2) %*% X2 + diag(rep(reg, ncol(X2))))
        H <- P %*% X2TX2_reg_inv %*% t(P)

        lb_neg <- matrix(0, nrow = m - mB)
        ub_neg1 <- matrix(Ck[1], nrow = mA)
        ub_neg2 <- matrix(Ck[2], nrow = mC)
        ub_neg <- rbind(ub_neg1, ub_neg2)

        qp2_solver <- clip_dcd_optimizer(H, e5, lb_neg, ub_neg,
                                         tol, max.steps, rcpp)
        gammas_neg <- as.matrix(qp2_solver$x)
        u_neg <- X2TX2_reg_inv %*% (
                 t(X1) %*% (gammas_neg[0:mA] - delta_neg) +
                 t(X3) %*% (gammas_neg[(mA+1):length(gammas_neg)] - theta_neg))

        idx_delta_pos <- which(- X2 %*% u_pos < sk[1])
        idx_delta_neg <- which(X1 %*% u_neg < sk[3])

        idx_theta_pos <- which(- X3 %*% u_pos < sk[2])
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

        cnt <- 0
        if(sum(delta_pos_new == delta_pos)!=mB){
          delta_pos <- delta_pos_new
          cnt <- cnt + 1
        }
        if(sum(delta_neg_new == delta_neg)!=mA){
          delta_neg <- delta_neg_new
          cnt <- cnt + 1
        }
        if(sum(theta_pos_new == theta_pos)!=mC){
          theta_pos <- theta_pos_new
          cnt <- cnt + 1
        }
        if(sum(theta_neg_new == theta_neg)!=mC){
          theta_neg <- theta_neg_new
          cnt <- cnt + 1
        }
        if(cnt == 0){
          break
        }
      }

      coef_list_pos[, idx] <- u_pos[1:coef_dim, ]
      intercept_list_pos[idx] <- u_pos[coef_dim+1]

      coef_list_neg[, idx] <- u_neg[1:coef_dim, ]
      intercept_list_neg[idx] <- u_neg[coef_dim+1]
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


#' Predict Method for Twin K Support Vector Machine
#'
#' @author Zhang Jiaqi
#' @param object Object of class `twinKsvm`.
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

  if(object$kernel == 'linear'){
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

  for(i in 1:class_num){
    for(j in i:class_num){
      if(i == j){
        next
      }
      idx <- idx + 1
      A <- kernelX %*% object$coef_pos[, idx] + object$intercept_pos[idx]
      B <- kernelX %*% object$coef_neg[, idx] + object$intercept_neg[idx]

      idxA <- which(A > object$eps - 1)
      idxB <- which(B < 1 - object$eps)


      idx_uni <- unique(idxA, idxB)
      vote_mat[idxA, i] <- vote_mat[idxA, i] + 1
      vote_mat[idxB, j] <- vote_mat[idxB, j] + 1
      if(length(idx_uni)!=0){
        vote_mat[-idx_uni, i] <- vote_mat[-idx_uni, i] + 1
        vote_mat[-idx_uni, j] <- vote_mat[-idx_uni, j] + 1
      }
    }
  }

  idx <- apply(vote_mat, 1, which.max)

  pred <- rep(0, m)
  for(i in 1:m){
    pred[i] <- object$class_set[idx[i]]
  }
  acc <- sum(pred == y) / m
  cat('kernel type :', object$kernel, '\n')
  cat(paste("total accuracy :", acc*100, '% \n'))
  predlist <- list("accuracy"= acc,
                   'vote_mat' = vote_mat, 'predict' = pred)
  return(predlist)
}


#' Computes K-fold Cross-Validation Accuracy for ramp Twin K svm
#'
#' @author Zhang Jiaqi
#' @param X a new data frame for predicting.
#' @param y a label data frame corresponding to X.
#' @param K number of folds.
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
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param cccp.steps the number of iterations of Concave–Convex Procedure (CCCP).
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param rcpp speed up your code with Rcpp, default \code{rcpp = TRUE}.
#' @param shuffer if set \code{shuffer==TRUE}, This function will shuffle the dataset.
#' @param seed random seed for \code{shuffer} option.
#' @export

cv.ramptwinKsvm <- function(X, y, K = 5,
                            Ck = rep(1, 4),
                            sk = rep(0.5, 4),
                            kernel = c('linear', 'rbf', 'poly'),
                            gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                            reg = 1, kernel_rect = 1,
                            eps = 0.1,
                            tol = 1e-5, cccp.steps = 30,max.steps = 300,
                            rcpp = TRUE,
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
    ramptwinKsvm_model <- ramptwinKsvm(X, y, Ck = rep(1, 4),
                              sk = rep(0.5, 4),
                              kernel = c('linear', 'rbf', 'poly'),
                              gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                              reg = 1, kernel_rect = 1, eps = 0.1,
                              tol = tol, cccp.steps = cccp.steps,
                              max.steps = max.steps, rcpp = rcpp)
    pred <- predict(ramptwinKsvm_model, test_X, test_y)
    accuracy_list <- append(accuracy_list, pred$accuracy)
  }
  avg_acc <- mean(accuracy_list)
  cat('average accuracy in ',K, 'fold cross validation :', 100*avg_acc, '%\n')
  return(avg_acc)
}
