#' Twin SVM for Multi-classification by Using Ones versus Rest Strategy
#'
#' @author Zhang Jiaqi
#' @param X,y dataset and label
#' @param Ck plenty term vector
#' @param kernel kernel function
#' @param gamma rbf kernel parameter
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param  reg regularization tern
#' @param kernel_rect set kernel size. \code{0<= kernel_rect <= 1}
#' @param tol the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param rcpp speed up your code with Rcpp, default \code{rcpp = TRUE}.
#' @return return twinsvm object
#' @export
#' @examples
#'
#' library(manysvms)
#' data('iris')
#'
#' m <- 100
#' X <- iris[0:m, 1:2]
#' y <- iris[0:m, 5]
#' model <- twinsvm_ovr(X, y, kernel = 'linear')
#'
#' coef(model)
#' pred <-predict(model, X, y)

twinsvm_ovr<- function(X, y,
                       Ck = rep(1, length(unique(y))),
                       kernel = c('linear', 'rbf', 'poly'),
                       gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                       reg = 1, kernel_rect = 1,
                       tol = 1e-5, max.steps = 300, rcpp = TRUE){
  kernel <- match.arg(kernel)

  m <- nrow(X)
  n <- ncol(X)

  X <- as.matrix(X)
  y <- as.matrix(y)

  class_set <- unique(as.matrix(y))
  class_num <- length(class_set)

  if(kernel == 'linear'){
    coef_dim <- n
  }else if(kernel == 'rbf'){
    coef_dim <- m*kernel_rect
  }

  coef_list <- matrix(0, nrow = coef_dim, ncol = class_num)
  intercept_list <- matrix(0, ncol = class_num)


  for(k in 1:class_num){
    idx_Ak <- which(y == class_set[k])
    idx_Bk <- which(y != class_set[k])

    A <- X[idx_Ak, ]
    B <- X[idx_Bk, ]

    mA <- nrow(A)
    mB <- nrow(B)

    e1 <- matrix(1, nrow = mA)
    e2 <- matrix(1, nrow = mB)

    if(kernel == 'linear'){
      S <- A
      R <- B
    }else{
      kernel_m <- round(m*kernel_rect, 0)

      S <- kernel_function(A, X[1:kernel_m, ],
                           kernel.type = kernel,
                           gamma = gamma, degree = degree, coef0 = coef0,
                           rcpp = rcpp)
      R <- kernel_function(B, X[1:kernel_m, ],
                           kernel.type = kernel,
                           gamma = gamma, degree = degree, coef0 = coef0,
                           rcpp = rcpp)
    }
    S <- cbind(S, e1)
    R <- cbind(R, e2)

    STS_reg_inv <- solve(t(S) %*% S + diag(reg, ncol(S)))
    H <- R %*% STS_reg_inv %*% t(R)
    lbA <- matrix(0, nrow = mB)
    ubA <- matrix(Ck[k], nrow = mB)
    qp1_solver <- clip_dcd_optimizer(H, e2, lbA, ubA, tol, max.steps, rcpp)
    alphas <- as.matrix(qp1_solver$x)
    Z1 <- -STS_reg_inv %*% t(R) %*% alphas

    coef_list[, k] <- Z1[1:coef_dim]
    intercept_list[k] <- Z1[coef_dim+1]
  }

  twinsvm_ovr <- list('X' = X, 'y' = y,
                      'class_set' = class_set,'class_num' = class_num,
                      'coef' = coef_list, 'intercept' = intercept_list,
                      'kernel' = kernel,
                      'gamma' = gamma,
                      'Rcpp' = rcpp,
                      'kernel_rect' = kernel_rect
                      )


  class(twinsvm_ovr)<-"twinsvm_ovr"
  return(twinsvm_ovr)
}

#' Predict Method for Twin Support Vector Machines (O v R)
#'
#' @author Zhang Jiaqi
#' @param object Object of class `twinsvm_ovr`.
#' @param X A new data frame for predicting.
#' @param y A label data frame corresponding to X.
#' @param ... unused parameter.
#' @importFrom stats predict
#' @export

predict.twinsvm_ovr <- function(object, X, y, ...){
  X <- as.matrix(X)
  m <- nrow(X)
  dis_mat <- matrix(0, nrow = m, ncol = object$class_num)
  X <- as.matrix(X)
  km <- nrow(object$X)
  if(object$kernel == 'linear'){
    kernelX <- X
  }else{
    kernel_m <- round(km*object$kernel_rect, 0)
    print(km)
    kernelX <- kernel_function(X, object$X[1:kernel_m, ],
                               kernel.type = object$kernel,
                               gamma = object$gamma,
                               degree = object$degree,
                               coef0 = object$coef0,
                               rcpp = object$Rcpp)
  }

  class_num <- object$class_num
  for(i in 1:class_num){
    dis_mat[, i] <- abs(kernelX %*% object$coef[, i] + object$intercept[i]) /
      norm(as.matrix(object$coef[, i]), type = "2")
  }
  pred <- rep(0, m)
  for(i in 1:m){
    idx <- which.min(dis_mat[i, ])
    pred[i] <- object$class_set[idx]
  }
  y <- as.matrix(y)
  acc <- sum(pred == y)/length(y)
  cat('kernel type :', object$kernel, '\n')
  cat(paste("total accuracy :", acc*100, '% \n'))
  predlist <- list("accuracy"= acc,
                   'distance' = dis_mat, 'predict' = pred)
  return(predlist)
}
