#' Twin K SVM for Multi-classification
#'
#' @author Zhang Jiaqi
#' @param X,y dataset and label
#' @param Ck plenty term vector
#' @param kernel kernel function
#' @param gamma rbf kernel parameter
#' @param reg regularization tern
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
                     C = rep(1, 4),
                     kernel = c('linear', 'rbf', 'poly'),
                     gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                     reg = 1, kernel_rect = 1,
                     eps = 0.1,
                     tol = 1e-5, max.steps = 300, rcpp = TRUE){
  kernel <- match.arg(kernel)

  X <- as.matrix(X)
  y <- as.matrix(y)

  m <- nrow(X)
  n <- ncol(X)

  class_set <- unique(as.matrix(y))
  class_num <- length(class_set)

  if(kernel == 'linear'){
    coef_dim <- n
  }else{
    coef_dim <- m * kernel_rect
  }

  coef_list <- matrix(0, nrow = coef_dim, ncol = class_num * (class_num - 1))
  intercept_list <- matrix(0, ncol = class_num * (class_num - 1))

  # solve K * K-1 models
  for(i in 1:class_num){
    for(j in 1:class_num){
      if(i == j){
        next
      }
      idxA <- which(y == class_set[i])
      idxB <- which(y == class_set[j])
      idxC <- which(y != class_set[i] & y != class_set[j])

      A <- X[idxA, ]
      B <- X[idxB, ]
      C <- X[idxC, ]

      mA <- nrow(A)
      mB <- nrow(B)
      mC <- m - mA -mB
      print(mB)
      e1 <- matrix(1, nrow = mA)
      e2 <- matrix(1, nrow = mB)
      e3 <- matrix(1, nrow = mC)

      if(kernel == 'linear'){
        S <- A
        R <- B
        W <- C
      }else{
        kernel_m <- round(m*kernel_rect, 0)

        H <- kernel_function(A, X[1:kernel_m, ],
                             kernel.type = kernel,
                             gamma = gamma, degree = degree, coef0 = coef0,
                             rcpp = rcpp)
        R <- kernel_function(B, X[1:kernel_m, ],
                             kernel.type = kernel,
                             gamma = gamma, degree = degree, coef0 = coef0,
                             rcpp = rcpp)
        W <- kernel_function(C, X[1:kernel_m, ],
                             kernel.type = kernel,
                             gamma = gamma, degree = degree, coef0 = coef0,
                             rcpp = rcpp)
      }
      S <- cbind(S, e1)
      R <- cbind(R, e2)
      W <- cbind(W, e3)

      # slove plane 1
      STS_reg_inv <- solve(t(S) %*% S + diag(rep(reg, ncol(S))))
      V <- rbind(R, W)
      H <- V %*% STS_reg_inv %*% t(V)
      e4 <- rbind(e2, e3 * (1 - eps))

      lbA <- matrix(0, nrow = m - mA)
      ubB1 <- matrix(C[1], nrow = mB)
      ubB2 <- matrix(C[2], nrow = mC)
      ubA <- rbind(ubB1, ubB2)
      qp1_solver <- clip_dcd_optimizer(H, e4, lbA, ubA, tol, max.steps, rcpp)
      gammas <- as.matrix(qp1_solver$x)
      Z1 <- - STS_reg_inv %*% (t(R) %*% gammas[0:mB] +
                               t(W) %*% gammas[(mB+1):length(gammas)])
      print(Z1)
    }
  }
}
