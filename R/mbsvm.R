#' Multiple Birth Support Vector Machine
#'
#' \code{mbsvm} is an R implementation of multiple birth SVM
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param Ck plenty term list.
#' @param kernel kernel function.
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param reg regularization term to take care of problems due to ill-conditioning in dual problem.
#' @param kernel_rect set kernel size. \code{0<= kernel_rect <= 1}
#' @param tol the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param rcpp speed up your code with Rcpp, default \code{rcpp = TRUE}.
#' @return return mbsvm object.
#' @export
#' @examples
#' library(manysvms)
#' data('iris')
#'
#'
#' X <- iris[1:150, 1:4]
#' y <- iris[1:150, 5]
#'
#' model1 <- mbsvm(X, y, kernel = 'linear')
#' pred <-predict(model1, X, y)

mbsvm <- function(X, y,
                  Ck = rep(1, length(unique(y))),
                  kernel = c('linear', 'rbf', 'poly'),
                  gamma = 1 / ncol(X), reg = 1, kernel_rect = 1,
                  tol = 1e-6, max.steps = 300,
                  rcpp = TRUE){
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
    coef_dim <- round(kernel_rect*m)
  }
  coef_list <- rep(0, coef_dim*class_num)
  dim(coef_list) <- c(coef_dim, class_num)
  intercept_list <- rep(0, class_num)
  coef_list <- as.matrix(coef_list)
  intercept_list <- as.matrix(intercept_list)


  for(k in 1:class_num){
    idx_Ak <- which(y == class_set[k])
    idx_Bk <- which(y != class_set[k])

    A <- X[idx_Ak, ]
    B <- X[idx_Bk, ]

    mA <- nrow(A)
    mB <- nrow(B)

    e1 <- as.matrix(rep(1, mA))
    dim(e1) <- c(mA, 1)
    e2 <- as.matrix(rep(1, mB))
    dim(e2) <- c(mB, 1)

    if(kernel == 'linear'){
      S <- cbind(A, e1)
      R <- cbind(B, e2)
    }else{
      kernel_m <- round(m*kernel_rect, 0)
      if(rcpp == TRUE){
        S = cpp_rbf_kernel(A, X[1:kernel_m, ], gamma = gamma)
        R = cpp_rbf_kernel(B, X[1:kernel_m, ], gamma = gamma)
      }else{
        S <- matrix(0, nrow = mA, ncol = kernel_m)
        R <- matrix(0, nrow = mB, ncol = kernel_m)
        for(i in 1:mA){
          for(j in 1:kernel_m){
              S[i, j] <-  rbf_kernel(A[i, ], X[j, ], gamma = gamma)
          }
        }

        for(i in 1:mB){
          for(j in 1:kernel_m){
              R[i, j] <-  rbf_kernel(B[i, ], X[j, ], gamma = gamma)
          }
        }
      }
      S <- cbind(S, e1)
      R <- cbind(R, e2)
    }
    RTR_reg_inv <- solve(t(R) %*% R + diag(rep(reg, ncol(R))))
    H <- S %*% RTR_reg_inv %*% t(S)
    lbB <- matrix(0, nrow = mA)
    ubB <- matrix(Ck[k], nrow = mA)

    x <- clip_dcd_optimizer(H, e1, lbB, ubB, tol, max.steps, rcpp = rcpp)$x

    gammas <- as.matrix(x)
    Z2 <- RTR_reg_inv %*% t(S) %*% gammas

    coef_list[, k] <- Z2[1:coef_dim]
    intercept_list[k] <- Z2[coef_dim+1]
  }

  mbsvm <- list('X' = X, 'y' = y,
                      'class_set' = class_set,'class_num' = class_num,
                      'coef' = coef_list, 'intercept' = intercept_list,
                      'kernel' = kernel,
                      'gamma' = gamma,
                      'kernel_rect' = kernel_rect,
                      'Rcpp' = rcpp
                )
  class(mbsvm)<-"mbsvm"
  return(mbsvm)
}


#' Predict Method for Multiple Birth Support Vector Machines
#'
#' @author Zhang Jiaqi
#' @param object A fitted object of class inheriting from \code{mbsvm}.
#' @param X A new data frame for predicting.
#' @param y A label data frame corresponding to X.
#' @param show.info If you set \code{show.info = TRUE}, this function
#'  will print the prediction results and other information.
#' @param ... unused parameter.
#' @importFrom stats predict
#' @export

predict.mbsvm <- function(object, X, y = NULL, show.info = TRUE, ...){
  X <- as.matrix(X)
  m <- nrow(X)
  dis_mat <- matrix(0, nrow = m, ncol = object$class_num)
  X <- as.matrix(X)
  m1 <- nrow(X)
  m2 <- round(nrow(object$X)*object$kernel_rect, 0)
  # get Kernel X
  if(object$kernel == 'linear'){
    kernelX <- X
  }else{
    if(object$Rcpp == TRUE){
      kernelX <- cpp_rbf_kernel(X, object$X[1:m2,], gamma = object$gamma)
    }else{
      if(object$kernel != 'linear'){
        kernelX <- matrix(0, nrow = m1, ncol = m2)
        for(i in 1:m1){
          for(j in 1:m2){
            if(object$kernel == 'rbf'){
              kernelX[i, j] <- rbf_kernel(X[i, ], object$X[j, ], gamma = object$gamma)
            }
          }
        }
      }
    }
  }


  class_num <- object$class_num
  for(i in 1:class_num){
    dis_mat[, i] <- abs(kernelX %*% object$coef[, i] + object$intercept[i]) /
      norm(as.matrix(object$coef[, i]), type = "2")
  }
  pred <- rep(0, m)
  for(i in 1:m){
    idx <- which.max(dis_mat[i, ])
    pred[i] <- object$class_set[idx]
  }
  if(is.null(y) == FALSE){
    y <- as.matrix(y)
    acc <- sum(pred == y)/length(y)
    if(show.info == TRUE){
      cat('kernel type :', object$kernel, '\n')
      cat(paste("total accuracy :", acc*100, '% \n'))
    }
    predlist <- list("accuracy"= acc,
                     'distance' = dis_mat,
                     'predict' = pred)
  }else{
    predlist <- list('distance' = dis_mat,
                     'predict' = pred)
  }
  return(predlist)
}


#' Computes K-fold Cross-Validated Accuracy for MBSVM
#'
#' @author Zhang Jiaqi
#' @param mbsvm a fitted object of class inheriting from \code{mbsvm}.
#' @param X A new data frame for predicting.
#' @param y A label data frame corresponding to X.
#' @param K Number of folds.
#' @param kernel kernel type. Default value \code{kernel = 'linear'}.
#' @param kernel_rect set kernel size. \code{0<= kernel_rect <= 1}
#' @param gamma parameter needed for rbf kernel.
#' @param reg regularization term to take care of problems due to ill-conditioning in dual problem.
#' @param shuffer if set \code{shuffer==TRUE}, This function will shuffle the dataset.
#' @param seed random seed for \code{shuffer} option.
#' @export

cv.mbsvm <- function(X, y , K = 5,
                           kernel = c('linear', 'rbf', 'poly'),
                           reg = 1, gamma = 1/ncol(X), kernel_rect = 1,
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
    mbsvm_model <- mbsvm(train_X, train_y, kernel = kernel, reg = reg,
                         gamma = gamma, kernel_rect = kernel_rect)
    pred <- predict(mbsvm_model, test_X, test_y)
    accuracy_list <- append(accuracy_list, pred$accuracy)
  }
  avg_acc <- mean(accuracy_list)
  cat('average accuracy in ',K, 'fold cross validation :', 100*avg_acc, '%\n')
  return(avg_acc)
}


#' Plot Method for Multiple Birth Support Vector Machines
#'
#' If there are only two variables in the training data set,
#' then you can call this function to draw the decision boundary
#' for your \code{mbsvm}.
#'
#' @author Zhang Jiaqi
#' @param x Object of class \code{mbsvm}.
#' @param xlab,ylab your x-axis name and y-axis name.
#' @param ... unused parameter.
#' @importFrom graphics plot
#' @import ggplot2
#' @export
#' @examples
#'
#' library(manysvms)
#' library(MASS)
#' set.seed(112)
#' n <- 300
#' n_c <- 3
#' sig <- diag(c(0.02, 0.03))
#' x1 <- mvrnorm(n/n_c, mu = c(-0.6, 0), Sigma = sig)
#' x2 <- mvrnorm(n/n_c, mu = c(0.6, 0), Sigma = sig)
#' x3 <- mvrnorm(n/n_c, mu = c(0, 0.6), Sigma = sig)
#' X <- rbind(x1, x2, x3)
#' y <- rep(c(1,2,3), rep(n/n_c, n_c))
#' mbsvm_model = mbsvm(X, y, reg = 1e-3)
#' plot(mbsvm_model)

plot.mbsvm <- function(x, xlab = 'x1', ylab = 'x2', ...){
  x1 <- NULL
  x2 <- NULL
  Class <- NULL
  X <- x$X
  y <- x$y
  grid <- expand.grid(seq(min(X[, 1]), max(X[, 1]),length.out=20),
                      seq(min(X[, 2]), max(X[, 2]),length.out=20))
  Z <- predict(x, grid, show.info = FALSE)$predict
  contour_data <- data.frame(grid, Z)
  colnames(contour_data) <- c("x1", "x2", "Class")

  contour_data$Class <- as.factor(contour_data$Class)
  Xydata = cbind(X, y)
  Xydata <- as.data.frame(Xydata)
  colnames(Xydata) <- c("x1", "x2", "Class")
  Xydata$Class <- as.factor(Xydata$Class)

  p <- ggplot2::ggplot(contour_data, aes(x = x1, y = x2)) +
       ggplot2::geom_tile(aes(fill = Class,
                                     alpha = 0.2), show.legend = FALSE) +
       ggplot2::geom_point(data = Xydata,
                           aes(x = x1, y= x2,
                             color = Class, shape = Class), size = 3) +
       ggplot2::labs(x=xlab,y=ylab) +
       ggplot2::theme_bw()
  return(p)
}
