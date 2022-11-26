#' Multiple Birth Support Vector Machine
#'
#' \code{mbsvm} is an R implementation of multiple birth SVM
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param Ck plenty term list.
#' @param kernel kernel function. The definitions of various kernel functions are as follows:
#' \describe{
#'     \item{linear:}{\eqn{u'v}{u'*v}}
#'     \item{poly:}{\eqn{(\gamma u'v + coef0)^{degree}}{(gamma*u'*v + coef0)^degree}}
#'     \item{rbf:}{\eqn{e^{(-\gamma |u-v|^2)}}{exp(-gamma*|u-v|^2)}}
#' }
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
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
                  gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                  reg = 1, kernel_rect = 1,
                  tol = 1e-6, max.steps = 300, rcpp = TRUE){

  kernel <- match.arg(kernel)

  m <- nrow(X)
  n <- ncol(X)

  X <- as.matrix(X)
  y <- as.matrix(y)

  class_set <- unique(as.matrix(y))
  class_num <- length(class_set)
  kernel_m <- round(m*kernel_rect, 0)

  if (kernel == 'linear') {
    coef_dim <- n
  }else if (kernel == 'rbf') {
    coef_dim <- round(kernel_rect*m)
  }
  coef_list <- rep(0, coef_dim*class_num)
  dim(coef_list) <- c(coef_dim, class_num)
  intercept_list <- rep(0, class_num)
  coef_list <- as.matrix(coef_list)
  intercept_list <- as.matrix(intercept_list)

  if (kernel != 'linear') {
    kernel_m <- round(m*kernel_rect, 0)
    # x_norm <- as.matrix(apply(X, 1, norm, type = '2')^2)
    # KernelX <- exp(gamma*(- x_norm %*% t(e) - e %*% t(x_norm) + 2* X%*%t(X)))
    KernelX <- kernel_function(X, X[1:kernel_m, ],
                               kernel.type = kernel,
                               gamma = gamma, degree = degree, coef0 = coef0,
                               rcpp = rcpp)
    }else{
      KernelX <- X
    }

  for (k in 1:class_num) {
    idx_Ak <- which(y == class_set[k])
    idx_Bk <- which(y != class_set[k])

    mA <- length(idx_Ak)
    mB <- length(idx_Bk)

    S <- as.matrix(KernelX[idx_Ak, ])
    dim(S)  <- c(mA, coef_dim)
    R <- as.matrix(KernelX[idx_Bk, ])
    dim(R)  <- c(mB, coef_dim)

    e1 <- matrix(1, nrow = mA)
    e2 <- matrix(1, nrow = mB)

    S <- cbind(S, e1)
    R <- cbind(R, e2)

    inv_mat <- t(R) %*% R + diag(rep(reg, ncol(R)))
    RTR_reg_inv_S <- chol2inv(chol(inv_mat)) %*% t(S)
    H <- S %*% RTR_reg_inv_S
    lbB <- matrix(0, nrow = mA)
    ubB <- matrix(Ck[k], nrow = mA)

    x <- clip_dcd_optimizer(H, e1, lbB, ubB, tol, max.steps, rcpp = rcpp)$x

    gammas <- as.matrix(x)
    Z2 <- RTR_reg_inv_S %*% gammas

    coef_list[, k] <- Z2[1:coef_dim]
    intercept_list[k] <- Z2[coef_dim + 1]
  }

  mbsvm <- list('X' = X, 'y' = y,
                      'class_set' = class_set,'class_num' = class_num,
                      'coef' = coef_list, 'intercept' = intercept_list,
                      'kernel' = kernel,
                      'gamma' = gamma,
                      'kernel_rect' = kernel_rect,
                      'Rcpp' = rcpp
                )
  class(mbsvm) <- "mbsvm"
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
  km <- nrow(object$X)
  dis_mat <- matrix(0, nrow = m, ncol = object$class_num)
  X <- as.matrix(X)
  # get Kernel X
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


#' Computes K-fold Cross-Validation Accuracy for MBSVM
#'
#' @author Zhang Jiaqi
#' @param X A new data frame for predicting.
#' @param y A label data frame corresponding to X.
#' @param C plenty term.
#' @param K Number of folds.
#' @param kernel kernel type. Default value \code{kernel = 'linear'}.
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param reg regularization term to take care of problems due to ill-conditioning in dual problem.
#' @param kernel_rect set kernel size. \code{0<= kernel_rect <= 1}
#' @param tol the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param rcpp speed up your code with Rcpp, default \code{rcpp = TRUE}.
#' @param shuffle if set \code{shuffle==TRUE}, This function will shuffle the dataset.
#' @param seed random seed for \code{shuffer} option.
#' @param threads.num The number of threads used for parallel execution.
#' @import foreach
#' @import doParallel
#' @import doSNOW
#' @import stats
#' @export

cv.mbsvm <- function(X, y , K = 5, C = 1,
                     kernel = c('linear', 'rbf', 'poly'),
                     gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                     reg = 1e-7, kernel_rect = 1, tol = 1e-5,
                     max.steps = 200, rcpp = TRUE, shuffle = TRUE, seed = NULL,
                     threads.num = parallel::detectCores() - 1){

  X <- as.matrix(X)
  y <- as.matrix(y)

  param <- expand.grid(C, gamma, degree, coef0)
  m <- nrow(X)
  if(shuffle == TRUE){
    if(is.null(seed) == FALSE){
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
                          .packages = c('manysvms', 'Rcpp'), .options.snow = opts)%dopar%{
    indx_cv <- 1
    accuracy_list <- c()
    for(i in 1:K){
      new_idx_k <- new_idx[indx_cv:(indx_cv+v_size - 1)] #get test dataset
      indx_cv <- indx_cv + v_size
      test_X <- X[new_idx_k, ]
      train_X <- X[-new_idx_k, ]
      test_y <- y[new_idx_k]
      train_y <- y[-new_idx_k]
      mbsvm_model <- mbsvm(train_X, train_y,
                           Ck = param[j, 1]*rep(1, length(unique(y))),
                           kernel = kernel, reg = reg,
                           gamma = param[j, 2], kernel_rect = kernel_rect)
      pred <- predict(mbsvm_model, test_X, test_y)
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
  cat("\nCall:", deparse(call, 0.8 * getOption("width")), "\n", sep="\n")
  cat("Total Parameters:", nrow(param), "\n")
  cat("Best Parameters :",
      "C = ", param[max_idx, 1],
      "\n",
      "gamma = ", param[max_idx, 2],
      "degree = ",param[max_idx, 3],
      "coef0 =", param[max_idx, 4],
      "\n")
  cat("Accuracy :", as.numeric(res[max_idx, 1]),
      "Sd :", as.numeric(res[max_idx, 2]), "\n")
  return(res)
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
