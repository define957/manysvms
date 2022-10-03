#' Clipping Dual Coordinate Decent Optimizer
#'
#' \code{mbsvm} is an R implementation of multiple birth SVM
#'
#' @author Zhang Jiaqi.
#' @param H,q input matrices.
#' @param lb,ub lower bound and upper bound.
#' @param eps,max.steps error and maximum iterations.
#' @param rcpp speed up your code with Rcpp, default \code{rcpp = TRUE}.
#' @export
clip_dcd_optimizer <- function(H, q, lb, ub, eps = 1e-5,
                               max.steps = 200, rcpp = TRUE){
  if(rcpp == TRUE){
    x <- cpp_clip_dcd_optimizer(H, q, lb, ub, eps = 1e-5, max.steps)
  }else if(rcpp == FALSE){
    x <- r_clip_dcd_optimizer(H, q, lb, ub, eps, max.steps)
  }
  return(x)
}

r_clip_dcd_optimizer <- function(H, q, lb, ub, eps = 1e-5, max.steps = 200){

  # Clipping dual coordinate decent optimizer
  # solve quadratic programming like:
  #       min t(x) %*% H %*% x - t(q) %*% x
  #       s.t lb < x < ub

  H <- as.matrix(H)
  q <- as.matrix(q)
  lb <- as.matrix(lb)
  ub <- as.matrix(ub)

  u <- (lb + ub) / 2
  for(i in 1:max.steps){
    numerator <- (q - t(t(u)%*%H))
    L_idx_val <- numerator / diag(H)
    L_val <- numerator^(2) / diag(H)

    if(max(L_val) < eps){
      break
    }

    idx1 = which(u > lb & L_idx_val < 0)
    idx2 = which(u < ub & L_idx_val >0)
    idx <- unique(c(idx1, idx2))

    if(length(idx) == 0){
      break
    }

    L_val[-idx] = -Inf

    k <- which.max(t(L_val))
    lambda_max <- L_idx_val[k]
    lambda_opt <- max(lb[k] - u[k], min(lambda_max, ub[k] - u[k]))

    u[k] <- u[k] + lambda_opt
  }
  obj_val <- 0.5 * t(u) %*% H %*% u - t(q) %*% u
  clip_dcd_res <- list('x' = u,'iterations' = i, 'objectiv.value' = obj_val)
  return(clip_dcd_res)
}
