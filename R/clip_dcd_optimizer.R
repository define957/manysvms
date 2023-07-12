#' Clipping Dual Coordinate Decent Optimizer
#'
#' @author Zhang Jiaqi.
#' @param H,q input matrices.
#' @param lb,ub lower bound and upper bound.
#' @param eps,max.steps error and maximum iterations.
#' @param u init solution.
#' @param ... redundant parameters.
#' @return return results list
#' @export
clip_dcd_optimizer <- function(H, q, lb, ub,
                                 eps = 1e-5, max.steps = 200,
                                 u = (lb + ub) / 2, ...) {

  # Clipping dual coordinate decent optimizer
  # solve quadratic programming like:
  #       min t(x) %*% H %*% x - t(q) %*% x
  #       s.t lb < x < ub

  H <- as.matrix(H)
  q <- as.matrix(q)
  lb <- as.matrix(lb)
  ub <- as.matrix(ub)
  diagH <- diag(H)
  n <- nrow(H)
  u_mat <- t(matrix(u, n, n))
  Hui <- H*u_mat
  Hu <- H%*%u
  for (t in 1:max.steps) {
    numerator <- q - Hu
    L_idx_val <- numerator / diagH
    L_val <- numerator*L_idx_val
    idx <- which((u > lb & L_idx_val < 0) | (u < ub & L_idx_val > 0))
    L_val_max <- max(L_val[idx])
    if (L_val_max < eps) {
      break
    }
    if (length(idx) == 0) {
      break
    }
    k <- which(L_val == L_val_max)[1]
    lambda_max <- L_idx_val[k]
    lambda_opt <- max(lb[k] - u[k], min(lambda_max, ub[k] - u[k]))
    u[k] <- u[k] + lambda_opt
    # for (i in 1:n) {
    #   Hu[i] <- Hu[i] - Hui[i, k] + H[i, k]*u[k]
    #   Hui[i, k] <- H[i, k]*u[k]
    # }
    Huik <- H[, k]*u[k]
    Hu <- Hu - Hui[, k] + Huik
    Hui[, k] <- Huik
  }
  obj_val <- 0.5 * t(u) %*% H %*% u - t(q) %*% u
  clip_dcd_res <- list('x' = u,'iterations' = t, 'objectiv.value' = obj_val)
  return(clip_dcd_res)
}
