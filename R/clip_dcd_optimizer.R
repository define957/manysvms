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
  Hui <- H*matrix(u, n, n, byrow = T)
  Hu <- H %*% u
  ub_u <- ub - u
  lb_u <- lb - u
  for (t in 1:max.steps) {
    numerator <- q - Hu
    L_idx_val <- numerator / diagH
    L_val <- numerator*L_idx_val
    idx <- which((u > lb & L_idx_val < 0) | (u < ub & L_idx_val > 0))
    if (length(idx) == 0) {
      break
    }
    L_val_max <- max(L_val[idx])
    if (L_val_max < eps) {
      break
    }
    k_list <- which(L_val == L_val_max)
    lambda_max_list <- L_idx_val[k_list]
    lambda_opt <- 0
    k <- 1
    for (i in 1:length(k_list)) {
      ktemp <- k_list[i]
      lambda_opt_temp <- max(lb_u[ktemp],
                             min(lambda_max_list[i], ub_u[ktemp]))
      if (abs(lambda_opt) < abs(lambda_opt_temp)) {
        k <- ktemp
        lambda_opt <- lambda_opt_temp
      }
    }
    u[k] <- u[k] + lambda_opt
    Huik <- H[, k]*u[k]
    Hu <- Hu - Hui[, k] + Huik
    Hui[, k] <- Huik
    lb_u[k] <- lb[k] - u[k]
    ub_u[k] <- ub[k] - u[k]
  }
  obj_val <- 0.5 * t(u) %*% Hu - t(q) %*% u
  clip_dcd_res <- list('x' = u,'iterations' = t, 'objectiv.value' = obj_val)
  return(clip_dcd_res)
}
