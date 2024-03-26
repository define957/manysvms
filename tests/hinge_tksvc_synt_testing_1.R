library(MASS)
library(ggplot2)
library(manysvms)

set.seed(100)
n <- 150
X1 <- mvrnorm(n, mu = c(3, 0), Sigma = diag(0.2, nrow = 2))
X2 <- mvrnorm(n, mu = c(-3, 0), Sigma = diag(0.2, nrow = 2))
X3 <- mvrnorm(n, mu = c(0, 3), Sigma = diag(0.2, nrow = 2))
X <- rbind(X1, X2, X3)

y <- rep(c(1, 2, 3), rep(n, 3))

model <- hinge_tksvc(X, y, C1 = 2^(8), C2 = 2^(8), epsilon = 0, max.steps = 80000)
accuracy(y, predict(model, X))
id <- 1
w_pos <- model$coef_pos[, id]
w_neg <- model$coef_neg[, id]
plot(X[, 1], X[, 2], col = "red")
abline(-w_pos[3]/w_pos[2], -w_pos[1]/w_pos[2],
       lty = 1, col = "red")
abline(-(w_pos[3] + 1)/w_pos[2], -w_pos[1]/w_pos[2],
       lty = 2, col = "red")
abline(-w_neg[3]/w_neg[2], -w_neg[1]/w_neg[2],
       lty = 1, col = "blue")
abline(-(w_neg[3] - 1)/w_neg[2], -w_neg[1]/w_neg[2],
       lty = 2, col = "blue")

model <- hinge_tksvc(X, y, C1 = 2^(8), C2 = 2^(8), epsilon = 0, max.steps = 80000)

accuracy(y, predict(model, X))


model_params <- list("C1" = 1, "C2" = 4, "max.steps" = 4000, "eps" = 0)
cross_validation(hinge_tksvc, X, y, K = 5,
                 metrics = list(accuracy),
                 model_settings = model_params)

param_list <- list("C1" = 1:4, "C2" = 1:4, "epsilon" = c(0, 0.1, 0.2))

grid_search_cv(hinge_tksvc, X, y, 5, accuracy, param_list,
               model_settings = list("max.steps" = 4000, "eps" = 1e-5),
               threads.num = 2)

cross_validation_noisy(hinge_tksvc, X, y, y, K = 5,
                       metrics = accuracy,
                       model_settings = list("max.steps" = 4000, "eps" = 1e-5))

grid_search_cv_noisy(hinge_tksvc, X, y, y, 5, accuracy, param_list,
                     model_settings = list("max.steps" = 4000, "eps" = 1e-5),
                     threads.num = 2)

# Pin-TKSVC

model_settings <- list("C1" = 1, "C2" = 4, "tau1" = 0.1, "tau2" = 0.1,
                     "max.steps" = 4000, "eps" = 0)
cross_validation(pin_tksvc, X, y, K = 5,
                 metrics = list(accuracy),
                 model_settings = model_settings)

param_list <- list("C1" = 1:4, "epsilon" = c(0, 0.1, 0.2),
                   "tau1" = c(0.1, 0.2, 0.3))

grid_search_cv(pin_tksvc, X, y, 5, accuracy, param_list,
               model_settings = list("max.steps" = 4000, "eps" = 1e-5),
               threads.num = 2)

cross_validation_noisy(pin_tksvc, X, y, y, K = 5,
                       metrics = accuracy,
                       model_settings = list("max.steps" = 4000, "eps" = 1e-5))

grid_search_cv_noisy(pin_tksvc, X, y, y, 5, accuracy, param_list,
                     model_settings = list("max.steps" = 4000, "eps" = 1e-5),
                     threads.num = 2)
