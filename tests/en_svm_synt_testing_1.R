library(MASS)
library(manysvms)

set.seed(100)
n <- 150
X1 <- mvrnorm(n, mu = c(-3, -3), Sigma = diag(1, nrow = 2))
X2 <- mvrnorm(n, mu = c(3, 3), Sigma = diag(1, nrow = 2))
X <- rbind(X1, X2)
y <- rep(c(-1, 1), rep(n, 2))
# model <- en_svm(X, y, solver = "primal", max.steps = 15000)
# print(coef(model))
# plot(model)

model <- en_svm(X, y, solver = "dual", max.steps = 4000)
print(coef(model))
plot(model)
y <- generate_character_label(y)
model_params <- list("C1" = 1, "max.steps" = 8000, "eps" = 0)
cross_validation(sh_svm, X, y, K = 5,
                 metrics = list(accuracy, binaryf1score),
                 metrics_params = list(NULL, list(positive = "Class-1")),
                 model_settings = list("max.steps" = 4000, "eps" = 0))
param_list <- list("C1" = 1:4)

grid_search_cv(sh_svm, X, y, 5, binaryf1score, param_list,
               metrics_params = list("positive" = "Class-2"),
               model_settings = list("max.steps" = 4000, "eps" = 1e-5),
               threads.num = 2)

cross_validation_noisy(en_svm, X, y, y, K = 5,
                       metrics = list(binaryf1score, binaryf1score),
                       metrics_params = list(list(positive = "Class-1"),
                                             list(positive = "Class-2")),
                       model_settings = list("max.steps" = 4000, "eps" = 1e-5))

grid_search_cv_noisy(en_svm, X, y, y, 5, binaryf1score, param_list,
                     metrics_params = list("positive" = "Class-2"),
                     model_settings = list("max.steps" = 4000, "eps" = 1e-5),
                     threads.num = 2)
