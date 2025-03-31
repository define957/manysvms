library(MASS)
library(manysvms)

set.seed(100)
n <- 150
X1 <- mvrnorm(n, mu = c(-3, -3), Sigma = diag(1, nrow = 2))
X2 <- mvrnorm(n, mu = c(3, 3), Sigma = diag(1, nrow = 2))
X <- rbind(X1, X2)
y <- rep(c(-1, 1), rep(n, 2))
model <- roboss_svm(X, y, solver = "primal", max.steps = 1500)
print(coef(model))
plot(model)

y <- generate_character_label(y)
model_params <- list("C" = 1, "max.steps" = 8000, "eps" = 0)
cross_validation(roboss_svm, X, y, K = 5,
                 metrics = list(accuracy, binaryf1score),
                 metrics_params = list(NULL, list(positive = "Class-1")),
                 model_settings = list("max.steps" = 4000, "eps" = 0))
param_list <- list("C" = 1:4)

grid_search_cv(roboss_svm, X, y, 5, binaryf1score, param_list,
               metrics_params = list("positive" = "Class-2"),
               model_settings = list("max.steps" = 4000, "eps" = 1e-5),
               threads.num = 2)

cross_validation_noisy(roboss_svm, X, y, y, K = 5,
                       metrics = list(binaryf1score, binaryf1score),
                       metrics_params = list(list(positive = "Class-1"),
                                             list(positive = "Class-2")),
                       model_settings = list("max.steps" = 4000,
                                             "solver" = "primal",
                                             "kernel" = "linear"))

cross_validation_noisy(roboss_svm, X, y, y, K = 5,
                       metrics = list(binaryf1score, binaryf1score),
                       metrics_params = list(list(positive = "Class-1"),
                                             list(positive = "Class-2")),
                       model_settings = list("max.steps" = 4000,
                                             "solver" = "primal",
                                             "kernel" = "rbf"))

reduce_set <- cbind(X[sample(1:nrow(X), 30), ], 1)

cross_validation_noisy(roboss_svm, X, y, y, K = 5,
                       metrics = list(binaryf1score, binaryf1score),
                       metrics_params = list(list(positive = "Class-1"),
                                             list(positive = "Class-2")),
                       model_settings = list("max.steps" = 4000,
                                             "solver" = "primal",
                                             "kernel" = "rbf",
                                             "reduce_set" = reduce_set))

grid_search_cv_noisy(roboss_svm, X, y, y, 5, binaryf1score, param_list,
                     metrics_params = list("positive" = "Class-2"),
                     model_settings = list("max.steps" = 4000, "eps" = 1e-5),
                     threads.num = 2)
