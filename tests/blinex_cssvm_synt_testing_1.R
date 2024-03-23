library(MASS)
library(ggplot2)
library(manysvms)

set.seed(100)
n <- 150
X1 <- mvrnorm(3*n, mu = c(-3, -3), Sigma = diag(1, nrow = 2))
X2 <- mvrnorm(2*n, mu = c(3, 3), Sigma = diag(1, nrow = 2))
X <- rbind(X1, X2)

y <- rep(c(1, 2), c(3*n, 2*n))
s <- Sys.time()
model1 <- blinex_cssvm(X, y , 0.2, kernel = "linear", max.steps = 2000, a = 1, optimizer = nesterov,
                       batch_size = 100)
e <- Sys.time()
print(e - s)
res <- predict(model1, X)
accuracy(res, y)

C <- 2^seq(-8, 8, 2)
param_list <- list("C" = C, "gamma " = 1, "a" = c(1, 2, 3))

s <- Sys.time()
grid_search_cv(blinex_cssvm, X, y, metrics = accuracy,
               param_list = param_list, seed = 1234, K = 5, threads.num = 2,
               model_settings = list("max.steps" = 500, "solver" = "primal",
                                     "randx" = 0.1, "batch_size" = 100,
                                     "kernel" = "rbf"))
e <- Sys.time()
print(e - s)
