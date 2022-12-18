library(manysvms)
library(MASS)
library(ggplot2)


set.seed(112)
n <- 300
n_c <- 3
sig <- diag(c(0.02, 0.03))

x1 <- mvrnorm(n/n_c, mu = c(0, 1), Sigma = sig)
x2 <- mvrnorm(n/n_c, mu = c(1, 1), Sigma = sig)
x3 <- mvrnorm(n/n_c, mu = c(0.5, 0.5), Sigma = diag(c(0.01,0.01)))

X <- rbind(x1, x2, x3)
y <- rep(c(1,2,3), rep(n/n_c, n_c))

s <- Sys.time()
model <- ramptwinKsvm(X, y, kernel = 'rbf',kernel_rect = 1,  reg = 0.2,
                      eps = 0.1, step_cccp = 10, max.steps = 300)
e <- Sys.time()
print(e - s)
pred <- predict(model, X, y)

cv.ramptwinKsvm(X, y, K = 5, threads.num = 2)
