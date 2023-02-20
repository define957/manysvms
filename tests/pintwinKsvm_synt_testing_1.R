library(manysvms)
library(MASS)
library(ggplot2)


set.seed(112)
n <- 300
n_c <- 3
sig <- diag(c(0.02, 0.03))

x1 <- MASS::mvrnorm(n/n_c, mu = c(0, 1), Sigma = sig)
x2 <- MASS::mvrnorm(n/n_c, mu = c(1, 1), Sigma = sig)
x3 <- MASS::mvrnorm(n/n_c, mu = c(0.5, 0.5), Sigma = diag(c(0.01,0.01)))

X <- rbind(x1, x2, x3)
y <- rep(c(1,2,3), rep(n/n_c, n_c))

data("iris")
X <- iris[, 1:4]
y <- iris[, 5]
s <- Sys.time()
cv.pintwinKsvm(X, y,tau = 0.5, gamma = 1 / 4, kernel = 'linear',
            eps = 0.1, kernel_rect = 1, threads.num = 2)
e <- Sys.time()
print(e - s)
