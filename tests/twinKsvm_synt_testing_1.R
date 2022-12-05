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

s <- Sys.time()
model <- twinKsvm(X, y, kernel = 'rbf', eps = 0.3, kernel_rect = 1)
e <- Sys.time()
print(e - s)
pred <- predict(model, X, y)

s <- Sys.time()
cv.twinKsvm(X, y,  kernel = 'rbf', eps = 0.3, kernel_rect = 1, threads.num = 2)
e <- Sys.time()
print(e - s)

X <- manysvms::glass[,1:10]
y <- manysvms::glass[,10]
s <- Sys.time()
cv.twinKsvm(X, y, K = 10, gamma = 1 / 4, kernel = 'rbf',
            eps = 0.2, kernel_rect = 0.5, threads.num = 2)
e <- Sys.time()
print(e - s)

data("iris")
X <- iris[, 1:4]
y <- iris[, 5]
s <- Sys.time()
cv.twinKsvm(X, y, gamma = 1 / 4, kernel = 'linear',
            eps = 0.1, kernel_rect = 1, threads.num = 2)
e <- Sys.time()
print(e - s)
