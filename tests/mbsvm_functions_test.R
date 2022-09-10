library(manysvms)
library(MASS)
set.seed(112)
n <- 300
n_c <- 3
sig <- diag(c(0.02, 0.03))

x1 <- mvrnorm(n/n_c, mu = c(-0.6, 0), Sigma = sig)
x2 <- mvrnorm(n/n_c, mu = c(0.6, 0), Sigma = sig)
x3 <- mvrnorm(n/n_c, mu = c(0, 0.6), Sigma = sig)
#x4 <- mvrnorm(n/n_c, mu = c(1, 0.6), Sigma = sig)

X <- rbind(x1, x2, x3)
y <- rep(c(1,2,3), rep(n/n_c, n_c))
s <- Sys.time()
mbsvm_model = mbsvm(X, y, reg = 1e-3, kernel = 'rbf', kernel_rect = 0.1, rcpp = TRUE)
e <- Sys.time()
print(e-s)
pred <- predict(mbsvm_model, X, y)

plot(mbsvm_model, xlab = 'x1', ylab = 'x2')
