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

system.time(model <- hinge_mbsvm(X, y, 1, kernel = "linear"))
pred <- predict(model, X)

plot(model)
