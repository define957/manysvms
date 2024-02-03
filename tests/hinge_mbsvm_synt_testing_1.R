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


system.time(model <- hinge_mbsvm(X, y, 1, kernel = "linear", randx = 1, max.steps = 8000))
pred <- predict(model, X)
accuracy(y, pred)
plot(model)

system.time(model <- pin_mbsvm(X, y, 1, kernel = "linear", randx = 1, max.steps = 8000))
pred <- predict(model, X)
accuracy(y, pred)
plot(model)

system.time(model <- ls_mbsvm(X, y, 1, kernel = "linear"))
pred <- predict(model, X)
accuracy(y, pred)
plot(model)

system.time(model <- ramp_mbsvm(X, y, 1, kernel = "linear", randx = 1, max.steps = 8000,
                                cccp.steps = 10, s = 0.01))
pred <- predict(model, X)
accuracy(y, pred)
plot(model)
