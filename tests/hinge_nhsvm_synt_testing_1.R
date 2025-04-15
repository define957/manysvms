library(MASS)
library(manysvms)

set.seed(100)
n <- 150
X1 <- mvrnorm(n, mu = c(-3, -3), Sigma = diag(1, nrow = 2))
X2 <- mvrnorm(n, mu = c(3, 3), Sigma = diag(1, nrow = 2))
X <- rbind(X1, X2)
y <- rep(c(2, 1), rep(n, 2))

model <- hinge_nhsvm(X, y, kernel = "linear")
coef(model)
accuracy(y, predict(model, X))
plot(model)

model <- hinge_nhsvm(X, y, kernel = "rbf",
                     reduce_set = X[sample(1:(2*n), 5), , drop = FALSE])
coef(model)
accuracy(y, predict(model, X))
