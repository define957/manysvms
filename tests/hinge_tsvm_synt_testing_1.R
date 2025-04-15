library(MASS)
library(ggplot2)
library(manysvms)

set.seed(100)
n <- 150
X1 <- mvrnorm(n, mu = c(-3, -3), Sigma = diag(1, nrow = 2))
X2 <- mvrnorm(n, mu = c(3, 3), Sigma = diag(1, nrow = 2))
X <- rbind(X1, X2)

y <- rep(c(1, 2), c(n, n))

s <- Sys.time()
model <- hinge_tsvm(X, y, C1 = 2, C2 = 2, max.steps = 8000)
accuracy(y, predict(model, X))
e <- Sys.time()
print(e - s)
plot(model)


s <- Sys.time()
model <- ls_tsvm(X, y, C1 = 2, C2 = 2)
accuracy(y, predict(model, X))
e <- Sys.time()
print(e - s)
plot(model)

s <- Sys.time()
model <- sh_tsvm(X, y, C1 = 2, C2 = 2, max.steps = 8000)
accuracy(y, predict(model, X))
e <- Sys.time()
print(e - s)
plot(model)
