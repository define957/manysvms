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
model1 <- blinex_cssvm(X, y, 1, "linear", max.steps = 200)
e <- Sys.time()
print(e - s)
res <- predict(model1, X)
accuracy(res, y)

dataXy <- data.frame(X, y)
dataXy$y <- as.factor(dataXy$y)

ggplot(dataXy, aes(x = X1, y = X2, color = y)) +
  geom_point() +
  geom_abline(slope = -model1$coef[1] / model1$coef[2],
              intercept = -model1$coef[3] / model1$coef[2]) +
  geom_abline(slope = -model1$coef[1] / model1$coef[2],
              intercept = -(-1 + model1$coef[3])/model1$coef[2]) +
  geom_abline(slope = -model1$coef[1] / model1$coef[2],
              intercept = -(1 + model1$coef[3])/model1$coef[2]) +
  theme_bw()

cat(model1$coef, "\n")

C <- 2^seq(-8, 8, 2)
param_list <- list("C" = C, "gamma " = C, "a" = c(1, 2, 3))

s <- Sys.time()
grid_search_cv(blinex_cssvm, X, y, metrics = accuracy,
               param_list = param_list, seed = 1234, K = 5,
               max.steps = 500, threads.num = 2,
               solver = "primal", randx = 1, batch_size = 1,
               kernel = "rbf")
e <- Sys.time()
print(e - s)
