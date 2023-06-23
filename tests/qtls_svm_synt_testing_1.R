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
model1 <- qtls_svm(X, y , 0.2, "linear", max.steps = 2000, a = 1, optimizer = pegasos,
                       batch_size = 50)
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
param_list <- list("C" = C, "gamma " = 1, "a" = -(1:50))

s <- Sys.time()
grid_search_cv(qtls_svm, X, y, metrics = accuracy,
               param_list = param_list, seed = 1234, K = 5,
               max.steps = 500, threads.num = 2,
               solver = "primal", randx = 0.1,
               kernel = "rbf")
e <- Sys.time()
print(e - s)
