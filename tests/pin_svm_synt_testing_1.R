library(MASS)
library(ggplot2)
library(manysvms)

set.seed(100)
n <- 150
X1 <- mvrnorm(n, mu = c(-3, -3), Sigma = diag(1, nrow = 2))
X2 <- mvrnorm(n, mu = c(3, 3), Sigma = diag(1, nrow = 2))
X <- rbind(X1, X2)

y <- rep(c(1, 2), c(n, n))

model1 <- pin_svm(X, y, C = 0.5, tau = 0.5, max.steps = 500)
w <- t(cbind(X, 1)) %*% model1$coef
dataXy <- data.frame(X, y)
dataXy$y <- as.factor(dataXy$y)

ggplot(dataXy, aes(x = X1, y = X2, color = y)) +
  geom_point() +
  geom_abline(slope = -w[1] / w[2],
              intercept = -w[3] / w[2]) +
  geom_abline(slope = -w[1] / w[2],
              intercept = -(-1 + w[3])/w[2]) +
  geom_abline(slope = -w[1] / w[2],
              intercept = -(1 + w[3])/w[2]) +
  theme_bw()

cat(w, "\n")

s <- Sys.time()
model1 <- pin_svm(X, y, C = 1, tau = 0.5,
                          solver = "primal", max.steps = 5000)
e <- Sys.time()
print(e - s)
dataXy <- data.frame(X, y)
dataXy$y <- as.factor(dataXy$y)

ggplot(dataXy, aes(x = X1, y = X2, color = y)) +
  geom_point() +
  geom_abline(slope = -model1$coef[1] / model1$coef[2],
              intercept = -model1$coef[3]/model1$coef[2]) +
  geom_abline(slope = -model1$coef[1] / model1$coef[2],
              intercept = -(-1 + model1$coef[3])/model1$coef[2]) +
  geom_abline(slope = -model1$coef[1] / model1$coef[2],
              intercept = -(1 + model1$coef[3])/model1$coef[2]) +
  theme_bw()

cat(model1$coef, "\n")

tau <- seq(0.1, 0.9, 0.2)
C <- matrix(0, nrow = 17)
for (i in -8:8) {
  C[i] <- 2^(i)
}
param_list <- list("C" = C, "gamma " = C)

s <- Sys.time()
grid_search_cv(pin_svm, X, y, metric = accuracy,
               param_list = param_list, seed = 1234, K = 5,
               max.steps = 500, threads.num = 11,
               solver = "dual", randx = 0.1, batch_size = 1,
               kernel = "rbf")
e <- Sys.time()
print(e - s)
