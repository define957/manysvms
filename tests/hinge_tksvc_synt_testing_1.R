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

model <- hinge_tksvc(X, y, C = 2^(8), epsilon = 0, max.steps = 80000)
accuracy(y, predict(model, X))
id <- 1
w_pos <- model$coef_pos[, id]
w_neg <- model$coef_neg[, id]
plot(X[, 1], X[, 2], col = "red")
abline(-w_pos[3]/w_pos[2], -w_pos[1]/w_pos[2],
       lty = 1, col = "red")
abline(-(w_pos[3]+1)/w_pos[2], -w_pos[1]/w_pos[2],
       lty = 2, col = "red")
abline(-w_neg[3]/w_neg[2], -w_neg[1]/w_neg[2],
       lty = 1, col = "blue")
abline(-(w_neg[3]-1)/w_neg[2], -w_neg[1]/w_neg[2],
       lty = 2, col = "blue")

model <- hinge_tksvc(X, y, C = 2^(8), epsilon = 0, max.steps = 80000)

accuracy(y, predict(model, X))
