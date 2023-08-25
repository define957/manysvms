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

system.time(model <- hinge_mbsvm(X, y, 1, kernel = "linear", randx = 1, max.steps = 800000))
pred <- predict(model, X)
accuracy(y, pred)
plot(model)

X <- glass[, 1:9]
X <- apply(X, 2, scale)
y <- glass[, 10]

X <- Ecoli[, 1:7]
X <- apply(X, 2, scale)
y <- Ecoli[, 8]

C <- 2^seq(-2, 12)
param_list <- list("C" = C, "gamma" = 2^seq(-10, 4))
grid_search_cv(hinge_mbsvm, X, y, 5, accuracy, param_list, kernel = "rbf", randx = 1, seed = 1213,
               max.steps = 15000)

grid_search_cv(OVR_Classifier, X, y, 5, accuracy, param_list, kernel = "rbf", randx = 1, seed = 1213,
               max.steps = 800, bin_model = hinge_svm, values = TRUE)
