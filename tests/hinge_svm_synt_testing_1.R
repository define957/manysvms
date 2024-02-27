library(MASS)
library(manysvms)

set.seed(100)
n <- 150
X1 <- mvrnorm(n, mu = c(-3, -3), Sigma = diag(1, nrow = 2))
X2 <- mvrnorm(n, mu = c(3, 3), Sigma = diag(1, nrow = 2))
X <- rbind(X1, X2)
y <- rep(c(-1, 1), rep(n, 2))
model <- hinge_svm(X, y, solver = "primal", max.steps = 800)
print(coef(model))
plot(model)

model <- hinge_svm(X, y, solver = "dual")
print(coef(model))
plot(model)
y <- generate_character_label(y)
model_params <- list("C" = 1, "max.steps" = 8000, "eps" = 0)
cross_validation(hinge_svm, X, y, K = 5,
                 metrics = list(accuracy, binaryf1score),
                 metrics_params = list(NULL, positive = "Class-1"))
