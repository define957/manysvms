library(manysvms)

set.seed(43)
n <- 100
x <- sort(rnorm(n))
y <- 3*x + 1 + rnorm(n, 0, 3)

epsilon <- 30
model <- hinge_tsvr(x, y, kernel = "linear", epsilon1 = epsilon, max.steps = 80000,
                    C1 = 800000, eps = 0)

lag1 <- model$solver.res$lag1
lag2 <- model$solver.res$lag2

idx1 <- which(lag1 != 0)
idx2 <- which(lag2 != 0)

plot(x, y, xlim = c(-3, 3), ylim = c(-50, 50))
abline(a = model$coef1[2], b = model$coef1[1], col = "green")
abline(a = model$coef2[2], b = model$coef2[1], col = "blue")

abline(a = model$coef1[2] + epsilon, b = model$coef1[1], col = "green", lty = 2)
abline(a = model$coef2[2] - epsilon, b = model$coef2[1], col = "blue", lty = 2)

points(x[idx1], y[idx1], col = "green")
points(x[idx2], y[idx2], col = "blue")


points(x, y + epsilon)
points(x, y - epsilon)
