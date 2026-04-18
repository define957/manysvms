library(manysvms)

set.seed(1234)
n <- 100
x <- sort(rnorm(n))
y <- 3*x + 1 + rnorm(n, 0, 3)

epsilon <- 1.5
model <- hinge_tsvr(x, y, kernel = "linear", epsilon1 = epsilon, max.steps = 80000,
                    C1 = 1, eps = 0)
plot(model)
mean_squared_error(y, predict(model, x))

model <- hinge_eps_tsvr(x, y, kernel = "linear", epsilon1 = epsilon, max.steps = 80000,
                        C1 = 1, eps = 0)
plot(model)
mean_squared_error(y, predict(model, x))

model <- sh_eps_tsvr(x, y, kernel = "linear", epsilon1 = epsilon, max.steps = 80000,
                        C1 = 1, eps = 0)
plot(model)
mean_squared_error(y, predict(model, x))
