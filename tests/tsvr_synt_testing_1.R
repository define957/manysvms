library(manysvms)

set.seed(43)
n <- 100
x <- sort(rnorm(n))
y <- 3*x + 1 + rnorm(n, 0, 3)

epsilon <- 0.5
model <- hinge_tsvr(x, y, kernel = "linear", epsilon1 = epsilon, max.steps = 80000,
                    C1 = 1, eps = 0)
plot(model)

model <- hinge_eps_tsvr(x, y, kernel = "linear", epsilon1 = epsilon, max.steps = 80000,
                        C1 = 1, eps = 0)
plot(model)

model <- sh_eps_tsvr(x, y, kernel = "linear", epsilon1 = epsilon, max.steps = 80000,
                        C1 = 1, eps = 0)
plot(model)
