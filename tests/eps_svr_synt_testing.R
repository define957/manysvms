library(manysvms)

set.seed(1234)
n <- 100
x <- sort(rnorm(n))
y <- 3*x + 1 + rnorm(n, 0, 0.1)

model <- ls_svr(x, y, C = 1, kernel = "linear", max.steps = 8000, solver = "primal", batch_size = 100,
                optimizer = pegasos, projection = FALSE)
pred <- predict(model, x)
plot(x, y, col = "red")
lines(x, pred)

model <- eps_svr(x, y, C = 1, epsilon = 0.1, kernel = "linear", max.steps = 8000)
pred <- predict(model, x)
plot(x, y, col = "red")
lines(x, pred)

model <- bls_svr(x, y, C = 1, b = 1, kernel = "linear", max.steps = 8000)
pred <- predict(model, x)
plot(x, y, col = "red")
lines(x, pred)
