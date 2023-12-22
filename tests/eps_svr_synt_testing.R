library(manysvms)

set.seed(1234)
n <- 100
x <- sort(rnorm(n))
y <- 3*x + 1 + rnorm(n, 0, 1)

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

model <- hinge_tsvr(x, y, C1 = 256,
                    epsilon1 = 5, epsilon2 = 5,max.steps = 800000, eps = -0.1)
c1 <- model$coef1
c2 <- model$coef2
plot(x, y, col = "red")
pred <- predict(model, x)
lines(x, pred)
abline(c1[2], c1[1])
abline(c2[2], c2[1])

abline(c1[2] + 5, c1[1], col = "red")
abline(c2[2] - 5, c2[1], col = "red")
