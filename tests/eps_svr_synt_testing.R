library(manysvms)


set.seed(1234)
n <- 100
x <- sort(rnorm(n, 0, 2))
y <- sin(3*x)/x + rnorm(n, 0, 0.1)

model <- eps_svr(x, y, kernel = "rbf", max.steps = 800)
pred <- predict(model, x)

plot(x, y, col = "red")
lines(x, pred)
