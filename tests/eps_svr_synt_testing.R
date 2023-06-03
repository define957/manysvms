library(manysvms)


set.seed(1234)
n <- 100
x <- sort(runif(n, -5, 5))
y <- sin(3*x)/x + rnorm(n, 0, 0.1)

model <- eps_svr(x, y, kernel = "rbf", max.steps = 800)
pred <- predict(model, x)

param_list <- list("C" = 1:4,
                   "epsilon" = c(0.1, 0.2, 0.5))

# grid_search_cv(eps_svr, x, y, metric = mean_squared_error,
#                param_list = param_list, max.steps = 1000)

plot(x, y, col = "red")
lines(x, pred)
