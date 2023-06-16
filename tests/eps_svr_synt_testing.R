library(manysvms)


set.seed(1234)
n <- 100
x <- sort(runif(n, -5, 5))
y <- sin(3*x)/x + rnorm(n, 0, 0.1)

model <- eps_svr(x, y, kernel = "rbf", max.steps = 800)
pred <- predict(model, x)

param_list <- list("C" = 1:4,
                   "epsilon" = c(0.1, 0.2, 0.5))

met <- list("mse" = mean_squared_error, "mae" = mean_absolute_error)
res <- grid_search_cv(eps_svr, x, y, metrics = met, threads.num = 2,
               param_list = param_list, max.steps = 1000, solver = "dual"
               )
res
plot(x, y, col = "red")
lines(x, pred)

