library(manysvms)
data("iris")


X <- iris[, 1:4]
y <- iris[, 5]
gamma <- c(8, 4, 2, 1, 1/4)

res <- cv.twinsvm_ovr(X, y, K = 5, kernel = 'rbf', reg = 1e-3,
                      gamma = gamma, threads.num = 2)
