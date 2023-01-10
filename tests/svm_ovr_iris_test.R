library(manysvms)

data('iris')
X <- iris[, 1:4]
y <- iris[, 5]

model <- svm_ovr(X, y, max.steps = 1000)
pred <- predict(model, X[1:10, ], y[1:10])

C <- c(2^(-4), 2^(-2), 2^(0), 2^(2), 2^(4), 2^(6))
rbf_gamma <- c(2^(-4), 2^(-2), 2^(0), 2^(2), 2^(4), 2^(6))
res <- cv.svm_ovr(X, y, K = 5, C = C, gamma = rbf_gamma, kernel = 'rbf',
                  max.steps = 1000, threads.num = 2)
