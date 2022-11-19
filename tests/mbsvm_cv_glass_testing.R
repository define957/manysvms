library(manysvms)


X <- glass[, 1:9]
y <- glass[, 10]

model1 <- mbsvm(X, y, kernel = 'linear',)
pred <-predict(model1, X, y)

model <- mbsvm(X, y, kernel = 'rbf', kernel_rect = 1, gamma = 1)
pred <-predict(model, X, y)

cv.mbsvm(X, y, kernel = 'linear', reg = 0.05, threads.num = 1)
