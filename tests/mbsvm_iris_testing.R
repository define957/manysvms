library(manysvms)
data('iris')


X <- iris[1:150, 1:4]
y <- iris[1:150, 5]

model1 <- mbsvm(X, y, kernel = 'linear')
pred <-predict(model1, X, y)

model <- mbsvm(X, y, kernel = 'rbf', kernel_rect = 0.45, gamma = 1)
pred <-predict(model, X, y)

cv.mbsvm(X, y, kernel = 'linear', reg = 0.05, threads.num = 1)

cv.mbsvm(X, y, kernel = 'rbf', reg = 0.05, threads.num = 1)
