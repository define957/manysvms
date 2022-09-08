library(manysvms)
data("iris")


X <- iris[, 1:4]
y <- iris[, 5]

cv.twinsvm_ovr(X, y, K = 5, kernel = 'rbf', gamma = 1/8)
