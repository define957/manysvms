library(manysvms)

data("iris")

X <- iris[, 1:4]
y <- iris[, 5]

model <- twinsvm_ovr(X, y, kernel = 'linear')
pred_twin_svm_model <- predict(model, X, y)

model <- twinsvm_ovr(X, y, kernel = 'rbf')
pred_twin_svm_model <- predict(model, X, y)
