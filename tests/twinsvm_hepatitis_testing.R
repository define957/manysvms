library(manysvms)

hepatitis_n<-na.omit(hepatitis)

X <- hepatitis_n[, 1:19]
y <- hepatitis_n[, 20]

model <- twinsvm(X, y, kernel = 'linear')
pred <-predict(model, X, y)

model <- twinsvm(X, y, kernel = 'rbf')
pred <-predict(model, X, y)
