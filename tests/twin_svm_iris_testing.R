library(manysvms)
library(ggplot2)
data("iris")

m <- 100
X <- iris[0:m, 1:2]
y <- iris[0:m, 5]

s <- Sys.time()
model <- twinsvm(X, y, kernel = 'linear')
e <- Sys.time()

print(e-s)

coef(model)

pred <-predict(model, X, y)

slopes <- c(-model$u1[1]/ model$u1[2], -model$u2[1]/model$u2[2])
intercepts <- c(-model$b1/model$u1[2], -model$b2/model$u2[2])

plane <- factor(c('1', '2'))
abline_df <- data.frame('slope' = slopes, 'intercept'=intercepts, 'plane' = plane)

X <- as.data.frame(X)
dataXy <- cbind(X, y)

ggplot(dataXy, aes(x=Sepal.Length, y = Sepal.Width, colour = y))+
  geom_abline(data = abline_df, aes(slope = slopes, intercept = intercepts))+
  geom_point()+
  theme_bw()
