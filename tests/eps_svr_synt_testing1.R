library(manysvms)
library(MASS)
library(ggplot2)

set.seed(1235)
x1 <- seq(-5, 5, length.out=200)#mvrnorm(200, mu = 0, Sigma = 0.1)
x2 <- sin(x1)  + mvrnorm(200, mu = 0, Sigma = 0.03)

X <- as.matrix(x1)
y <- as.matrix(x2)

m <- eps.svr(X,y, eps=0.1, kernel = 'rbf', C = 1, gamma = 1, max.steps = 3000)
dataXy <- as.data.frame(cbind(X, y))
ggplot(data = dataXy, aes(x = X, y = y))+
  geom_point()+
  geom_abline(slope = m$coef, intercept = m$intercept)+
  geom_abline(slope = m$coef, intercept = m$intercept + m$epsilon, linetype=2, color = 'red')+
  geom_abline(slope = m$coef, intercept = m$intercept - m$epsilon, linetype=2, color = 'red')+
  theme_classic()

ggplot(data = dataXy, aes(x = X, y = y))+
  geom_point()+
  geom_line(aes(x=X, y=m$fitted))+
  geom_line(aes(x=X, y=m$fitted + m$epsilon, color = 'red'))+
  geom_line(aes(x=X, y=m$fitted - m$epsilon, color = 'red'))+
  theme_classic()
