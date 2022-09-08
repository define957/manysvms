library(manysvms)
library(MASS)
library(ggplot2)

set.seed(1235)
x1 <- mvrnorm(200, mu = 0, Sigma = 3)
x2 <- 2*x1  + mvrnorm(200, mu = 0, Sigma = 1)

X <- as.matrix(x1)
y <- as.matrix(x2)

m <- eps.svr(X,y, eps=2, C = 1)
dataXy <- as.data.frame(cbind(X, y))
ggplot(data = dataXy, aes(x = X, y = y))+
  geom_point()+
  geom_abline(slope = m$coef, intercept = m$intercept)+
  geom_abline(slope = m$coef, intercept = m$intercept + m$epsilon, linetype=2, color = 'red')+
  geom_abline(slope = m$coef, intercept = m$intercept - m$epsilon, linetype=2, color = 'red')+
  theme_classic()
