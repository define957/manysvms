library(ggplot2)
library(manysvms)
library(MASS)

set.seed(1234)
obs <- 300
x1 <- seq(-5, 5, length.out = obs)
x2 <- sin(x1)  + mvrnorm(obs, mu = 0, Sigma = 0.03)

X <- as.matrix(x1)
y <- as.matrix(x2)

coef_ <- lssvr(X, y, kernel = 'linear')
