library(manysvms)
library(MASS)
library(ggplot2)

set.seed(1235)
obs <- 200
x1 <- seq(-5, 5, length.out = obs)
x2 <- sin(x1)  + mvrnorm(obs, mu = 0, Sigma = 0.03)

X <- as.matrix(x1)
y <- as.matrix(x2)

s <- Sys.time()
m <- eps.svr(X,y, eps=0.3, kernel = 'rbf', C = 1,
             gamma = 1, max.steps = 500, rcpp = FALSE)
e <- Sys.time()
print(e - s)
dataXy <- as.data.frame(cbind(X, y))

ggplot(data = dataXy, aes(x = X, y = y))+
  geom_point()+
  geom_line(aes(x=X, y=m$fitted))+
  geom_line(aes(x=X, y=m$fitted + m$epsilon, color = 'red'), show.legend = FALSE)+
  geom_line(aes(x=X, y=m$fitted - m$epsilon, color = 'red'), show.legend = FALSE)+
  theme_bw()
