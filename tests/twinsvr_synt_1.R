library(ggplot2)
library(manysvms)
library(MASS)

set.seed(1234)
obs <- 300
x1 <- seq(-5, 5, length.out = obs)
x2 <- sin(x1)  + mvrnorm(obs, mu = 0, Sigma = 0.03)

X <- as.matrix(x1)
y <- as.matrix(x2)

model <- twinsvr(X, y, reg = 1e-7, kernel = 'rbf', kernel_rect = 1)

df1 <- as.data.frame(cbind(x1, y))
df2 <- as.data.frame(cbind(x1, model$fitted))
colnames(df2) <- c("x", "y")
colnames(df1) <- c("x", "y")
ggplot(data = df1, aes(x = x, y = y))+
  geom_point()+
  geom_line(data = df2, aes(x = x, y = y))+
  theme_bw()
