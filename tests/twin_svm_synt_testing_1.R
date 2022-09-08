# test file

library(manysvms)
library(MASS)
library(ggplot2)


set.seed(112)
n <- 100
sig <- diag(c(0.02, 0.03))
x1 <- mvrnorm(n/2, mu = c(-0.1, 0.1), Sigma = sig)
x2 <- mvrnorm(n/2, mu = c(0.1, -0.1), Sigma = sig)
X <- rbind(x1, x2)
y <- rep(c(1,2), rep(n/2, 2))

twinsvm_model = twinsvm(X, y, reg = 1e-3)

coefs <- coef(twinsvm_model)
pred <- predict(twinsvm_model, X, y)
pred$accuracy
slopes <- c(-twinsvm_model$u1[2]/ twinsvm_model$u1[1], -twinsvm_model$u2[2]/ twinsvm_model$u2[1])
intercepts <- c(-twinsvm_model$b1/twinsvm_model$u1[1], -twinsvm_model$b2/twinsvm_model$u2[1])
plane <- factor(c('1', '2'))
abline_df <- data.frame('slope' = slopes, 'intercept'=intercepts, 'plane' = plane)

xdata <- as.data.frame(cbind(X, y))
colnames(xdata) <- c("x1", "x2", "Class")
xdata$Class <- as.factor(xdata$Class)

acctext <- data.frame(acc = paste('accuracy',pred$accuracy, '%'))
p <- ggplot(data = xdata, aes(x = x1, y = x2, colour = Class)) +
  geom_point(aes(colour = Class), size = 2) +
  geom_abline(data = abline_df, aes(slope = slopes, intercept = intercepts))+
  annotate(geom="text", x=-0.3, y=-0.2, label=paste('Twin-SVM accuracy :',pred$accuracy, '%'))+
  theme_bw()+
  theme(legend.position = c(0.8,0.8))

p
