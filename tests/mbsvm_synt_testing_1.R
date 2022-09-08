# test file
library(manysvms)
library(MASS)
library(ggplot2)


set.seed(112)
n <- 300
n_c <- 3
sig <- diag(c(0.02, 0.03))

x1 <- mvrnorm(n/n_c, mu = c(-0.6, 0), Sigma = sig)
x2 <- mvrnorm(n/n_c, mu = c(0.6, 0), Sigma = sig)
x3 <- mvrnorm(n/n_c, mu = c(0, 0.6), Sigma = sig)

X <- rbind(x1, x2, x3)
y <- rep(c(1,2,3), rep(n/n_c, n_c))

mbsvm_model = mbsvm(X, y, reg = 1e-3)

coefs <- mbsvm_model$coef
intercept <- mbsvm_model$intercept
slopes <-  - coefs[1, ] / coefs[2, ]
intercepts <- - intercept / coefs[2, ]
plane <- factor(c('1', '2', '3'))

grid <- expand.grid(seq(min(X[, 1]), max(X[, 1]),length.out=20),
                    seq(min(X[, 2]), max(X[, 2]),length.out=20))
Z <- predict(mbsvm_model, grid, rep(1, 400))$predict
contour_data <- data.frame(grid, Z)
colnames(contour_data) <- c("x1", "x2", "Class")
contour_data$Class <- as.factor(contour_data$Class)
abline_df <- data.frame('slope' = slopes, 'intercept'=intercepts, 'plane' = plane)

xdata <- as.data.frame(cbind(X, y))
colnames(xdata) <- c("x1", "x2", "Class")
xdata$Class <- as.factor(xdata$Class)

pred <- predict(mbsvm_model, X, y)

p <- ggplot(contour_data, aes(x = x1, y = x2)) +
  geom_tile(aes(fill = Class, alpha = 0.2), show.legend = FALSE) +
  geom_point(data = xdata, aes(x = x1, y=x2, color = Class, shape=Class), size = 3) +
  geom_abline(data = abline_df, aes(slope = slopes, intercept = intercepts, colour = plane), show.legend = FALSE) +
  annotate(geom="text", x=-0.7, y=0.75, label=paste('MBSVM accuracy :',100*round(pred$accuracy, 4), '%')) +
  theme_bw()
p

