library(MASS)
library(ggplot2)
library(manysvms)

set.seed(112)
n <- 300
n_c <- 3
sig <- diag(c(0.02, 0.03))

x1 <- mvrnorm(n/n_c, mu = c(-0.6, 0), Sigma = sig)
x2 <- mvrnorm(n/n_c, mu = c(0.6, 0), Sigma = sig)
x3 <- mvrnorm(n/n_c, mu = c(0, 0.6), Sigma = sig)

X <- rbind(x1, x2, x3)
y <- rep(c(1, 2, 3), rep(n/n_c, n_c))
s <- Sys.time()
svm_ovr_model <- OVR_Classifier(X, y, hinge_svm, C = 100, solver = "primal",
                                max.steps = 5000, batch_size = 50)
e <- Sys.time()
print(e - s)
model <- svm_ovr_model$classifier_list[[3]]

res <- predic.OVR_Classifier(svm_ovr_model, X, y, predict.hinge_svm)
dataXy <- data.frame(X, y)
dataXy$y <- as.factor(dataXy$y)

p <- ggplot(dataXy, aes(x = X1, y = X2, color = y)) + geom_point()
color_set <- c("red", "green", "blue")
for (i in 1:n_c) {
  model <- svm_ovr_model$classifier_list[[i]]
  p <- p + geom_abline(slope = -model$coef[1] / model$coef[2],
                       intercept = -model$coef[3]/model$coef[2], colour = color_set[i]) +
            geom_abline(slope = -model$coef[1] / model$coef[2],
                        intercept = -(-1 + model$coef[3])/model$coef[2], colour = color_set[i]) +
            geom_abline(slope = -model$coef[1] / model$coef[2],
                        intercept = -(1 + model$coef[3])/model$coef[2], colour = color_set[i])
}
p <- p + theme_bw()
p
