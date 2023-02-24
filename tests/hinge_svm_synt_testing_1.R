library(MASS)
library(ggplot2)
library(manysvms)

set.seed(100)
X1 <- mvrnorm(100, mu = c(-3, -3), Sigma = diag(1, nrow = 2))
X2 <- mvrnorm(100, mu = c(3, 3), Sigma = diag(1, nrow = 2))
X <- rbind(X1, X2)

y <- rep(c(1, 2), c(100, 100))
s <- Sys.time()
svm_ovr_model <- OVR_Classifier(X, y, hinge_svm, C = 200, solver = "primal",
                                max.steps = 10000, batch_size = 200)
e <- Sys.time()
print(e - s)
res <- predic.OVR_Classifier(svm_ovr_model, X, y, predict.hinge_svm)
model1 <- svm_ovr_model$classifier_list[[1]]

dataXy <- data.frame(X, y)
dataXy$y <- as.factor(dataXy$y)

ggplot(dataXy, aes(x = X1, y = X2, color = y)) +
  geom_point() +
  geom_abline(slope = -model1$coef[1] / model1$coef[2],
              intercept = -model1$coef[3] / model1$coef[2]) +
  geom_abline(slope = -model1$coef[1] / model1$coef[2],
              intercept = -(-1 + model1$coef[3])/model1$coef[2]) +
  geom_abline(slope = -model1$coef[1] / model1$coef[2],
              intercept = -(1 + model1$coef[3])/model1$coef[2]) +
  theme_bw()

cat(model1$coef, "\n")

s <- Sys.time()
svm_ovr_model1 <- OVR_Classifier(X, y, hinge_svm, C = 1,
                                 solver = "dual", max.steps = 2000)
e <- Sys.time()
print(e - s)
model1 <- svm_ovr_model1$classifier_list[[1]]
model1$coef <- t(cbind(X, 1)) %*% model1$coef
dataXy <- data.frame(X, y)
dataXy$y <- as.factor(dataXy$y)

ggplot(dataXy, aes(x = X1, y = X2, color = y)) +
  geom_point() +
  geom_abline(slope = -model1$coef[1] / model1$coef[2],
              intercept = -model1$coef[3]/model1$coef[2]) +
  geom_abline(slope = -model1$coef[1] / model1$coef[2],
              intercept = -(-1 + model1$coef[3])/model1$coef[2]) +
  geom_abline(slope = -model1$coef[1] / model1$coef[2],
              intercept = -(1 + model1$coef[3])/model1$coef[2]) +
  theme_bw()

cat(model1$coef, "\n")
