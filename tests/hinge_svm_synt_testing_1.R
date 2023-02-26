library(MASS)
library(ggplot2)
library(manysvms)

set.seed(100)
X1 <- mvrnorm(100, mu = c(-3, -3), Sigma = diag(1, nrow = 2))
X2 <- mvrnorm(100, mu = c(3, 3), Sigma = diag(1, nrow = 2))
X <- rbind(X1, X2)

y <- rep(c(1, 2), c(100, 100))
s <- Sys.time()
svm_ovr_model <- OVR_Classifier(X, y, hinge_svm, C = 0.1, solver = "primal",
                                max.steps = 5000, batch_size = 10)
e <- Sys.time()
print(e - s)
res <- predict(svm_ovr_model, X)
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
                                 solver = "dual", max.steps = 500)
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

cross_validation(OVR_Classifier, X, y, bin_model = hinge_svm, shuffle = TRUE,
                 metric = accuracy, K = 5, max.steps = 50)

C <- rep(0, 4)
for (i in 1:4) {
  C[i] <- 2^(i)
}
param_list <- list("C" = C, "gamma " = C)

data("iris")
X <- iris[ ,1:4]
y <- iris[ ,5]
grid_search_cv(OVR_Classifier, X, y, metric = accuracy,
               param_list = param_list, seed = 1234, K = 5,
               max.steps = 5000, bin_model = hinge_svm, threads.num = 2,
               solver = "primal",
               kernel = "rbf")
