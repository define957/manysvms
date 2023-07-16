library(MASS)
library(manysvms)

set.seed(100)
n <- 150
X1 <- mvrnorm(n, mu = c(-3, -3), Sigma = diag(1, nrow = 2))
X2 <- mvrnorm(n, mu = c(3, 3), Sigma = diag(1, nrow = 2))
X <- rbind(X1, X2)
y <- rep(c(-1, 1), rep(n, 2))
model <- hinge_svm(X, y, solver = "primal")
print(coef(model))
plot(model)

model <- hinge_svm(X, y, solver = "dual")
print(coef(model))
plot(model)

cross_validation(OVR_Classifier, X, y, bin_model = hinge_svm, shuffle = TRUE,
                 metrics = list(accuracy), K = 5, max.steps = 200, values = T)

C <- rep(-8, 8)
for (i in 0:17) {
  C[i] <- 2^(i)
}
param_list <- list("C" = C, "gamma " = C)

s <- Sys.time()
grid_search_cv(hinge_svm, X, y, metrics = accuracy,
               param_list = param_list, seed = 1234, K = 5,
               max.steps = 500, threads.num = 2,
               solver = "primal", randx = 0.1, batch_size = 1,
               kernel = "rbf")
e <- Sys.time()
print(e - s)
