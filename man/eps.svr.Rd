% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svr.R
\name{eps.svr}
\alias{eps.svr}
\title{epsilon - Support Vector Regression}
\usage{
eps.svr(
  X,
  y,
  eps = 0.1,
  kernel = c("linear", "rbf", "poly"),
  C = 1,
  gamma = 1/ncol(X),
  degree = 3,
  coef0 = 0,
  max.steps = 1000,
  tol = 1e-05,
  rcpp = TRUE
)
}
\arguments{
\item{X, y}{dataset and explained variable.}

\item{eps}{epsilon in the insensitive-loss function (default \code{eps = 0.1}).}

\item{kernel}{kernel function.}

\item{C}{plenty term (default \code{C = 1}).}

\item{gamma}{parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.}

\item{degree}{parameter for polynomial kernel, default: \code{degree = 3}.}

\item{coef0}{parameter for polynomial kernel,  default: \code{coef0 = 0}.}

\item{max.steps}{the number of iterations to solve the optimization problem.}

\item{tol}{tolerance of termination criterion, default: \code{1e-5}.}

\item{rcpp}{speed up your code with Rcpp, default \code{rcpp = TRUE}.}
}
\value{
return eps.svr object.
}
\description{
\code{eps.svr} is an R implementation of epsilon - support vector regression
}
\examples{
library(manysvms)
library(MASS)
library(ggplot2)

set.seed(1235)
obs <- 200
x1 <- seq(-5, 5, length.out = obs)
x2 <- sin(x1) + mvrnorm(obs, mu = 0, Sigma = 0.03)

X <- as.matrix(x1)
y <- as.matrix(x2)

s <- Sys.time()
m <- eps.svr(X, y, eps=0.3, kernel = 'rbf', C = 1,
             gamma = 1, max.steps = 500, rcpp = TRUE)
e <- Sys.time()
print(e - s)
dataXy <- as.data.frame(cbind(X, y))

ggplot(data = dataXy, aes(x = X, y = y))+
  geom_point()+
  geom_line(aes(x=X, y=m$fitted))+
  geom_line(aes(x=X, y=m$fitted + m$epsilon, color = 'red'), show.legend = FALSE)+
  geom_line(aes(x=X, y=m$fitted - m$epsilon, color = 'red'), show.legend = FALSE)+
  theme_bw()
}
\author{
Zhang Jiaqi.
}
