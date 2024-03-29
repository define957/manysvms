% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bin_to_multi_classifier.R
\name{predict.Classifiers_OVR}
\alias{predict.Classifiers_OVR}
\title{Predict Method of OVR Classifier}
\usage{
\method{predict}{Classifiers_OVR}(object, X, predict_func = predict, ...)
}
\arguments{
\item{object}{a fitted object of class inheriting from \code{Classifiers_OVR}.}

\item{X}{new data for predicting.}

\item{predict_func}{predict function of your model.}

\item{...}{parameters of your model.}
}
\value{
return predict results.
}
\description{
Predict Method of OVR Classifier
}
\author{
Zhang Jiaqi.
}
