% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_selection.R
\name{cross_validation_noisy}
\alias{cross_validation_noisy}
\title{K-Fold Cross Validation with Noisy (Simulation Only)}
\usage{
cross_validation_noisy(
  model,
  X,
  y,
  y_noisy,
  K = 5,
  metrics,
  predict_func = predict,
  pipeline = NULL,
  metrics_params = NULL,
  predict_params = NULL,
  model_settings = NULL,
  transy = FALSE,
  model_seed = NULL
)
}
\arguments{
\item{model}{your model.}

\item{X, y}{dataset and label.}

\item{y_noisy}{label with label noise.}

\item{K}{number of folds.}

\item{metrics}{this parameter receive a metric function.}

\item{predict_func}{this parameter receive a function for predict.}

\item{pipeline}{preprocessing pipline.}

\item{metrics_params}{set parameters for each metrics (need a list).}

\item{predict_params}{set parameters for each predict method (need a list).}

\item{model_settings}{set parameters for model (need a list).}

\item{transy}{apply transforms defined in `pipeline` on y, default FALSE.}

\item{model_seed}{random_seed for model.}
}
\value{
return a metric matrix
}
\description{
\code{cross_validation_noisy} function use noisy data for training,
then calculates the average and standard deviation of your metric
using clean samples
}
\author{
Zhang Jiaqi.
}
