% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forest.R
\name{isolation_forest}
\alias{isolation_forest}
\title{Detect anomalies using an ensemble of decision trees.}
\usage{
isolation_forest(df, n_trees = 10, subsample_size = 256)
}
\arguments{
\item{df}{A \code{data.frame} with observations.}

\item{n_trees}{The number of trees to fit.}

\item{subsample_size}{The number of examples to use in fitting each tree.}
}
\value{
ret An object with class \code{isolation_forest}.
}
\description{
Detect anomalies using an ensemble of decision trees.
}

