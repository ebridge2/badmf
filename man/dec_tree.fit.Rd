% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dec_tree.R
\name{dec_tree.fit}
\alias{dec_tree.fit}
\title{Bayesian Decision Tree Fit}
\usage{
dec_tree.fit(formuler, data = NULL, d = NULL, alpha = NULL,
  method = "classification", depth.max = 5L, size = 1L,
  debug = FALSE, ...)
}
\arguments{
\item{formuler}{ravioli ravioli give me the formuoli.}

\item{data}{the data associated with the formuler. Note: if you want an intercept, you must
add it ahead of time.}

\item{d}{the number of features to subsample at each node. Defaults to \code{NULL}, which tries every feature.}

\item{alpha}{the prior parameters for the feature probabilities. A \code{[p]} vector. If \code{NULL}, samples uniformly.
Defaults to \code{NULL}.}

\item{method}{whether you want "classification" or "regression". Defaults to \code{"classification"}.}

\item{depth.max}{the maximum allowed tree depth. Defaults to \code{5L}.}

\item{size}{the minimum allowed number of samples for an individual node. Defaults to \code{1L}.}

\item{debug}{whether to save the predictors and responses that are categorized. Defaults to \code{FALSE}.}

\item{...}{trailing arguments.}
}
\value{
an object of class \code{dec.tree.class} containing the following:
\item{\code{tree}}{the decision tree.}
\item{\code{X}}{The training predictors.}
\item{\code{Y}}{the training responses.}
\item{\code{d}}{d the number of features subsampled at each node.}
\item{\code{alpha}}{the sampling distribution for the features. A \code{[p]} vector.}
\item{\code{depth.max}}{the maximum allowed tree depth.}
\item{\code{size}}{the maximum allowed tree depth.}
\item{\code{debug}}{whether to save the predictors and responses that are categorized.}
}
\description{
Fit a Bayesian Decision Tree with a \code{stats}-like formula frontend interface.
}
\author{
Eric Bridgeford
}
