% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sims.R
\name{badmf.sims.exp}
\alias{badmf.sims.exp}
\title{Exponential Simulation}
\usage{
badmf.sims.exp(n, d, K, signal.scale = 1, signal.lshift = 1,
  non.scale = 1, rotate = FALSE, class.equal = TRUE, ind = FALSE)

badmf.sims.exp(n, d, K, signal.scale = 1, signal.lshift = 1,
  non.scale = 1, rotate = FALSE, class.equal = TRUE, ind = FALSE)
}
\arguments{
\item{n}{the number of samples.}

\item{d}{the number of dimensions. The first dimension will be the signal dimension; the remainders noise.}

\item{K}{the number of classes in the dataset.}

\item{signal.scale}{the scaling for the signal dimension. Defaults to \code{1}.}

\item{signal.lshift}{the location shift for the signal dimension between the classes. Defaults to \code{1}.}

\item{non.scale}{the scaling for the non-signal dimensions. Defaults to \code{1}.}

\item{class.equal}{whether the number of samples/class should be equal, with each
class having a prior of 1/K, or inequal, in which each class obtains a prior
of k/sum(K) for k=1:K. Defaults to \code{TRUE}.}

\item{ind}{whether to sample x and y independently. Defaults to \code{FALSE}.}

\item{n}{the number of samples.}

\item{d}{the number of dimensions. The first dimension will be the signal dimension; the remainders noise.}

\item{K}{the number of classes in the dataset.}

\item{signal.scale}{the scaling for the signal dimension. Defaults to \code{1}.}

\item{signal.lshift}{the location shift for the signal dimension between the classes. Defaults to \code{1}.}

\item{non.scale}{the scaling for the non-signal dimensions. Defaults to \code{1}.}

\item{class.equal}{whether the number of samples/class should be equal, with each
class having a prior of 1/K, or inequal, in which each class obtains a prior
of k/sum(K) for k=1:K. Defaults to \code{TRUE}.}

\item{ind}{whether to sample x and y independently. Defaults to \code{FALSE}.}
}
\description{
A function to simulate multi-class data with an Exponential class-mean trend.

A function to simulate multi-class data with an Exponential class-mean trend.
}
\author{
Eric Bridgeford

Eric Bridgeford
}
