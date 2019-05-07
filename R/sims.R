#' Mean Difference Simulation
#'
#' A function for simulating data in which a difference in the means is present only in a subset of dimensions, and equal covariance.
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param rotate whether to apply a random rotation to the mean and covariance. With random rotataion matrix \code{Q}, \code{mu = Q*mu}, and \code{S = Q*S*Q}. Defaults to \code{FALSE}.
#' @param priors the priors for each class. If \code{NULL}, class priors are all equal. If not null, should be \code{|priors| = K}, a length \code{K} vector for \code{K} classes. Defaults to \code{NULL}.
#' @param K the number of classes. Defaults to \code{2}.
#' @param md the magnitude of the difference in the means in the specified subset of dimensions. Ddefaults to \code{1}.
#' @param subset the dimensions to have a difference in the means. Defaults to only the first dimension. \code{max(subset) < d}. Defaults to \code{c(1)}.
#' @param offdiag the off-diagonal elements of the covariance matrix. Should be < 1. \code{S_{ij} = offdiag} if \code{i != j}, or 1 if \code{i == j}. Defaults to \code{0}.
#' @param s the scaling parameter of the covariance matrix. S_{ij} = scaling*1 if i == j, or scaling*offdiag if i != j. Defaults to \code{1}.
#' @return A list of class \code{simulation} with the following:
#' \item{X}{\code{[n, d]} the \code{n} data points in \code{d} dimensions as a matrix.}
#' \item{Y}{\code{[n]} the \code{n} labels as an array.}
#' \item{mus}{\code{[d, K]} the \code{K} class means in \code{d} dimensions.}
#' \item{Sigmas}{\code{[d, d, K]} the \code{K} class covariance matrices in \code{d} dimensions.}
#' \item{priors}{\code{[K]} the priors for each of the \code{K} classes.}
#' \item{simtype}{The name of the simulation.}
#' \item{params}{Any extraneous parameters the simulation was created with.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("sims", package = "badmf")}
#'
#' @author Eric Bridgeford
#' @examples
#' library(badmf)
#' data <- badmf.sims.mean_diff(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
badmf.sims.mean_diff <- function(n, d, rotate=FALSE, priors=NULL, K=2, md=1, subset=c(1), offdiag=0, s=1) {
  if (is.null(priors)) {
    priors <- array(1/K, dim=c(K))
  } else if (length(priors) != K) {
    stop(sprintf("You have specified %d priors for %d classes.", length(priors), K))
  } else if (sum(priors) != 1) {
    stop(sprintf("You have passed invalid priors. The sum(priors) should be 1; yours is %.3f", sum(priors)))
  }
  if (max(subset) > d) {
    stop(sprintf("Specified a difference in dimension %d; maximum should be %d.", max(subset), d))
  }
  mus <- array(0, dim=c(d, K))
  for (i in 1:K) {
    mus[subset] <- (i - 1)*md  # scale the difference in the means accordingly
  }
  S <- array(offdiag, dim=c(d, d))
  diag(S) <- 1  # identity variance
  S <- s*S  # scale accordingly

  S <- array(unlist(replicate(K, S, simplify=FALSE)), dim=c(d, d, K))

  if (rotate) {
    res <- badmf.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }
  # simulate from GMM
  sim <- badmf.sims.sim_gmm(mus, S, n, priors)
  return(structure(list(X=sim$X, Y=as.factor(sim$Y), mus=mus, Sigmas=S, priors=sim$priors, simtype="Mean Difference",
                        params=list(K=K, md=md, subset=subset, offdiag=offdiag, s=s)), class="simulation"))
}

#' Toeplitz Simulation
#'
#' A function for simulating data in which the covariance is a non-symmetric toeplitz matrix.
#'
#' @importFrom abind abind
#' @importFrom stats toeplitz
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param rotate whether to apply a random rotation to the mean and covariance. With random rotataion matrix \code{Q}, \code{mu = Q*mu}, and \code{S = Q*S*Q}. Defaults to \code{FALSE}.
#' @param priors the priors for each class. If \code{NULL}, class priors are all equal. If not null, should be \code{|priors| = K}, a length \code{K} vector for \code{K} classes. Defaults to \code{NULL}.
#' @param D1 the dimensionality for the non-equal covariance terms. Defaults to \code{10}.
#' @param b a scaling parameter for the means. Defaults to \code{0.4}.
#' @param rho the scaling of the covariance terms, should be < 1. Defaults to \code{0.5}/
#' @return A list of class \code{simulation} with the following:
#' \item{X}{\code{[n, d]} the \code{n} data points in \code{d} dimensions as a matrix.}
#' \item{Y}{\code{[n]} the \code{n} labels as an array.}
#' \item{mus}{\code{[d, K]} the \code{K} class means in \code{d} dimensions.}
#' \item{Sigmas}{\code{[d, d, K]} the \code{K} class covariance matrices in \code{d} dimensions.}
#' \item{priors}{\code{[K]} the priors for each of the \code{K} classes.}
#' \item{simtype}{The name of the simulation.}
#' \item{params}{Any extraneous parameters the simulation was created with.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("sims", package = "badmf")}
#'
#' @author Eric Bridgeford
#' @examples
#' library(badmf)
#' data <- badmf.sims.toep(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
badmf.sims.toep <- function(n, d, rotate=FALSE, priors=NULL, D1=10, b=0.4, rho=0.5) {
  K <- 2
  if (is.null(priors)) {
    priors <- array(1/K, dim=c(K))
  } else if (length(priors) != K) {
    stop(sprintf("You have specified %d priors for %d classes.", length(priors), K))
  } else if (sum(priors) != 1) {
    stop(sprintf("You have passed invalid priors. The sum(priors) should be 1; yours is %.3f", sum(priors)))
  }
  if (rho >= 1) {
    stop(sprintf("rho should be < 1; user specified %.3f", rho))
  }
  cT <- rho^(0:(D1 - 1))
  A <- toeplitz(cT)
  K1 <- sum(A)

  cT <- rho^(0:(d-1))
  A <- toeplitz(cT)
  K <- sum(A)

  mudelt <- (K1 * b^2/K)^0.5/2

  mu0 <- array(mudelt*c(1, -1), dim=c(d))
  mus <- abind(mu0, -mu0, along=2)
  S <- abind(A, A, along=3)

  if (rotate) {
    res <- badmf.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }
  # simulate from GMM
  sim <- badmf.sims.sim_gmm(mus, S, n, priors)
  return(structure(list(X=sim$X, Y=as.factor(sim$Y), mus=mus, Sigmas=S, priors=sim$priors, simtype="Toeplitz",
                        params=list(D1=D1, b=b, rho=rho)), class="simulation"))
}

#' Quadratic Discriminant Toeplitz Simulation
#'
#' A function for simulating data generalizing the Toeplitz setting, where each class has a different covariance matrix. This results in a Quadratic Discriminant.
#'
#' @importFrom abind abind
#' @importFrom stats toeplitz
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param rotate whether to apply a random rotation to the mean and covariance. With random rotataion matrix \code{Q}, \code{mu = Q*mu}, and \code{S = Q*S*Q}. Defaults to \code{FALSE}.
#' @param priors the priors for each class. If \code{NULL}, class priors are all equal. If not null, should be \code{|priors| = K}, a length \code{K} vector for \code{K} classes. Defaults to \code{NULL}.
#' @param D1 the dimensionality for the non-equal covariance terms. Defaults to \code{10}.
#' @param b a scaling parameter for the means. Defaults to \code{0.4}.
#' @param rho the scaling of the covariance terms, should be < 1. Defaults to \code{0.5}.
#' @return A list of class \code{simulation} with the following:
#' \item{X}{\code{[n, d]} the \code{n} data points in \code{d} dimensions as a matrix.}
#' \item{Y}{\code{[n]} the \code{n} labels as an array.}
#' \item{mus}{\code{[d, K]} the \code{K} class means in \code{d} dimensions.}
#' \item{Sigmas}{\code{[d, d, K]} the \code{K} class covariance matrices in \code{d} dimensions.}
#' \item{priors}{\code{k
#'
#' A simulation for the reversed random trunk experiment, in which the maximal covariant directions are the same as the directions with the maximal mean diffe[K]} the priors for each of the \code{K} classes.}
#' \item{simtype}{The name of the simulation.}
#' \item{params}{Any extraneous parameters the simulation was created with.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("sims", package = "badmf")}
#'
#' @author Eric Bridgeford
#' @examples
#' library(badmf)
#' data <- badmf.sims.qdtoep(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
badmf.sims.qdtoep <- function(n, d, rotate=FALSE, priors=NULL, D1=10, b=0.4, rho=0.5) {
  K <- 2
  if (is.null(priors)) {
    priors <- array(1/K, dim=c(K))
  } else if (length(priors) != K) {
    stop(sprintf("You have specified %d priors for %d classes.", length(priors), K))
  } else if (sum(priors) != 1) {
    stop(sprintf("You have passed invalid priors. The sum(priors) should be 1; yours is %.3f", sum(priors)))
  }
  if (rho >= 1) {
    stop(sprintf("rho should be < 1; user specified %.3f", rho))
  }
  if (D1 > d) {
    stop(sprintf("User has specified more dimensions for non-equal cov terms. D1 is %d, yet d is %d", D1, d))
  }
  tR <- rho^(0:(D1 - 1))
  A <- toeplitz(tR)
  K1 <- sum(A)

  tR <- rho^(0:(d-1))
  A0 <- toeplitz(tR)
  K <- sum(A0)

  mudelt <- (K1 * b^2/K)^0.5/2

  mu0 <- array(mudelt*c(1, -1), dim=c(d))
  Q <- badmf.sims.rotation(d)
  mu1 <- -Q %*% (mu0 + 0.1)
  mus <- abind(mu0, mu1, along=2)
  A1 <- Q %*% A0 %*% t(Q)
  S <- abind(A0, A1, along=3)

  if (rotate) {
    res <- badmf.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }
  # simulate from GMM
  sim <- badmf.sims.sim_gmm(mus, S, n, priors)
  return(structure(list(X=sim$X, Y=as.factor(sim$Y), mus=mus, Sigmas=S, priors=sim$priors, simtype="Quadratic Toeplitz",
                        params=list(D1=D1, b=b, rho=rho)), class="simulation"))
}

#' Random Trunk
#'
#' A simulation for the random trunk experiment, in which the maximal covariant dimensions are the reverse of the maximal mean differences.
#' @importFrom abind abind
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param rotate whether to apply a random rotation to the mean and covariance. With random rotataion matrix \code{Q}, \code{mu = Q*mu}, and \code{S = Q*S*Q}. Defaults to \code{FALSE}.
#' @param priors the priors for each class. If \code{NULL}, class priors are all equal. If not null, should be \code{|priors| = K}, a length \code{K} vector for \code{K} classes. Defaults to \code{NULL}.
#' @param b scalar for mu scaling. Default to \code{4}.
#' @param K number of classes, should be <4. Defaults to \code{2}.
#' @param maxvar the maximum covariance between the two classes. Defaults to \code{100}.
#' @param maxvar.outlier the maximum covariance for the outlier points. Defaults to \code{maxvar*5}.
#' @return A list of class \code{simulation} with the following:
#' \item{X}{\code{[n, d]} the \code{n} data points in \code{d} dimensions as a matrix.}
#' \item{Y}{\code{[n]} the \code{n} labels as an array.}
#' \item{mus}{\code{[d, K]} the \code{K} class means in \code{d} dimensions.}
#' \item{Sigmas}{\code{[d, d, K]} the \code{K} class covariance matrices in \code{d} dimensions.}
#' \item{priors}{\code{[K]} the priors for each of the \code{K} classes.}
#' \item{simtype}{The name of the simulation.}
#' \item{params}{Any extraneous parameters the simulation was created with.}
#' \item{robust}{If robust is not false, a list containing \code{inlier} a boolean array indicating which points are inliers, \code{s.outlier} the covariance structure of outliers, and \code{mu.outlier} the means of the outliers.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("sims", package = "badmf")}
#'
#' @author Eric Bridgeford
#' @examples
#' library(badmf)
#' data <- badmf.sims.rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
badmf.sims.rtrunk <- function(n, d, rotate=FALSE, priors=NULL, b=4, K=2, maxvar=100) {
  if (is.null(priors)) {
    priors <- array(1/K, dim=c(K))
  } else if (length(priors) != K) {
    stop(sprintf("You have specified %d priors for %d classes.", length(priors), K))
  } else if (sum(priors) != 1) {
    stop(sprintf("You have passed invalid priors. The sum(priors) should be 1; yours is %.3f", sum(priors)))
  }
  mu1 <- b/sqrt(0:(d-1)*2 + 1)
  if (K == 2) {
    mus <- abind(mu1, -mu1, along=2)
  } else if (K == 3) {
    mus <- abind(mu1, 0*mu1, -mu1, along=2)
  } else if (K == 4) {
    mus <- abind(mu1, b/(seq(from=d, to=1, by=-1)), b/(seq(from=1, to=d, by=1)), -mu1, along=2)
  }
  s <- diag(d)
  diag(s) <- maxvar/sqrt(seq(from=d, to=1, by=-1))

  S <- array(unlist(replicate(K, s, simplify=FALSE)), dim=c(d, d, K))

  if (rotate) {
    res <- badmf.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }

  # simulate from GMM
  sim <- badmf.sims.sim_gmm(mus, S, n, priors)
  structure(list(X=sim$X, Y=as.factor(sim$Y), mus=mus, Sigmas=S, priors=sim$priors, simtype="Random Trunk",
                 params=list(b=b, K=K)), class="simulation")
}

#' Linear Simulation
#'
#' A function to simulate multi-class data with a linear class-mean trend. The signal dimension is the dimension carrying all of the
#' between-class difference, and the non-signal dimensions are noise.
#'
#' @import abind
#' @param n the number of samples.
#' @param d the number of dimensions. The first dimension will be the signal dimension; the remainders noise.
#' @param K the number of classes in the dataset.
#' @param signal.scale the scaling for the signal dimension. Defaults to \code{1}.
#' @param signal.lshift the location shift for the signal dimension between the classes. Defaults to \code{1}.
#' @param non.scale the scaling for the non-signal dimensions. Defaults to \code{1}.
#' @param class.equal whether the number of samples/class should be equal, with each
#' class having a prior of 1/K, or inequal, in which each class obtains a prior
#' of k/sum(K) for k=1:K. Defaults to \code{TRUE}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @author Eric Bridgeford
#' @export
badmf.sims.linear <- function(n, d, K, signal.scale=1, signal.lshift=1, non.scale=1, rotate=FALSE, class.equal=TRUE, ind=FALSE) {
  priors <- gen.sample.labels(K, class.equal=class.equal)
  S <- diag(d)
  S[1, 1] <- signal.scale
  S[-c(1), -c(1)] <- non.scale
  S <- abind(lapply(1:K, function(i) {
    S
  }), along=3)
  mu <- c(1, rep(0, d-1))
  mu[1] <- signal.lshift
  mus <- abind(lapply(1:K, function(i) {
    mu*i*signal.lshift
  }), along=2)

  if (rotate) {
    res <- badmf.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
    rotate=res$Q
  }

  sim <- badmf.sims.sim_gmm(mus, S, n, priors=priors)
  X <- sim$X; Y <- factor(sim$Y)
  return(list(X=X, Y=as.factor(Y), mus=mus, Sigmas=S, priors=priors, simtype="Linear",
              params=list(signal.scale=signal.scale, signal.lshift=signal.lshift, non.scale=non.scale,
                          rotate=rotate, class.equal=class.equal, ind=ind)))
}

#' Exponential Simulation
#'
#' A function to simulate multi-class data with an Exponential class-mean trend.
#'
#' @import abind
#' @param n the number of samples.
#' @param d the number of dimensions. The first dimension will be the signal dimension; the remainders noise.
#' @param K the number of classes in the dataset.
#' @param signal.scale the scaling for the signal dimension. Defaults to \code{1}.
#' @param signal.lshift the location shift for the signal dimension between the classes. Defaults to \code{1}.
#' @param non.scale the scaling for the non-signal dimensions. Defaults to \code{1}.
#' @param class.equal whether the number of samples/class should be equal, with each
#' class having a prior of 1/K, or inequal, in which each class obtains a prior
#' of k/sum(K) for k=1:K. Defaults to \code{TRUE}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @author Eric Bridgeford
#' @export
badmf.sims.exp <- function(n, d, K, signal.scale=1, signal.lshift=1, non.scale=1, rotate=FALSE, class.equal=TRUE, ind=FALSE) {
  priors <- gen.sample.labels(K, class.equal=class.equal)
  S <- diag(d)
  S[1, 1] <- signal.scale; S[-c(1), -c(1)] <- non.scale
  S <- abind(lapply(1:K, function(i) {
    S
  }), along=3)
  mu <- c(1, rep(0, d-1))
  mu[1] <- signal.lshift
  mus <- abind(lapply(1:K, function(i) {
    mu*exp(i)*signal.lshift
  }), along=2)

  if (rotate) {
    res <- badmf.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
    rotate=res$Q
  }

  sim <- badmf.sims.sim_gmm(mus, S, n, priors=priors)
  X <- sim$X; Y <- factor(sim$Y)
  return(list(X=X, Y=as.factor(Y), mus=mus, Sigmas=S, priors=priors, simtype="Linear",
              params=list(signal.scale=signal.scale, signal.lshift=signal.lshift,
                          non.scale=non.scale, rotate=rotate, class.equal=class.equal,
                          ind=ind)))
}

#' Spread Simulation
#'
#' A function to simulate data with the same mean that spreads as class id increases.
#'
#' @import abind
#' @param n the number of samples.
#' @param d the number of dimensions.
#' @param K the number of classes in the dataset. Defaults to \code{2}.
#' @param signal.scale the scaling for the signal dimension. Defaults to \code{1}.
#' @param class.equal whether the number of samples/class should be equal, with each
#' class having a prior of 1/K, or inequal, in which each class obtains a prior
#' of k/sum(K) for k=1:K. Defaults to \code{TRUE}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @examples
#' library(badmf)
#' sim <- badmf.sims.spread(100, 3, 2)
#' @author Eric Bridgeford
#' @export
badmf.sims.fat_tails <- function(n, d, K=2, signal.scale=1, rotate=FALSE, class.equal=TRUE, ind=FALSE) {
  priors <- gen.sample.labels(K, class.equal=class.equal)
  S <- signal.scale*diag(d)
  S <- abind(lapply(1:K, function(i) {
    i^2*S
  }), along=3)
  mus <- array(0, dim=c(d, K))

  if (rotate) {
    res <- badmf.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }

  sim <- badmf.sims.sim_gmm(mus, S, n, priors=priors)
  X <- sim$X; Y <- factor(sim$Y)

  if (ind) {
    X <- badmf.sims.sim_gmm(mus, S, n, priors=priors)$X
  }

  return(list(X=X, Y=as.factor(Y), mus=mus, Sigmas=S, priors=priors, simtype="Fat Tails",
              params=list(signal.scale=signal.scale, rotate=rotate, class.equal=class.equal,
                          ind=ind)))
}

#' Cross Simulation
#'
#' A function to simulate data with a crossed layout.
#'
#' @import abind
#' @param n the number of samples.
#' @param d the number of dimensions.
#' @param K the number of classes in the dataset.
#' @param signal.scale the scaling for the signal dimension. Defaults to \code{10}.
#' @param non.scale the scaling for the non-signal dimensions. Defaults to \code{1}.
#' @param class.equal whether the number of samples/class should be equal, with each
#' class having a prior of 1/K, or inequal, in which each class obtains a prior
#' of k/sum(K) for k=1:K. Defaults to \code{TRUE}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @examples
#' library(badmf)
#' sim <- badmf.sims.cross(100, 3)
#' X <- sim$X; Y <- sim$Y
#' @author Eric Bridgeford
#' @export
badmf.sims.cross <- function(n, d, K=2, signal.scale=10, non.scale=1, rotate=FALSE, class.equal=TRUE, ind=FALSE) {
  priors <- gen.sample.labels(K, class.equal=class.equal)
  S <- diag(d)
  S <- abind(lapply(1:K, function(i) {
    S.class <- S
    S.class[i,i] <- signal.scale
    S.class
  }), along=3)
  mus <- array(0, dim=c(d, K))

  if (rotate) {
    res <- badmf.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }

  sim <- badmf.sims.sim_gmm(mus, S, n, priors=priors)
  X <- sim$X; Y <- factor(sim$Y)

  if (ind) {
    X <- badmf.sims.sim_gmm(mus, S, n, priors=priors)$X
  }

  return(list(X=X, Y=as.factor(Y), mus=mus, Sigmas=S, priors=priors, simtype="Cross",
              params=list(signal.scale=signal.scale, rotate=rotate, class.equal=class.equal,
                          ind=ind)))
}
#' Beta Simulation
#'
#' Signal Dimension is beta-distributed; remainder are Uniform 0, 1
#'
#' @import abind
#' @param n the number of samples.
#' @param d the number of dimensions.
#' @param class.equal whether the number of samples/class should be equal, with each
#' class having a prior of 1/K, or inequal, in which each class obtains a prior
#' of k/sum(K) for k=1:K. Defaults to \code{TRUE}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @author Eric Bridgeford
#' @export
badmf.sims.beta <- function(n, d, class.equal=TRUE, ind=FALSE) {
  K <- 5
  priors <- gen.sample.labels(K, class.equal=class.equal)
  class <- 1:K
  # sample n points from 1:K for assignment of class label
  Y <- sample(class, size=n, replace=TRUE, prob=priors)
  # initialize X
  X <- array(NaN, dim=c(n, d))
  alpha <-c(5, 5, 5, 2, 1)
  beta <- c(1, 2, 5, 5, 5)
  for (i in 1:K) {
    # fill the top dimension with beta distributed data
    X[Y == i,1] <- rbeta(sum(Y == i), alpha[i], beta[i])
  }
  if (d > 1) {
    X[, 2:d] <- runif(n*(d - 1))  # fill the rest with uniforms
  }
  if (ind) {
    Y <- sample(class, size=n, replace=TRUE, prob=priors)
  }
  return(list(X=X, Y=factor(Y), priors=priors, simtype="Beta",
              params=list(alpha=alpha, beta=beta)))
}

#' A helper function for simulating sample labels
#' @param K the number of classes
#' @param class.equal whether the number of samples/class should be equal, with each
#' class having a prior of 1/K, or inequal, in which each class obtains a prior
#' of k/sum(K) for k=1:K. Defaults to \code{TRUE}.
gen.sample.labels <- function(K, class.equal=TRUE) {
  if (isTRUE(class.equal)) {
    priors <- array(1/K, dim=c(K)) # prior is just 1/K for all k
  } else {
    priors <- (1:K)/sum(1:K) # prior is k/sum(K) for all k
  }
  # return class priors
  return(priors)
}

#' Reverse Random Trunk
#'
#' A simulation for the reversed random trunk experiment, in which the maximal covariant directions are the same as the directions with the maximal mean difference.
#' @importFrom abind abind
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param robust the number of outlier points to add, where outliers have opposite covariance of inliers. Defaults to \code{FALSE}, which will not add any outliers.
#' @param rotate whether to apply a random rotation to the mean and covariance. With random rotataion matrix \code{Q}, \code{mu = Q*mu}, and \code{S = Q*S*Q}. Defaults to \code{FALSE}.
#' @param priors the priors for each class. If \code{NULL}, class priors are all equal. If not null, should be \code{|priors| = K}, a length \code{K} vector for \code{K} classes. Defaults to \code{NULL}.
#' @param b scalar for mu scaling. Default to \code{4}.
#' @param K number of classes, should be <4. Defaults to \code{2}.
#' @param maxvar the maximum covariance between the two classes. Defaults to \code{100}.
#' @return A list of class \code{simulation} with the following:
#' \item{X}{\code{[n, d]} the \code{n} data points in \code{d} dimensions as a matrix.}
#' \item{Y}{\code{[n]} the \code{n} labels as an array.}
#' \item{mus}{\code{[d, K]} the \code{K} class means in \code{d} dimensions.}
#' \item{Sigmas}{\code{[d, d, K]} the \code{K} class covariance matrices in \code{d} dimensions.}
#' \item{priors}{\code{[K]} the priors for each of the \code{K} classes.}
#' \item{simtype}{The name of the simulation.}
#' \item{params}{Any extraneous parameters the simulation was created with.}
#' \item{robust}{If robust is not false, a list containing \code{inlier} a boolean array indicating which points are inliers, \code{s.outlier} the covariance structure of outliers, and \code{mu.outlier} the means of the outliers.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("sims", package = "badmf")}
#'
#' @author Eric Bridgeford
#' @examples
#' library(badmf)
#' data <- badmf.sims.rev_rtrunk(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
badmf.sims.rev_rtrunk <- function(n, d, robust=FALSE, rotate=FALSE, priors=NULL, b=4, K=2, maxvar=b^3, maxvar.outlier=maxvar^3) {
  if (is.null(priors)) {
    priors <- array(1/K, dim=c(K))
  } else if (length(priors) != K) {
    stop(sprintf("You have specified %d priors for %d classes.", length(priors), K))
  } else if (sum(priors) != 1) {
    stop(sprintf("You have passed invalid priors. The sum(priors) should be 1; yours is %.3f", sum(priors)))
  }
  mu1 <- b/sqrt(0:(d-1)*2 + 1)
  if (K == 2) {
    mus <- abind(mu1, -mu1, along=2)
  } else if (K == 3) {
    mus <- abind(mu1, 0*mu1, -mu1, along=2)
  } else if (K == 4) {
    mus <- abind(mu1, b/(seq(from=d, to=1, by=-1)), b/(seq(from=1, to=d, by=1)), -mu1, along=2)
  }
  if (!identical(robust, FALSE)) {
    mus[(d/2):d,] <- 0
  }
  S <- diag(d)
  diag(S) <- maxvar/sqrt(seq(from=1, to=d, by=1))

  S <- array(unlist(replicate(K, S, simplify=FALSE)), dim=c(d, d, K))

  if (rotate) {
    res <- badmf.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }
  # simulate from GMM
  sim <- badmf.sims.sim_gmm(mus, S, n, priors)
  return.list <- list(mus=mus, Sigmas=S, priors=sim$priors, simtype="Reverse Random Trunk",
                      params=list(b=b, K=K))
  if (identical(robust, FALSE)) {
    inlier <- !logical(n)
    return.list$X <- sim$X; return.list$Y <- sim$Y
  } else {
    s.outlier <- diag(d)
    diag(s.outlier) <- maxvar.outlier/sqrt(seq(from=d, to=1, by=-1))
    s.outlier <- array(unlist(replicate(K, s.outlier, simplify=FALSE)), dim=c(d, d, K))
    mu.outlier <- -mus
    sim.outlier <- badmf.sims.sim_gmm(mu.outlier, s.outlier, robust, priors)

    X <- rbind(sim$X, sim.outlier$X)
    Y <- c(sim$Y, sim.outlier$Y)
    inlier <- c(!logical(n), logical(robust))
    # randomly reorder X and Y
    reord <- sample(1:length(Y))
    return.list$X <- X[reord,]; return.list$Y <- Y[reord]; return.list$robust <- list()
    return.list$robust$inlier <- inlier[reord]
    return.list$robust$sigma.outlier <- s.outlier; return.list$robust$mu.outlier <- mu.outlier
  }
  return(structure(return.list, class="simulation"))
}

#' Stacked Cigar
#'
#' A simulation for the stacked cigar experiment.
#' @importFrom abind abind
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param rotate whether to apply a random rotation to the mean and covariance. With random rotataion matrix \code{Q}, \code{mu = Q*mu}, and \code{S = Q*S*Q}. Defaults to \code{FALSE}.
#' @param priors the priors for each class. If \code{NULL}, class priors are all equal. If not null, should be \code{|priors| = K}, a length \code{K} vector for \code{K} classes. Defaults to \code{NULL}.
#' @param a scalar for all of the mu1 but 2nd dimension. Defaults to \code{0.15}.
#' @param b scalar for 2nd dimension value of mu2 and the 2nd variance term of S. Defaults to \code{4}.
#' @return A list of class \code{simulation} with the following:
#' \item{X}{\code{[n, d]} the \code{n} data points in \code{d} dimensions as a matrix.}
#' \item{Y}{\code{[n]} the \code{n} labels as an array.}
#' \item{mus}{\code{[d, K]} the \code{K} class means in \code{d} dimensions.}
#' \item{Sigmas}{\code{[d, d, K]} the \code{K} class covariance matrices in \code{d} dimensions.}
#' \item{priors}{\code{[K]} the priors for each of the \code{K} classes.}
#' \item{simtype}{The name of the simulation.}
#' \item{params}{Any extraneous parameters the simulation was created with.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("sims", package = "badmf")}
#'
#' @author Eric Bridgeford
#' @examples
#' library(badmf)
#' data <- badmf.sims.cigar(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
badmf.sims.cigar <- function(n, d, rotate=FALSE, priors=NULL, a=0.15, b=4) {
  K <- 2
  if (is.null(priors)) {
    priors <- array(1/K, dim=c(K))
  } else if (length(priors) != K) {
    stop(sprintf("You have specified %d priors for %d classes.", length(priors), K))
  } else if (sum(priors) != 1) {
    stop(sprintf("You have passed invalid priors. The sum(priors) should be 1; yours is %.3f", sum(priors)))
  }
  mu1 <- array(a, dim=c(d))
  mu1[2] <- b
  mus <- cbind(array(0, dim=c(d)), mu1)

  S <- diag(d)
  S[2,2] <- b

  S <- abind(diag(d), S, along=3)

  if (rotate) {
    res <- badmf.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }
  # simulate from GMM
  sim <- badmf.sims.sim_gmm(mus, S, n, priors)
  return(structure(list(X=sim$X, Y=as.factor(sim$Y), mus=mus, Sigmas=S, priors=sim$priors, simtype="Stacked Cigar",
                        params=c(a=a, b=b)), class="simulation"))
}


#' Xor Problem
#'
#' A function to simulate from the 2-class xor problem.
#' @param n the number of samples of the simulated data.
#' @param d the dimensionality of the simulated data.
#' @param priors the priors for each class. If \code{NULL}, class priors are all equal. If not null, should be \code{|priors| = K}, a length \code{K} vector for \code{K} classes. Defaults to \code{NULL}.
#' @param fall the falloff for the covariance structuring. Sigma declines by ndim/fall across the variance terms. Defaults to \code{100}.
#' @return A list of class \code{simulation} with the following:
#' \item{X}{\code{[n, d]} the \code{n} data points in \code{d} dimensions as a matrix.}
#' \item{Y}{\code{[n]} the \code{n} labels as an array.}
#' \item{mus}{\code{[d, K]} the \code{K} class means in \code{d} dimensions.}
#' \item{Sigmas}{\code{[d, d, K]} the \code{K} class covariance matrices in \code{d} dimensions.}
#' \item{priors}{\code{[K]} the priors for each of the \code{K} classes.}
#' \item{simtype}{The name of the simulation.}
#' \item{params}{Any extraneous parameters the simulation was created with.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("sims", package = "badmf")}
#'
#' @author Eric Bridgeford
#' @examples
#' library(badmf)
#' data <- badmf.sims.xor2(n=200, d=30)  # 200 examples of 30 dimensions
#' X <- data$X; Y <- data$Y
#' @export
badmf.sims.xor2 <- function(n, d, priors=NULL, fall=100) {
  K <- 2
  if (is.null(priors)) {
    priors <- array(1/K, dim=c(K))
  } else if (length(priors) != K) {
    stop(sprintf("You have specified %d priors for %d classes.", length(priors), K))
  } else if (sum(priors) != 1) {
    stop(sprintf("You have passed invalid priors. The sum(priors) should be 1; yours is %.3f", sum(priors)))
  }
  n1 <- ceiling(n/2)
  n2 <- floor(n/2)
  # first simulation set
  mus <- abind(array(0, dim=c(d)), array(c(1, 0), dim=c(d)), along=2)
  S <- sqrt(d/fall)*diag(d)
  S <- abind(S, S, along=3)

  # simulate from GMM for first set of training examples
  sim1 <- badmf.sims.sim_gmm(mus, S, n1, priors)

  # second simulation set
  mus <- abind(array(1, dim=c(d)), array(c(0, 1), dim=c(d)), along=2)

  # simulate from GMM for second set of training examples
  sim2 <- badmf.sims.sim_gmm(mus, S, n2, priors=priors)

  X <- abind(sim1$X, sim2$X, along=1)
  Y <- abind(sim1$Y, sim2$Y, along=1)

  reorder <- sample(n)
  return(structure(list(X=X[reorder,], Y=as.factor(Y[reorder]), mus=mus, Sigmas=S, priors=sim2$priors, simtype="Xor",
                        params=list(fall=fall)), class="simulation"))
}

#' Linear Simulation
#'
#' A function to simulate multi-class data with a linear class-mean trend. The signal dimension is the dimension carrying all of the
#' between-class difference, and the non-signal dimensions are noise.
#'
#' @import abind
#' @param n the number of samples.
#' @param d the number of dimensions. The first dimension will be the signal dimension; the remainders noise.
#' @param K the number of classes in the dataset.
#' @param signal.scale the scaling for the signal dimension. Defaults to \code{1}.
#' @param signal.lshift the location shift for the signal dimension between the classes. Defaults to \code{1}.
#' @param non.scale the scaling for the non-signal dimensions. Defaults to \code{1}.
#' @param class.equal whether the number of samples/class should be equal, with each
#' class having a prior of 1/K, or inequal, in which each class obtains a prior
#' of k/sum(K) for k=1:K. Defaults to \code{TRUE}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @author Eric Bridgeford
#' @export
badmf.sims.linear <- function(n, d, K, signal.scale=1, signal.lshift=1, non.scale=1, rotate=FALSE, class.equal=TRUE, ind=FALSE) {
  priors <- gen.sample.labels(K, class.equal=class.equal)
  S <- diag(d)
  S[1, 1] <- signal.scale
  S[-c(1), -c(1)] <- non.scale
  S <- abind(lapply(1:K, function(i) {
    S
  }), along=3)
  mu <- c(1, rep(0, d-1))
  mu[1] <- signal.lshift
  mus <- abind(lapply(1:K, function(i) {
    mu*i*signal.lshift
  }), along=2)

  if (rotate) {
    res <- badmf.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
    rotate=res$Q
  }

  sim <- badmf.sims.sim_gmm(mus, S, n, priors=priors)
  X <- sim$X; Y <- factor(sim$Y)
  return(list(X=X, Y=as.factor(Y), mus=mus, Sigmas=S, priors=priors, simtype="Linear",
              params=list(signal.scale=signal.scale, signal.lshift=signal.lshift, non.scale=non.scale,
                          rotate=rotate, class.equal=class.equal, ind=ind)))
}

#' Exponential Simulation
#'
#' A function to simulate multi-class data with an Exponential class-mean trend.
#'
#' @import abind
#' @param n the number of samples.
#' @param d the number of dimensions. The first dimension will be the signal dimension; the remainders noise.
#' @param K the number of classes in the dataset.
#' @param signal.scale the scaling for the signal dimension. Defaults to \code{1}.
#' @param signal.lshift the location shift for the signal dimension between the classes. Defaults to \code{1}.
#' @param non.scale the scaling for the non-signal dimensions. Defaults to \code{1}.
#' @param class.equal whether the number of samples/class should be equal, with each
#' class having a prior of 1/K, or inequal, in which each class obtains a prior
#' of k/sum(K) for k=1:K. Defaults to \code{TRUE}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @author Eric Bridgeford
#' @export
badmf.sims.exp <- function(n, d, K, signal.scale=1, signal.lshift=1, non.scale=1, rotate=FALSE, class.equal=TRUE, ind=FALSE) {
  priors <- gen.sample.labels(K, class.equal=class.equal)
  S <- diag(d)
  S[1, 1] <- signal.scale; S[-c(1), -c(1)] <- non.scale
  S <- abind(lapply(1:K, function(i) {
    S
  }), along=3)
  mu <- c(1, rep(0, d-1))
  mu[1] <- signal.lshift
  mus <- abind(lapply(1:K, function(i) {
    mu*exp(i)*signal.lshift
  }), along=2)

  if (rotate) {
    res <- badmf.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
    rotate=res$Q
  }

  sim <- badmf.sims.sim_gmm(mus, S, n, priors=priors)
  X <- sim$X; Y <- factor(sim$Y)
  return(list(X=X, Y=as.factor(Y), mus=mus, Sigmas=S, priors=priors, simtype="Exp",
              params=list(signal.scale=signal.scale, signal.lshift=signal.lshift,
                          non.scale=non.scale, rotate=rotate, class.equal=class.equal,
                          ind=ind)))
}

#' Radial Simulation
#'
#' A function to simulate data with the same mean with radial symmetry as class id increases.
#'
#' @import abind
#' @param n the number of samples.
#' @param d the number of dimensions.
#' @param K the number of classes in the dataset.
#' @param er.scale the scaling for the error of the samples. Defaults to \code{0.1}.
#' @param r the radial spacing between each class. Defaults to \code{1}.
#' @param class.equal whether the number of samples/class should be equal, with each
#' class having a prior of 1/K, or inequal, in which each class obtains a prior
#' of k/sum(K) for k=1:K. Defaults to \code{TRUE}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @examples
#' library(badmf)
#' sim <- badmf.sims.radial(100, 3, 2)
#' @author Eric Bridgeford
#' @export
badmf.sims.radial <- function(n, d, K, er.scale=0.1, r=1, class.equal=TRUE, ind=FALSE) {
  priors <- gen.sample.labels(K, class.equal=class.equal)

  class <- 1:K
  # sample n points from 1:K for assignment of class label
  Y <- sample(class, size=n, replace=TRUE, prob=priors)
  # initialize X
  X <- array(NaN, dim=c(n, d))

  # loop over class labels that we are sampling from
  for (i in class[class %in% Y]) {
    if (i == 1) {
      f <- badmf.sims.2ball
    } else {
      f <- badmf.sims.2sphere
    }
    repvec <- Y == i
    pts <- do.call(f, list(sum(repvec), d, r=r*i, cov.scale=er.scale))
    X[repvec,] <- pts
  }

  if (ind) {
    Y <- sample(class, size=n, replace=TRUE, prob=priors)
  }
  return(list(X=X, Y=factor(Y), priors=priors, simtype="Radial",
              params=list(er.scale=er.scale, r=r, class.equal=class.equal,
                          ind=ind)))
}

#' Beta Simulation
#'
#' Signal Dimension is beta-distributed; remainder are Uniform 0, 1
#'
#' @import abind
#' @param n the number of samples.
#' @param d the number of dimensions.
#' @param class.equal whether the number of samples/class should be equal, with each
#' class having a prior of 1/K, or inequal, in which each class obtains a prior
#' of k/sum(K) for k=1:K. Defaults to \code{TRUE}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @author Eric Bridgeford
#' @export
badmf.sims.beta <- function(n, d, class.equal=TRUE, ind=FALSE) {
  K <- 5
  priors <- gen.sample.labels(K, class.equal=class.equal)
  class <- 1:K
  # sample n points from 1:K for assignment of class label
  Y <- sample(class, size=n, replace=TRUE, prob=priors)
  # initialize X
  X <- array(NaN, dim=c(n, d))
  alpha <-c(5, 5, 5, 2, 1)
  beta <- c(1, 2, 5, 5, 5)
  for (i in 1:K) {
    # fill the top dimension with beta distributed data
    X[Y == i,1] <- rbeta(sum(Y == i), alpha[i], beta[i])
  }
  if (d > 1) {
    X[, 2:d] <- runif(n*(d - 1))  # fill the rest with uniforms
  }
  if (ind) {
    Y <- sample(class, size=n, replace=TRUE, prob=priors)
  }
  return(list(X=X, Y=factor(Y), priors=priors, simtype="Beta",
              params=list(alpha=alpha, beta=beta)))
}


#'
#' A helper function to generate a d-dimensional linear transformation matrix.
#' @param d the number of dimensions.
#' @return A \code{[d]} the coefficient vector.
#' @author Eric Bridgeford
gen.coefs <- function(d) {
  A = as.array(1/1:d, dim=c(d, 1))
  return(A)
}

#' A helper function to generate n samples of a d-dimensional uniform vector.
#' @param n the number of samples.
#' @param d the number of dimensions.
#' @param a the lower limit.
#' @param b the upper limit.
#' @param x \code{[n, d]} the simulated data matrix.
#' @author Eric Bridgeford
#' @importFrom stats runif
gen.x.unif <- function(n, d, a=-1, b=1) {
  x <- array(runif(n=(n*d), min=a, max=b), dim=c(n, d))
  return(x)
}

#' GMM Simulate
#'
#' A helper function for simulating from Gaussian Mixture.
#' @param mus \code{[d, K]} the mus for each class.
#' @param Sigmas \code{[d,d,K]} the Sigmas for each class.
#' @param n the number of examples.
#' @param priors \code{K} the priors for each class.
#' @return A list with the following:
#' \item{X}{\code{[n, d]} the simulated data.}
#' \item{Y}{\code{[n]} the labels for each data point.}
#' \item{priors}{\code{[K]} the priors for each class.}
#' @author Eric Bridgeford
#' @importFrom MASS mvrnorm
badmf.sims.sim_gmm <- function(mus, Sigmas, n, priors) {
  K <- dim(mus)[2]
  labs <- sample(1:K, size=n, prob=priors, replace=TRUE)
  ylabs <- as.vector(sort(unique(labs)))
  res <- sapply(ylabs, function(y) mvrnorm(n=sum(labs == y), mus[,y], Sigmas[,,y]), USE.NAMES=TRUE, simplify=FALSE)
  X <- array(0, dim=c(n, dim(Sigmas)[1]))
  for (y in ylabs) {
    X[labs == y,] <- res[[y]]
  }
  return(list(X=X, Y=labs, priors=priors))
}

#' Sample Random Rotation
#'
#' A helper function for estimating a random rotation matrix.
#' @importFrom stats rnorm
#' @param d dimensions to generate a rotation matrix for.
#' @return the rotation matrix
#' @author Eric Bridgeford
badmf.sims.rotation <- function(d) {
  Q <- qr.Q(qr(array(rnorm(d*d), dim=c(d, d))))
  if (det(Q) < -.99) {
    Q[,1] <- -Q[,1]
  }
  return(Q)
}

#' Random Rotation
#'
#' A helper function for applying a random rotation to gaussian parameter set.
#' @param mus means per class.
#' @param Sigmas covariances per class.
#' @param Q rotation to use, if any
#' @author Eric Bridgeford
badmf.sims.random_rotate <- function(mus, Sigmas, Q=NULL) {
  dimm <- dim(mus)
  K <- dimm[2]
  d <- dim(mus)[1]
  if (is.null(Q)) {
    Q <- badmf.sims.rotation(d)
  } else if (!isTRUE(all.equal(dim(Q), c(d, d)))) {
    stop(sprintf("You have specified a rotation matrix with dimensions (%d, %d), but should be (%d, %d).", dim(Q)[1], dim(Q)[2], d, d))
  }

  for (i in 1:K) {
    mus[,i] <- Q %*% mus[,i,drop=FALSE]
    Sigmas[,,i] <- Q %*% Sigmas[,,i] %*% t(Q)
  }
  return(list(mus=mus, S=Sigmas, Q=Q))
}

#' Sample from Unit 2-Ball
#'
#' Sample from the 2-ball in d-dimensions.
#'
#' @importFrom MASS mvrnorm
#' @param n the number of samples.
#' @param d the number of dimensions.
#' @param r the radius of the 2-ball. Defaults to \code{1}.
#' @param cov.scale if desired, sample from 2-ball with error sigma. Defaults to \code{NaN},
#' which has no noise.
#' @examples
#' library(badmf)
#' # sample 100 points from 3-d 2-ball with radius 2
#' X <- badmf.sims.rball(100, 3, 2)
#' @author Eric Bridgeford
#' @export
badmf.sims.2ball <- function(n, d, r=1, cov.scale=0) {
  Y <- mvrnorm(n=n, mu=array(0, dim=c(d, 1)), Sigma=diag(d))
  u <- runif(n)
  r <- r * u^(1/d)
  X <- r * Y/sqrt(apply(Y^2, 1, sum)) + mvrnorm(n=n, mu=array(0, dim=c(d,1)), Sigma=cov.scale*diag(d))
}

#' Sample from Unit 2-Sphere
#'
#' Sample from the 2-sphere in d-dimensions.
#'
#' @importFrom MASS mvrnorm ginv
#' @param n the number of samples.
#' @param d the number of dimensions.
#' @param r the radius of the 2-ball. Defaults to \code{1}.
#' @param cov.scale if desired, sample from 2-ball with error sigma. Defaults to \code{0},
#' which has no noise.
#' @examples
#' library(badmf)
#' # sample 100 points from 3-d 2-ball with radius 2
#' X <- badmf.sims.rball(100, 3, 2)
#' @author Eric Bridgeford
#' @export
badmf.sims.2sphere <- function(n, r, d, cov.scale=0) {
  u <- mvrnorm(n=n, mu=array(0, dim=c(d,1)), Sigma=diag(d))
  unorm <- diag(sqrt(apply(u^2, 1, sum)))
  pts <- r*(ginv(unorm) %*% u)
  pts <- pts + mvrnorm(n=n, mu=array(0, dim=c(d,1)), Sigma=cov.scale*diag(d))
}
