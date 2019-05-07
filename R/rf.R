#' Random Forest (RF)
#'
#' Fit a Random Forest with a `stats`-like formula frontend interface.
#'
#' @param formuler ravioli ravioli give me the formuoli.
#' @param data the data associated with the formuler. Note: if you want an intercept, you must
#' add it ahead of time.
#' @param d the number of features to subsample at a split node. Defaults to \code{as.integer(round(sqrt(p)))}, where \code{p}
#' is the dimensionality of the predictor feature space.
#' @param alpha the feature sampling prior. Should be a \code{[p]} vector, where \code{p} is the number of predictors.
#' Corresponds to alpha for a Dirichlet distribution.   If \code{NULL}, samples uniformly.
#' @param ntrees the number of trees to construct. Defaults to \code{10L}.
#' @param bagg the relative size of the subsamples for the training set. A numeric s.t.
#' \code{0 < bagg <= 1}. Each subsample will be \code{bagg*nsamples} elements. Defaults to \code{0.632}.
#' @param method whether you want "classification" or "regression".
#' @param depth.max the maximum allowed tree depth.
#' @param size the minimum allowed number of samples for an individual node.
#' @param debug whether to save the predictors and responses that are categorized
#' @param mc.cores the number of cores to use. Should be \code{0 < mc.cores <= parallel::detectCores()}.
#' Any unset parameters will default to the values provided above (or the corresponding defaults if unprovided).
#' @return an object of class \code{rf.class} containing the following:
#' \item{\code{forest}}{A list a decision trees.}
#' \item{\code{method}}{the method used to fit the forest.}
#' @author Eric Bridgeford
#' @export
rf.fit <- function(formuler, data=NULL, d=NULL, alpha=NULL, ntrees=10L, bagg=0.632, method="classification",
                   depth.max=5L, size=1L, debug=FALSE, mc.cores=1L, ...) {
  call <- match.call()

  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm))
      names(Y) <- nm
  }
  X <- if (!is.empty.model(mt)) {
    model.matrix(mt, mf, contrasts)
  }

  if (method == "classification") {
    fit <- do.call(badmf.class.fit, list(X[,-c(1)], Y, d, depth.max, size, debug, mc.cores))
    fit$formula <- formuler
  } else {
    stop("Not yet implemented!")
  }
  return(fit)
}

#' Random Forest (RF) for Classification
#'
#' Fit a Random Forest Classifier.
#'
#' @param X the predictors. A \code{[n, p]} matrix.
#' @param Y the responses. A \code{[n]} vector or, optionally, a factor.
#' @param d the number of features to subsample at a split node. Defaults to \code{as.integer(round(sqrt(p)))}, where \code{p}
#' is the dimensionality of the predictor feature space.
#' @param alpha the feature sampling prior. Corresponds to alpha for a Dirichlet distribution. If \code{NULL}, samples uniformly.
#' @param ntrees the number of trees to construct. Defaults to \code{10L}.
#' @param bagg the relative size of the subsamples for the training set. A numeric s.t.
#' \code{0 < bagg <= 1}. Each subsample will be \code{bagg*n} elements. Defaults to \code{0.632}.
#' @param depth.max the maximum allowed tree depth. Defaults to \code{5L}.
#' @param size the minimum allowed number of samples for an individual node. Defaults to \code{1L}.
#' @param debug whether to save the predictors and responses that are categorized. Defaults to \code{FALSE}.
#' @param mc.cores the number of cores to use. Should be \code{0 < mc.cores <= parallel::detectCores()}.
#' Defaults to \code{1L}.
#'
#' @return an object of class \code{rf.class} containing the following:
#' \item{\code{forest}}{A list a decision trees.}
#' \item{\code{method}}{the method used to fit the forest.}
#' \item{\code{alpha}}{the hyperparams for sampling distn of feature probabilities.}
#' @author Eric Bridgeford
#' @importFrom parallel mclapply detectCores
#' @export
rf.class.fit <- function(X, Y, d=NULL, alpha=NULL, ntrees=10L,  bagg=0.632, depth.max=5L,
                         size=1L, debug=FALSE, mc.cores=1L, ...) {
  Y <- as.factor(Y)
  n <- length(Y); p <- dim(X)[2]

  if (is.null(d)) {
    d <- as.integer(round(sqrt(p)))
  }

  if (!is.integer(mc.cores) || mc.cores > detectCores()) {
    stop("You have passed an invalid entry for mc.cores. Should be a positive integer <= detectCores().")
  }

  # if alpha is null, sample uniformly with Dir(1, 1, ...)
  if (!ifelse(is.numeric(alpha), all(alpha > 0), is.null(alpha))) {
    stop("You have not entered a valid Dirichlet prior. Should be NULL, or a p vector with all values should be > 0.")
  }

  rf.input.validator(d=d, p=p, ntrees=ntrees, bagg=bagg, depth.max=depth.max, size=size)

  fit <- structure(list(forest=mclapply(1:ntrees, function(i) {
    ss <- sample(1:n, round(bagg*n))
    Xs <- X[ss,,drop=FALSE]; Ys <- Y[ss]
    dec.tree.class.fit(Xs, Ys, d, alpha, depth.max, size, debug)
  }, mc.cores=mc.cores), method="rf.class.fit", alpha=alpha), class="rf.class")
  return(fit)
}

#' Input Validator
#' @param d the number of features to subsample at each node. Defaults to \code{sqrt(p)}.
#' @param p the number of features total.
#' @param ntrees the number of trees to construct. Defaults to \code{10}.
#' @param bagg the relative size of the subsamples for the training set. Defaults to \code{0.632}. A numeric s.t.
#' \code{0 < bagg <= 1}. Each subsample will be \code{bagg*n} elements.
#' @param depth.max the maximum allowed tree depth.
#' @param size the minimum allowed number of samples for an individual node.
#' @return nothing; throws an error if the input is invalid.
#' @author Eric Bridgeford
rf.input.validator <- function(d, p, ntrees, bagg, depth.max, depth, size) {
  if (!ifelse(is.integer(d), d <= p & d > 0, FALSE)) {
    stop("d should be a positive integer <= p, or NULL to indicate to sample every feature.")
  }

  if (!is.integer(ntrees) || ntrees < 0) {
    stop("You have passed an invalid entry for ntrees. Should be a positive integer.")
  }

  if (!ifelse(is.numeric(bagg), bagg <= 1 & bagg > 0, FALSE)) {
    stop("You have not entered a valid value for bagg. Should be a numeric 0 < bagg <= 1.")
  }

  if (!ifelse(is.integer(depth.max), depth.max > 0, FALSE)) {
    stop("You have not entered a valid value for depth.max. Should be a positive integer.")
  }

  if (!ifelse(is.integer(size), size > 0, FALSE)) {
    stop("You have not passed a valid value for size. Should be a positive integer.")
  }
}

#' Count Feature Utilization
#'
#' A function to count the number of times a feature is used
#' in a forest at split nodes.
#'
#' @param object a fit random forest of class \code{rf.class}.
#' @param ... trailing args.
#' @return the counts of all features in the forest.
#' @author Eric Bridgeford
#' @export
count.features.rf.class <- function(object, ...) {
  rf.ct <- Reduce("+", lapply(object$forest, function(tree) {
    count.features(tree)
  }))
  return(rf.ct)
}
