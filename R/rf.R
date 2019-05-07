#' Random Forest (RF)
#'
#' Fit a Random Forest with a `stats`-like formula frontend interface.
#'
#' @param formuler ravioli ravioli give me the formuoli.
#' @param data the data associated with the formuler. Note: if you want an intercept, you must
#' add it ahead of time.
#' @param d the number of features to subsample at a split node. Defaults to sqrt(nsamples).
#' @param alpha the feature sampling prior. Should be a `[p]` vector, where `p` is the number of predictors.
#' Corresponds to alpha for a Dirichlet distribution.   If `NULL`, samples uniformly.
#' @param ntrees the number of trees to construct. Defaults to 10.
#' @param bagg the relative size of the subsamples for the training set. Defaults to 0.632. A numeric s.t.
#' `0 < bagg <= 1`. Each subsample will be `bagg*nsamples`` elements.
#' @param method whether you want "classification" or "regression".
#' @param depth.max the maximum allowed tree depth.
#' @param size the minimum allowed number of samples for an individual node.
#' @param debug whether to save the predictors and responses that are categorized
#' @param no.cores the number of cores to use. Should be `0 < no.cores <= parallel::detectCores()`.
#' Any unset parameters will default to the values provided above (or the corresponding defaults if unprovided).
#' @return an object of class `rf` containing the following:
#' \item{`forest`}{A list a decision trees.}
#' \item{`method`}{the method used to fit the forest.}
#' @author Eric Bridgeford
#' @export
rf.fit <- function(formuler, data=NULL, d=NULL, alpha=NULL, ntrees=10L, bagg=0.632, method="classification",
                   depth.max=1L, size=1L, debug=FALSE, no.cores=1L) {
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
    fit <- do.call(badmf.class.fit, list(X[,-c(1)], Y, d, depth.max, size, debug, no.cores))
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
#' @param X the predictors. A `[n, p]` matrix.
#' @param Y the responses. A `[n]` vector or, optionally, a factor.
#' @param d the number of features to subsample at each node. Defaults to `sqrt(p)`.
#' @param alpha the feature sampling prior. Corresponds to alpha for a Dirichlet distribution. If `NULL`, samples uniformly.
#' @param ntrees the number of trees to construct. Defaults to 10.
#' @param bagg the relative size of the subsamples for the training set. Defaults to 0.632. A numeric s.t.
#' `0 < bagg <= 1`. Each subsample will be `bagg*n`` elements.
#' @param depth.max the maximum allowed tree depth.
#' @param size the minimum allowed number of samples for an individual node.
#' @param debug whether to save the predictors and responses that are categorized
#' @param no.cores the number of cores to use. Should be `0 < no.cores <= parallel::detectCores()`.
#' @param train.params if you wish to provide specialized parameters for training, a list containing the following:
#' \itemize{
#' \item{`d`}{the number of features to subsample at a split node.}
#' \item{`ntrees`}{the number of trees to construct.}
#' \item{`bagg`}{the relative size of the subsamples from the training set.}
#' \item{`depth.max`}{the maximum allowed tree depth.}
#' \item{`size`}{the minimum allowed number of samples for an individual node.}
#' }
#' Any unset parameters will default to the values provided above (or the corresponding defaults if unprovided).
#' @return an object of class `rf` containing the following:
#' \item{`forest`}{A list a decision trees.}
#' \item{`method`}{the method used to fit the forest.}
#' \item{`alpha`}{the hyperparams for sampling distn of feature probabilities.}
#' @author Eric Bridgeford
#' @importFrom parallel mclapply
#' @export
rf.class.fit <- function(X, Y, d=NULL, alpha=NULL, ntrees=10L,  bagg=0.632, depth.max=1L,
                         size=1L, debug=FALSE, no.cores=1L) {
  Y <- as.factor(Y)
  n <- length(Y); p <- dim(X)[2]

  if (is.null(d)) {
    d <- as.integer(round(sqrt(p)))
  }

  if (!is.integer(no.cores) || no.cores > detectCores()) {
    stop("You have passed an invalid entry for no.cores.")
  }

  # if alpha is null, sample uniformly with Dir(1, 1, ...)
  if (!ifelse(is.numeric(alpha), all(alpha > 0), is.null(alpha))) {
    stop("You have not entered a valid Dirichlet prior. All values should be > 0.")
  }
  rf.input.validator(d=d, ntrees=ntrees, bagg=bagg, depth.max=depth.max, size=size)

  fit <- structure(list(forest=mclapply(1:ntrees, function(i) {
    ss <- sample(1:n, round(bagg*n))
    Xs <- X[ss,,drop=FALSE]; Ys <- Y[ss]
    print(dim(Xs)); print(length(Ys)); print(Ys)
    tree.class.fit(Xs, Ys, d, alpha, depth.max, size, debug)
  }, mc.cores=no.cores), method="rf.class.fit", alpha=alpha), class="rf.class")
  return(fit)
}

#' Input Validator
#' @param d the number of features to subsample at each node. Defaults to `sqrt(p)`.
#' @param ntrees the number of trees to construct. Defaults to 10.
#' @param bagg the relative size of the subsamples for the training set. Defaults to 0.632. A numeric s.t.
#' `0 < bagg <= 1`. Each subsample will be `bagg*n`` elements.
#' @param depth.max the maximum allowed tree depth.
#' @param size the minimum allowed number of samples for an individual node.
#' @importFrom parallel detectCores
rf.input.validator <- function(d, ntrees, bagg, depth.max, depth, size) {
  if (!ifelse(is.integer(d), d <= p & d > 0, FALSE)) {
    stop("d should be a positive integer <= p, or NULL to indicate to sample every feature.")
  }

  if (!is.integer(ntrees) || ntrees < 0) {
    stop("You have passed an invalid entry for ntrees.")
  }

  if (!ifelse(is.numeric(bagg), bagg <= 1 & bagg > 0, FALSE)) {
    stop("You have not entered a valid value for bagg.")
  }

  if (!ifelse(is.integer(depth.max), depth.max > 0, FALSE)) {
    stop("You have not entered a valid value for depth.max.")
  }

  if (!ifelse(is.integer(size), size > 0, FALSE)) {
    stop("You have not passed a valid value for size.")
  }
}
