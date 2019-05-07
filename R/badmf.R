#' Bayesian Decision-Making Forest (BaD-MF)
#'
#' Fit a Bayesian Decision-Making Forest with a `stats`-like formula frontend interface.
#'
#' @param formuler ravioli ravioli give me the formuoli.
#' @param data the data associated with the formuler. Note: if you want an intercept, you must
#' add it ahead of time.
#' @param d the number of features to subsample at a split node. Defaults to sqrt(nsamples).
#' @param alpha the feature sampling prior. Should be a `[p]` vector, where `p` is the number of predictors.
#' Corresponds to alpha for a Dirichlet distribution. If `NULL`, samples uniformly for the initial
#' training iteration.
#' @param ntrees the number of trees to construct. Defaults to 10.
#' @param bagg the relative size of the subsamples for the training set. Defaults to 0.632. A numeric s.t.
#' `0 < bagg <= 1`. Each subsample will be `bagg*nsamples`` elements.
#' @param method whether you want "classification" or "regression".
#' @param depth the maximum allowed tree depth.
#' @param size the minimum allowed number of samples for an individual node.
#' @param debug whether to save the predictors and responses that are categorized
#' @param no.cores the number of cores to use. Should be `0 < no.cores <= parallel::detectCores()`.
#' @param train.params if you wish to provide specialized parameters for training, a list containing the following:
#' \itemize{
#' \item{`d`}{the number of features to subsample at a split node.}
#' \item{`ntrees`}{the number of trees to construct.}
#' \item{`bagg`}{the relative size of the subsamples from the training set.}
#' \item{`depth`}{the maximum allowed tree depth.}
#' \item{`size`}{the minimum allowed number of samples for an individual node.}
#' }
#' Any unset parameters will default to the values provided above (or the corresponding defaults if unprovided).
#' @return an object of class `rf` containing the following:
#' \item{`forest`}{A list a decision trees.}
#' \item{`method`}{the method used to fit the forest.}
#' @author Eric Bridgeford
#' @export
badmf.fit <- function(formuler, data=NULL, d=NULL, alpha=NULL, ntrees=10L, bagg=0.632, method="classification",
                      depth.max=1L, size=1L, debug=FALSE, no.cores=1L, train.params=NULL) {
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
    fit <- do.call(badmf.class.fit, list(X[,-c(1)], Y, d, alpha, depth.max, size, debug, no.cores, train.params))
    fit$formula <- formuler
  } else {
    stop("Not yet implemented!")
  }
  return(fit)
}

#' Bayesian Decision-Making Forest (BaD-MF) for Classification
#'
#' Fit a Bayesian Decision-Making Forest Classifier.
#'
#' @param X the predictors. A `[n, p]` matrix.
#' @param Y the responses. A `[n]` vector or, optionally, a factor.
#' @param d the number of features to subsample at each node. Defaults to `sqrt(p)`.
#' @param alpha the feature sampling prior. Corresponds to alpha for a Dirichlet distribution. If `NULL`, samples uniformly for the initial training iteration.
#' @param ntrees the number of trees to construct. Defaults to 10.
#' @param bagg the relative size of the subsamples for the training set. Defaults to 0.632. A numeric s.t.
#' `0 < bagg <= 1`. Each subsample will be `bagg*n`` elements.
#' @param depth the maximum allowed tree depth.
#' @param size the minimum allowed number of samples for an individual node.
#' @param debug whether to save the predictors and responses that are categorized
#' @param no.cores the number of cores to use. Should be `0 < no.cores <= parallel::detectCores()`.
#' @param train.params if you wish to provide specialized parameters for training, a list containing the following:
#' \itemize{
#' \item{`d`}{the number of features to subsample at a split node.}
#' \item{`ntrees`}{the number of trees to construct.}
#' \item{`bagg`}{the relative size of the subsamples from the training set.}
#' \item{`depth`}{the maximum allowed tree depth.}
#' \item{`size`}{the minimum allowed number of samples for an individual node.}
#' }
#' Any unset parameters will default to the values provided above (or the corresponding defaults if unprovided).
#' @return an object of class `rf` containing the following:
#' \item{`forest`}{A list a decision trees.}
#' \item{`method`}{the method used to fit the forest.}
#' @author Eric Bridgeford
#' @importFrom parallel detectCores, mclapply
badmf.class.fit <- function(X, Y, d=NULL, alpha=NULL, ntrees=10L, depth.max=1L,
                            size=1L, debug=FALSE, no.cores=1L, train.params=NULL) {
  Y <- as.factor(Y)
  n <- length(Y); p <- dim(X)[2]

  if (is.null(d)) {
    d <- round(sqrt(p))
  }

  # if alpha is null, sample uniformly with Dir(1, 1, ...)
  if (is.null(alpha)) {
    alpha <- rep(1, p)
  } else {
    if (!ifelse(is.numeric(alpha), all(alpha > 0), FALSE)) {
      stop("You have not entered a valid Dirichlet prior. All values should be > 0.")
    }
  }
  rf.input.validator(d=d, ntrees=ntrees, bagg=bagg, depth=depth, size=size)
  if (is.null(train.params)) {
    train.params <- list()
  }
  # initialize training vector
  train.params <- list(d=ifelse(!is.null(train.params$d), train.params$d, d),
                       ntrees=ifelse(!is.null(train.params$ntrees), train.params$ntrees, ntrees),
                       bagg=ifelse(!is.null(train.params$bagg), train.params$bagg, bagg),
                       depth=ifelse(!is.null(train.params$depth), train.params$depth, depth),
                       size=ifelse(!is.null(train.params$size), train.params$size, size))
  tryCatch({
      do.call(rf.input.validator, c(train.params, list(no.cores=no.cores)))
  }, error=function(e) {
    print("Your training parameters are invalid.")
    stop(e)
  })

  # train random forest classifier with specified prior
  fit.rf <- rf.class.fit(X, Y, d, alpha, ntrees, bagg, depth, size, no.cores=no.cores)

  # construct dirichlet posterior
}
