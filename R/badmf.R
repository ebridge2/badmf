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
    fit <- do.call(badmf.class.fit, list(X[,-c(1)], Y, d, alpha, depth.max,
                                         size, debug, no.cores, train.params))
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
#' \item{`prior`}{the hyperparameters of the Dirichlet Prior.}
#' \item{`post`}{the hyperparamaters of the Dirichlet Posterior.}
#' @author Eric Bridgeford
#' @export
badmf.class.fit <- function(X, Y, d=NULL, alpha=NULL, ntrees=10L, depth.max=1L,
                            size=1L, debug=FALSE, no.cores=1L, train.params=NULL) {
  Y <- as.factor(Y)
  n <- length(Y); p <- dim(X)[2]

  if (is.null(d)) {
    d <- as.integer(round(sqrt(p)))
  }

  # if alpha is null, sample uniformly with Dir(1, 1, ...)
  if (is.null(alpha)) {
    alpha <- rep(1, p)
  } else {
    if (!ifelse(is.numeric(alpha), all(alpha > 0), FALSE)) {
      stop("You have not entered a valid Dirichlet prior. All values should be > 0.")
    }
  }
  rf.input.validator(d=d, ntrees=ntrees, bagg=bagg, depth.max=depth.max, size=size)
  if (is.null(train.params)) {
    train.params <- list()
  }
  # initialize training vector
  train.params <- list(d=ifelse(!is.null(train.params$d), train.params$d, d),
                       ntrees=ifelse(!is.null(train.params$ntrees), train.params$ntrees, ntrees),
                       bagg=ifelse(!is.null(train.params$bagg), train.params$bagg, bagg),
                       depth.max=ifelse(!is.null(train.params$depth.max), train.params$depth.max, depth.max),
                       size=ifelse(!is.null(train.params$size), train.params$size, size))
  tryCatch({
      do.call(rf.input.validator, c(train.params))
  }, error=function(e) {
    print("Your training parameters are invalid. Revise train.params.")
    stop(e)
  })

  # train random forest classifier with specified prior
  fit.rf <- rf.class.fit(X, Y, train.params$d, alpha, train.params$ntrees, train.params$bagg,
                         train.params$depth.max, train.params$size, no.cores=no.cores)

  # construct dirichlet posterior
   alpha.post <- alpha + count.forest(fit.rf, p)

  # train random forest classifier with posterior
   fit.rf <- rf.class.fit(X, Y, d, alpha, ntrees, bagg, depth.max, size, no.cores=no.cores)

   fit.rf$method <- "badmf.class.fit"; fit.rf$alpha <- NULL; fit.rf$prior <- alpha; fit.rf$post <- alpha.post
   return(fit.rf)
}

count.forest <- function(rf.fit, p) {
  rf.ct <- do.call(c, lapply(rf.fit$forest, function(tree) {
    count.tree(tree$tree)
  })) %>%
    table
  counts <- rep(0, p)
  counts[as.numeric(names(rf.ct))] <- as.numeric(rf.ct)
  return(counts)
}
count.tree <- function(tree) {
  if (class(tree) == "leaf.node") {
    return(NULL)
  }
  return(c(tree$feature, count.tree(tree$left), count.tree(tree$right)))
}
