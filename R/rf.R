#' Decision Tree Frontend Interface
#'
#' @param form a valid formula object.
#' @param data the data associated with the formula.
#' @param family the family to use for the feature prior. If NULL, does not construct a prior.
#' @return A trained decision tree.
#' @author Eric Bridgeford
#' @export
rf.fit <- function(formula, data, family="dirichlet", method="forest.fit", prior=NULL) {
  mf <- Call <- match.call()
  print(mf)
  m <- match(c("formula", "data", "subset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(mf)
  Terms <- attr(mf, "terms")

  Y <- model.response(mf)
  X <- if (!is.empty.model(Terms))
    model.matrix(Terms, mf, contrasts)
  fit <- do.call(forest.fit, c(X, Y, family, prior))
  return(fit)
}

forest.fit <- function(X, Y, family, prior) {

}
