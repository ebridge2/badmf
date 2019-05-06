#' Decision Tree Frontend Interface
#'
#' @param form a valid formula object.
#' @param data the data associated with the formula.
#' @param family the family to use for the feature prior. If NULL, does not construct a prior.
#' @return A trained decision tree.
#' @author Eric Bridgeford
#' @export
dec_tree.fit <- function(formula, data, family="dirichlet", method="tree.fit", prior=NULL) {
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
  fit <- do.call(tree.fit, c(X, Y))
  return(fit)
}

tree.fit <- function(X, Y) {

}

#' Create a new split of the data.
#' @param X the predictors.
#' @param Y the responses.
#' @param i the feature index.
#' @param t the threshold for the feature.
#' @return the split, as a list.
#' @author Eric Bridgeford
create.split <- function(X, Y, i, t) {
  idx <- X[,i] < t
  split <- list(left=list(X=X[idx,], Y=Y[idx]),
                right=list(X=X[!idx,], Y=Y[!idx]))
  return(split)
}

#' A function for computing the importance index for a split.
#' @param groups a 2 element list, where the first element are the samples to the left
#' and the second element the samples to the right of a split.
#' @return the gini impurity index.
#' @author Eric Bridgeford
importance.idx <- function(groups) {
  groups <- lapply(groups, function(group) return(group$Y))
  # the total number of observations across groups
  N <- sum(sapply(groups, function(group) length(group)))
  sum(sapply(groups, function(group) {
    # the total number of observations within group
    n <- length(group)
    # if group is empty, do not factor into gini
    # computation
    if (n == 0) {
      return(0)
    }
    # the unique ys
    ys <- levels(group)
    # compute the squared proportion of each y in this particular group
    sc <- sum(sapply(ys, function(y) {
      mean(group == y)^2
    }))
    return((1 - sc)*n/N)
  }))
}
