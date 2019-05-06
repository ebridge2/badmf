#' Bayesian Decision Tree Fit
#'
#' Fit a Bayesian Decision Tree.
#'
#' @param formuler ravioli ravioli give me the formuoli
#' @param data the data associated with the formuler.
#' @param family the family to use for the feature prior. If NULL, does not construct a prior.
#' @return A trained decision tree.
#' @author Eric Bridgeford
#' @export
dec_tree.fit <- function(formuler, data, family="dirichlet", method="tree.fit", prior=NULL) {
  mf <- Call <- match.call()
  print(mf)
  m <- match(c("formuler", "data", "subset"), names(mf), 0)
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

#' build tree recursively
#' @param split a particular split node.
#' @param d the number of features to subsample.
#' @param depth.max the max tree depth.
#' @param size the minimum number of elements at a particular.
#' @importFrom forcats fct_c
build.tree <- function(split, d, depth.max, size, depth) {
  # if we have an empty child, create a leaf by merging
  # so that we don't have a totally empty child
  if (length(split$left$Y) == 0 || length(split$right$Y) == 0) {
    split$left <- split$right <- leaf.node(fct_c(split$left$Y,
                                                 split$right$Y))
    return(split)
  }

  # if we are at the max depth, create a leaf for the left and
  # right children
  if (depth >= depth.max) {
    split$left <- leaf.node(split$left)
    split$right <- leaf.node(split$right)
    return(split)
  }
  # process left child
  if (length(split$left$Y) > size) {
    # split the node if we can still do better
    split$right <- get.split(split$right$X, split$right$Y, d)
    split <- build.tree(split$right, d, depth.max, size, depth + 1)
    return(split)
  } else {
    split$right <- leaf.node(split$right$Y)
    return(split)
  }

  # process left child
  if (length(split$left$Y) > size) {
    # split the node if we can still do better
    split$left <- get.split(split$left$X, split$left$Y, d)
    split <- build.tree(split$left, d, depth.max, size, depth + 1)
    return(split)
  } else {
    split$left <- leaf.node(split$left$Y)
    return(split)
  }
  return(split)
}

#' Fit a tree.
#' @param X the predictors.
#' @param Y the responses.
#' @param d the number of features to subsample.
#' @param depth the maximum allowed tree depth.
#' @param size the minimum allowed number of samples for an individual node.
tree.fit <- function(X, Y, d, depth.max=1, size=1) {
  n <- length(Y); p <- dim(X)[2]
  if (!(ifelse(is.integer(p), d <= p & d > 0, FALSE) ||
      (ifelse(is.null(p), TRUE, FALSE)))) {
    stop("d should be a positive integer <= p, or NULL to indicate to sample every feature.")
  }
  tree <- build.tree(get.split(X, Y, d), d, depth.max, size, 1)
  return(tree)
}

#' Create a new split of the data.
#' @param X the predictors.
#' @param Y the responses.
#' @param i the feature index.
#' @param t the threshold for the feature.
#' @return the split, as a list.
#' @author Eric Bridgeford
create.split <- function(X, Y, i, t) {
  # split based on the feature of interest and threshold
  # of interest
  idx <- X[,i] < t
  # return split as a list
  split <- list(left=list(X=X[idx,], Y=Y[idx]),
                right=list(X=X[!idx,], Y=Y[!idx]),
                feature=i, threshold=t)
  return(split)
}

#' Find Best split in the data.
#'
#' @param X the predictors.
#' @param Y the responses.
#' @param d the number of features to subsample. Should be an integer 0 < d <= p.
#' @return the best split.
get.split <- function(X, Y, d) {
  n <- length(Y); p <- dim(X)[2]
  # sample features to check
  features <- sample(1:p, d)
  # initialize best split
  best_split <- list(t=NULL, feature=NULL, score=1.1, split=NULL)
  # loop over subsampled features
  for (feature in features) {
    # check all possible splits for feature of interest
    for (t in X[,feature]) {
      # split on feature, t
      sp <- create.split(X, Y, feature, t)
      # get importance score
      sp$score <- impurity.idx(sp)
      # save the best split as the one that minimizes
      # the impurity score
      if (sp$score < best_split$score) {
        best_split <- sp
      }
    }
  }
  return(best_split)
}

#' A function for computing the impurity of a particular split.
#' @param groups a 2 element list, where the first element are the samples to the left
#' and the second element the samples to the right of a split.
#' @return the gini impurity index.
#' @author Eric Bridgeford
impurity.idx <- function(groups) {
  groups <- lapply(groups[c("left", "right")], function(group) return(group$Y))
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

#' Return Most Probable Class at a Leaf Node
#'
#' @param Y the responses at a leaf node.
#' @return the most probable response at the leaf.
leaf.node <- function(Y) {
  gr.ct <- sapply(levels(Y), function(y) {
    sum(Y == y)
  })
  return(levels(Y)[which.max(gr.ct)])
}
