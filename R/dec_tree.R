#' Bayesian Decision Tree Fit
#'
#' Fit a Bayesian Decision Tree with a \code{stats}-like formula frontend interface.
#'
#' @param formuler ravioli ravioli give me the formuoli.
#' @param data the data associated with the formuler. Note: if you want an intercept, you must
#' add it ahead of time.
#' @param d the number of features to subsample at each node. Defaults to \code{NULL}, which tries every feature.
#' @param alpha the prior parameters for the feature probabilities. A \code{[p]} vector. If \code{NULL}, samples uniformly.
#' Defaults to \code{NULL}.
#' @param method whether you want "classification" or "regression". Defaults to \code{"classification"}.
#' @param depth.max the maximum allowed tree depth. Defaults to \code{5L}.
#' @param size the minimum allowed number of samples for an individual node. Defaults to \code{1L}.
#' @param debug whether to save the predictors and responses that are categorized. Defaults to \code{FALSE}.
#' @param ... trailing arguments.
#' @return an object of class \code{dec.tree.class} containing the following:
#' \item{\code{tree}}{the decision tree.}
#' \item{\code{X}}{The training predictors.}
#' \item{\code{Y}}{the training responses.}
#' \item{\code{d}}{d the number of features subsampled at each node.}
#' \item{\code{alpha}}{the sampling distribution for the features. A \code{[p]} vector.}
#' \item{\code{depth.max}}{the maximum allowed tree depth.}
#' \item{\code{size}}{the maximum allowed tree depth.}
#' \item{\code{debug}}{whether to save the predictors and responses that are categorized.}
#' @author Eric Bridgeford
#' @export
dec_tree.fit <- function(formuler, data=NULL, d=NULL, alpha=NULL, method="classification",
                         depth.max=5L, size=1L, debug=FALSE, ...) {
  # miscellaneous formula jargon
  call <- match.call()

  if (missing(data))
    data <- environment(formuler)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formuler", "data"), names(mf), 0L)
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
    fit <- do.call(tree.class.fit, list(X[,-c(1)], Y, d, alpha, depth.max, size, debug))
    fit$formula <- formuler
  } else {
    stop("Not yet implemented!")
  }
  return(fit)
}

#' Fit a decision tree classifier.
#' @param X the predictors. A \code{[n, p]} matrix.
#' @param Y the responses. A \code{[n]} vector or, optionally, a factor.
#' @param d the number of features to subsample at each node. Defaults to \code{NULL}, which tries every feature.
#' @param alpha the prior parameters for the feature probabilities. A \code{[p]} vector. If \code{NULL}, samples uniformly.
#' Defaults to \code{NULL}.
#' @param depth.max the maximum allowed tree depth. Defaults to \code{5L}.
#' @param size the minimum allowed number of samples for an individual node. Defaults to \code{1L}.
#' @param debug whether to save the predictors and responses that are categorized. Defaults to \code{FALSE}.
#' into a particular leaf node.
#' @return an object of class \code{dec.tree.class} containing the following:
#' \item{\code{tree}}{the decision tree.}
#' \item{\code{X}}{The training predictors.}
#' \item{\code{Y}}{the training responses.}
#' \item{\code{d}}{d the number of features subsampled at each node.}
#' \item{\code{alpha}}{the sampling distribution for the features. A \code{[p]} vector.}
#' \item{\code{depth.max}}{the maximum allowed tree depth.}
#' \item{\code{size}}{the maximum allowed tree depth.}
#' \item{\code{debug}}{whether to save the predictors and responses that are categorized.}
#' @author Eric Bridgeford
#' @export
dec.tree.class.fit <- function(X, Y, d=NULL, alpha=NULL, depth.max=5L, size=1L, debug=FALSE) {
  Y <- as.factor(Y)
  n <- length(Y); p <- dim(X)[2]

  if (!ifelse(is.integer(d), d <= p & d > 0, is.null(d))) {
    stop("d should be a positive integer <= p, or NULL to indicate to sample every feature.")
  }

  if (is.null(d)) {
    d <- as.integer(p)
  }

  if (!ifelse(is.integer(depth.max), depth.max > 0, FALSE)) {
    stop("You have not entered a valid value for depth.max.")
  }

  if (!ifelse(is.integer(size), size > 0, FALSE)) {
    stop("You have not passed a valid value for size.")
  }

  # if alpha is NULL, sample uniformly
  if (!ifelse(is.numeric(alpha), all(alpha > 0), is.null(alpha))) {
    stop("You have not entered a valid Dirichlet prior. Should be NULL, or a p vector with all values should be > 0.")
  }

  # build the tree with the Decision Tree Algorithm
  tree <- build.tree(get.split(X, Y, d, alpha), d, alpha, depth.max, size, 1, debug)
  return(structure(
    list(tree=tree, X=X, Y=Y, p=p, d=d, alpha=alpha, depth.max=depth.max, size=size, debug=debug),
    class="dec.tree.class"
  ))
}

#' Recursive Approach for Building Decision Tree
#' @param split a particular split node.
#' @param d the number of features to subsample.
#' @param alpha the sampling distribution for the features. A \code{[p]} vector. If \code{NULL}, uses a uniform.
#' @param depth.max the max tree depth.
#' @param size the minimum number of elements at a particular.
#' @param depth the current depth of the tree.
#' @param debug a boolean indicating whether to save the predictors and responses
#' that are categorized into leaf nodes.
#' @importFrom forcats fct_c
#' @return a layer of a decision tree.
#' @author Eric Bridgeford
build.tree <- function(split, d, alpha, depth.max, size, depth, debug=FALSE) {
  # if we have an empty child, create a leaf by merging
  # so that we don't have a totally empty child
  if (length(split$left$Y) == 0 || length(split$right$Y) == 0) {
    if (!debug) {
      return(leaf.node(fct_c(split$left$Y, split$right$Y)))
    } else {
      return(c(leaf.node(fct_c(split$left$Y, split$right$Y)),
               list(X=rbind(split$left$X, split$right$X),
                    Y=fct_c(split$left$Y, split$right$Y))))
    }
  }

  # if we are at the max depth, create a leaf for the left and
  # right children
  if (depth >= depth.max) {
    if (!debug) {
      split$left <- leaf.node(split$left$Y)
    } else {
      split$left <- c(leaf.node(split$left$Y),
                      list(X=split$left$X, Y=split$left$Y))
    }

    if (!debug) {
      split$right <- leaf.node(split$right$Y)
    } else {
      split$right <- c(leaf.node(split$right$Y),
                       list(X=split$right$X, Y=split$right$Y))
    }
    return(split)
  }
  # process right child
  if (length(split$right$Y) > size) {
    # split the node if we can still do better
      split$right <- build.tree(
        get.split(split$right$X, split$right$Y, d, alpha),
        d, alpha, depth.max, size, depth + 1, debug=debug
      )
  } else {
    if (!debug) {
      split$right <- leaf.node(split$right$Y)
    } else {
      split$right <- c(leaf.node(split$right$Y),
                       list(X=split$right$X, Y=split$right$Y))
    }
  }

  # process left child
  if (length(split$left$Y) > size) {
    # split the node if we can still do better
    split$left <- build.tree(
      get.split(split$left$X, split$left$Y, d, alpha),
      d, alpha, depth.max, size, depth + 1, debug=debug
    )
  } else {
    if (!debug) {
      split$left <- leaf.node(split$left$Y)
    } else {
      split$left <- c(leaf.node(split$left$Y),
                      list(X=split$left$X, Y=split$left$Y))
    }
  }
  return(split)
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
  split <- structure(list(left=list(X=X[idx,], Y=Y[idx]),
                     right=list(X=X[!idx,], Y=Y[!idx]),
                     feature=i, threshold=t, n=length(Y)),
                     class="split.node")
  return(split)
}

#' Find Best split in the data.
#'
#' @param X the predictors.
#' @param Y the responses.
#' @param d the number of features to subsample. Should be an integer 0 < d <= p.
#' @param alpha the sampling distribution for the features. A [p] vector. If \code{NULL}, uses a uniform.
#' @return the best split.
#' @importFrom MCMCpack rdirichlet
#' @author Eric Bridgeford
get.split <- function(X, Y, d, alpha) {
  n <- length(Y); p <- dim(X)[2]
  # sample features to check
  if (is.null(alpha)) {
    features <- sample(1:p, d, replace=FALSE)
  } else {
    feat.probs <- rdirichlet(1, alpha)
    features <- sample(1:p, d, replace=FALSE, prob=feat.probs)
  }
  # initialize best split
  best_split <- list(X=X, Y=Y, t=NULL, feature=NULL, score=1.1, split=NULL)
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
#' @param Y the responses at a leaf node, as a \code{[n]} factor.
#' @return the most probable response at the leaf.
#' @author Eric Bridgeford
leaf.node <- function(Y) {
  gr.ct <- sapply(levels(Y), function(y) {
    sum(Y == y)
  })
  return(structure(list(vote=factor(levels(Y)[which.max(gr.ct)],
                                    levels = levels(Y))), class="leaf.node"))
}

#' Count Features in a Decision Tree
#'
#' A function to count the feature utilization in a Decision Tree.
#'
#' @param object an object of class \code{dec.tree.class}.
#' @param ... trailing args.
#' @return a vector containing the number of times each feature
#' is used at a split node.
#' @author Eric Bridgeford
#' @export
count.features.dec.tree.class <- function(object, ...) {
  feat.usage <- table(count.features(object$tree))
  counts <- rep(0, object$p)
  names(counts) <- 1:object$p
  counts[names(feat.usage)] <- as.numeric(feat.usage)
  return(counts)
}

#' Count Features Leaf
#'
#' A default to return NULL and stop recursion when counting features
#' in a stump once we reach a leaf node.
#'
#' @param object an object of class \code{leaf.node}.
#' @param ... trailing args.
#' @return a vector containing the number of times each feature is used
#' at a split node.
#' @author Eric Bridgeford
count.features.leaf.node <- function(object, ...) {
  return(NULL)
}

#' Count Features Stump
#'
#' A function to count the Usage of Features At or Below a Split Node (a Stump)
#'
#' @param object an object of class \code{split.node}.
#' @param ... trailing args.
#' @return a vector containing the number of times each feature is used
#' at a split node.
#' @author Eric Bridgeford
count.features.split.node <- function(object, ...) {
  return(c(object$feature, count.features(object$left), count.features(object$right)))
}
