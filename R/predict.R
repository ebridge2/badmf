#' Decision Tree Prediction
#'
#' A method for Decision Tree prediction. Predicts by pushing samples down
#' to a leaf node based on the split node criterion.
#'
#' @param object a fit decision tree.
#' @param X an \code{[n, p]} predictor to predict the class of. Must be specified in an appropriate
#' 2d array, even if you only want a single sample predicted.
#' @return an \code{[n]} array of predictions.
#' @author Eric Bridgeford
#' @export
predict.dec.tree.class <- function(object, X, ...) {
  preds <- sapply(1:dim(X)[1], function(i) {
    return(predict(object$tree, X[i,]))
  })
  return(preds)
}

#' A function for split node prediction
#' pushes point down the tree recursively
#' until it ends at a leaf node.
#' @param object the stump to use.
#' @param Xi the poinbt.
#' @return the predic
predict.split.node <- function(object, Xi, ...) {
  if (Xi[object$feature] < object$threshold) {
    return(predict(object$left, Xi))
  } else {
    return(predict(object$right, Xi))
  }
}

#' A function for leaf node prediction
#' @param object a leaf.node
#' @return the vote for the leaf.
predict.leaf.node <- function(object, ...) {
  return(object$vote)
}

#' Random Forest Prediction
#'
#' A method for Random Forest prediction. Predicts by, for each tree in the forest,
#' doing a standard decision tree prediction. Then, the most likely vote across the trees
#' is used for predictions.
#'
#' If the two or more classes are tied for the maximum number of votes across the trees,
#' the predicted class is chosen amongst the classes tied for the maximal number of votes
#' at random.
#'
#' @param object a fit object of class rf.class.
#' @param X an \code{[n, p]} predictor to predict the class of. Must be specified in an appropriate
#' 2d array, even if you only want a single sample predicted.
#' @return an \code{[n]} array of predictions.
#' @author Eric Bridgeford
#' @export
predict.rf.class <- function(object, X, ...) {
  preds <- sapply(1:dim(X)[1], function(i) {
    votes.xi <- sapply(object$forest, function(tree) {
      predict(tree, X[i,,drop=FALSE])
    })
    vote.ct <- table(votes.xi)
    maxima <- which(vote.ct == max(vote.ct, na.rm=TRUE))
    if (length(maxima) > 1) {
      # if there is a tie, pick at random
      maxima <- sample(maxima, 1)
    }
    return(names(vote.ct)[maxima])
  })
  return(preds)
}
