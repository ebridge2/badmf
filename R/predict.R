#' Decision Tree Prediction
#'
#' A helper method for Decision Tree prediction.
#'
#' @param object a fit decision tree.
#' @param X an [n, p] predictor to predict the class of. Must be specified in an appropriate
#' 2d array, even if you only want a single sample predicted.
#' @return an [n] array of predictions.
#' @author Eric Bridgeford
#' @export
predict.dec.tree <- function(object, X, ...) {
  preds <- sapply(1:dim(X)[1], function(i) {
    return(predict(object$tree, X[i,]))
  })
  return(preds)
}

#' A function for split node prediction
predict.split.node <- function(object, Xi, ...) {
  if (Xi[object$feature] < object$threshold) {
    return(predict(object$left, Xi))
  } else {
    return(predict(object$right, Xi))
  }
}
#' A function for leaf node prediction
predict.leaf.node <- function(object, ...) {
  return(object$vote)
}
