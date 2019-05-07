#' Counting Feature Utilization
#'
#' A generic method for counting features in a
#' Tree or Forest.
#'
#' @param object a model object for which feature counting
#' is desired.
#' @param ... additional arguments impacting feature counting.
#' @return the counts of features used in the tree or forest.
#' @author Eric Bridgeford
#' @export
count.features <- function(object, ...) {
  UseMethod("count.features")
}
