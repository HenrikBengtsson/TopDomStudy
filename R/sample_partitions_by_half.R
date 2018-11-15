#' Generate Random, Non-Overlapping Partitions of Indices with Reference vs Rest
#'
#' @param n Number of elements to partition.
#'
#' @param fraction A numeric in (0,1/2] specifying the size of each partition
#'                 relative to `n`.
#'
#' @param warn If `TRUE`, a warning is produced if the partitions produced
#'             are not of equal size.

#' @return A list of disjoint random non-overlapping partitions where each
#' element holds indices in \eqn{{1, 2, ..., n}}.
#' The first element holds the `reference` partition \eqn{P_0.5} with \eqn{n/2}
#' indices and the remaining elements holds partitions with indices in
#' \eqn{{1, 2, ..., n} \ P_0.5} where the size of the partitions corresponds
#' to the sizes specified by `fraction`.
#'
#' @export
sample_partitions_by_half <- function(n, fraction, warn = TRUE) {
  stop_if_not(is.numeric(n), length(n) == 1L, !is.na(n), n > 0)
  stop_if_not(is.numeric(fraction), length(fraction) == 1L, !is.na(fraction),
              fraction > 0, fraction <= 1/2)

  parts <- sample_partitions(n = n, fraction = 1/2)

  partitions <- list(reference = parts[[1]])
  idxs <- parts[[2]]
  parts <- sample_partitions(n = length(idxs), fraction = 2 * fraction, warn = warn)
  parts <- lapply(parts, FUN = function(i) idxs[i])
  idxs <- NULL
  partitions <- c(partitions, parts)
  names(partitions)[1] <- "reference"
  parts <- NULL

  ## Sanity checks
  idxs <- unlist(partitions, use.names = FALSE)
  stop_if_not(length(idxs) == n)
  idxs <- sort(unique(idxs))
  stop_if_not(length(idxs) == n, all(idxs == seq_len(n)))
  
  partitions
}



#' Generate Random, Non-Overlapping Similarly-Weighted Partitions of Indices
#'
#' @param w Numeric vector of `n` non-negative, finite weights.
#'          Weights are normalized such that `sum(w)` equals one.
#'
#' @param fraction A numeric in (0,1] specifying the size of each partition
#'                 relative to `n`.
#'
#' @param w_tolerance Maximum allowed difference between target weight of
#'        each partition (e.g. `fraction * 1`) and the actual total weight
#'        of the partion (i.e. `sum(w[partition])`.  If _all_ partions are
#'        within the tolerance, the sample is accepted, otherwise rejected.
#'
#' @param max_rejections The maximum number of rejections before giving up.
#'
#' @param warn If `TRUE`, a warning is produced if the partitions produced
#'             are not of equal size.
#'
#' @return A list of random non-overlapping (disjoint) partitions where each
#' element holds indices in \eqn{{1, 2, ..., n}} and where the union of all
#' partitions is \eqn{{1, 2, ..., n}}.  Attribute `weights` gives the total
#' normalized weight of each partition.  Attribute `count` gives the number
#' of internal samples produced before arriving at an accepted sample.
#' If no accepted sample was found, the `NA` is returned
#' (with `count` attribute set).
#'
#' @importFrom parallel splitIndices
#' @export
sample_partitions_similar_weights_by_half <- function(w, fraction = NULL, w_tolerance = 0.01, max_rejections = 100L, warn = TRUE) {
  stop_if_not(is.numeric(w), length(w) > 0, !anyNA(w), all(w > 0))
  stop_if_not(is.numeric(w_tolerance), length(w_tolerance) == 1L,
              !is.na(w_tolerance), w_tolerance >= 0, w_tolerance <= 1)
  stop_if_not(is.numeric(fraction), length(fraction) == 1L, !is.na(fraction),
              fraction > 0, fraction <= 1/2)
  stop_if_not(is.numeric(max_rejections), length(max_rejections) == 1L,
              !is.na(max_rejections), max_rejections >= 0)

  ## Normalize weights
  w <- w / sum(w)
  n <- length(w)

  parts <- sample_partitions_similar_weights(w = w, fraction = 1/2, w_tolerance = w_tolerance, max_rejections = max_rejections, warn = warn)
  if (!is.list(parts) && is.na(parts)) return(parts)

  partitions <- list(reference = parts[[1]])
  idxs <- parts[[2]]
  w <- w[idxs]
  
  parts <- sample_partitions_similar_weights(w = w, fraction = 2 * fraction, w_tolerance = w_tolerance, max_rejections = max_rejections, warn = warn)
  if (!is.list(parts) && is.na(parts)) return(parts)
  
  parts <- lapply(parts, FUN = function(i) idxs[i])
  idxs <- NULL
  partitions <- c(partitions, parts)
  attributes(partitions) <- attributes(parts)
  names(partitions)[1] <- "reference"
  parts <- NULL

  ## Sanity checks
  idxs <- unlist(partitions, use.names = FALSE)
  stop_if_not(length(idxs) == n)
  idxs <- sort(unique(idxs))
  stop_if_not(length(idxs) == n, all(idxs == seq_len(n)))
  
  partitions
}
