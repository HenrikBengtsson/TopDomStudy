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
  idxs <- parts[[2]]
  nidxs <- length(idxs)
  size <- min(round(fraction * n), nidxs)
  parts[[2]] <- idxs[sample.int(nidxs, size = size)]
  names(parts) <- c("reference", sprintf("fraction=%g", fraction))
  attr(parts, "fraction") <- fraction
  attr(parts, "n") <- n
  
  parts
}



#' Generate Random, Non-Overlapping, Similarly-Weighted, Two-Set Partition of Indices
#'
#' @param w Numeric vector of `n` non-negative, finite weights.
#' Weights are normalized such that `sum(w)` equals one.
#'
#' @param fraction A numeric in (0,1/2] specifying the size of the second
#' set relative to `n`.  If `fraction = 1/2`, then there will be no down
#' sampling performed and the two 50%-sized sets are returned as is.
#'
#' @param w_tolerance Maximum allowed difference between target weight of
#' each of the two sets (e.g. `fraction * 1`) and the actual total weight
#' of the sets (i.e. `sum(w[set])`.  If _both_ sets are within the tolerance,
#' the partition is accepted, otherwise rejected.
#'
#' @param max_rejections The maximum number of rejections before giving up.
#'
#' @param warn If `TRUE`, a warning is produced if the partitions produced
#'             are not of equal size.
#'
#' @return A list of two random non-overlapping (disjoint) sets where each
#' element holds indices in \eqn{{1, 2, ..., n}} and where their union is
#' \eqn{{1, 2, ..., n}}.  Attribute `weights` gives the total normalized
#' weight of each set.
#' If no accepted sample was found, then `NA` is returned.
#' Attribute `count` gives the number of internal samples produced before
#' arriving at an accepted or rejected partition.
#'
#' @importFrom parallel splitIndices
#' @export
sample_partitions_similar_weights_by_half <- function(w, fraction, w_tolerance = 0.01, max_rejections = 100L, warn = TRUE) {
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

  ## Failed to find a solution?
  if (!is.list(parts) && is.na(parts)) return(parts)

  ## Sanity check
  stop_if_not(length(parts) == 2L)

  ## Down-sample test set?
  if (fraction < 1/2) {
    idxs <- parts[[2]]
    nidxs <- length(idxs)
    size <- min(round(fraction * n), nidxs)
    idxs <- idxs[sample.int(nidxs, size = size)]
    parts[[2]] <- idxs
  }
  
  names(parts) <- c("reference", sprintf("fraction=%g", fraction))
  attr(parts, "fraction") <- fraction
  attr(parts, "n") <- n
  
  parts
}
