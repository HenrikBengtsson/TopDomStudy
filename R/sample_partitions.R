#' Generate Random, Non-Overlapping Partitions of Indices
#'
#' @param n Number of elements to partition.
#'
#' @param fraction A numeric in (0,1] specifying the size of each partition
#'                 relative to `n`.  If `NULL`, argument `size` is used.
#'
#' @param size An integer in \eqn{{1, 2, ..., n}} specifying the size of each
#'             partition.  If `NULL`, argument `fraction` is used.
#'
#' @param warn If `TRUE`, a warning is produced if the partitions produced
#'             are not of equal size.
#'
#' @return A list of random non-overlapping (disjoint) partitions where each
#' element holds indices in \eqn{{1, 2, ..., n}} and where the union of all
#' partitions is \eqn{{1, 2, ..., n}}.
#'
#' @importFrom parallel splitIndices
#' @export
sample_partitions <- function(n, fraction = NULL, size = NULL, warn = TRUE) {
  stop_if_not(is.numeric(n), length(n) == 1L, !is.na(n), n > 0)
  stop_if_not(!is.null(fraction) || !is.null(size))

  ## Number of partitions
  if (!is.null(fraction)) {
    stop_if_not(is.numeric(fraction), length(fraction) == 1L, !is.na(fraction),
                fraction > 0, fraction <= 1)
    npartitions <- as.integer(ceiling(1 / fraction))
  } else if (!is.null(size)) {
    stop_if_not(is.numeric(size), length(size) == 1L, !is.na(size),
                size >= 1, size <= n)
    npartitions <- as.integer(ceiling(n / size))
  }

  ## Sanity checks
  stop_if_not(is.integer(npartitions), length(npartitions) == 1L,
              !is.na(npartitions), npartitions > 0L, npartitions <= n)

  ## Produce shuffled indices in 1:n
  idxs <- sample.int(n, size = n, replace = FALSE)

  partitions <- splitIndices(nx = n, ncl = npartitions)
  partitions <- lapply(partitions, FUN = function(partition) idxs[partition])

  if (warn) {
    sizes <- unique(lengths(partitions))
    if (length(sizes) > 1L) {
      sizes <- sort(sizes)
      warning(sprintf("sample_partitions(n = %g, ...) produced %d partitions of unequal sizes: %d and %d",
              n, length(partitions), sizes[1], sizes[2]))
    }
  }
  
  partitions
}



#' Generate Random, Non-Overlapping Similarly-Weighted Partitions of Indices
#'
#' @param w Numeric vector of `n` non-negative, finite weights.
#'          Weights are normalized such that `sum(w)` equals one.
#'
#' @param fraction A numeric in (0,1] specifying the size of each partition
#'                 relative to `n`.  If `NULL`, argument `size` is used.
#'
#' @param size An integer in \eqn{{1, 2, ..., n}} specifying the size of each
#'             partition.  If `NULL`, argument `fraction` is used.
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
sample_partitions_similar_weights <- function(w, fraction = NULL, size = NULL, w_tolerance = 0.01, max_rejections = 100L, warn = TRUE) {
  stop_if_not(is.numeric(w), length(w) > 0, !anyNA(w), all(w > 0))
  stop_if_not(is.numeric(w_tolerance), length(w_tolerance) == 1L,
              !is.na(w_tolerance), w_tolerance >= 0, w_tolerance <= 1)
  stop_if_not(!is.null(fraction) || !is.null(size))
  stop_if_not(is.numeric(max_rejections), length(max_rejections) == 1L,
              !is.na(max_rejections), max_rejections >= 0)

  ## Normalize weights
  w <- w / sum(w)
  n <- length(w)
  
  ## Number of partitions
  if (!is.null(fraction)) {
    stop_if_not(is.numeric(fraction), length(fraction) == 1L, !is.na(fraction),
                fraction > 0, fraction <= 1)
    npartitions <- as.integer(ceiling(1 / fraction))
  } else if (!is.null(size)) {
    stop_if_not(is.numeric(size), length(size) == 1L, !is.na(size),
                size >= 1, size <= n)
    npartitions <- as.integer(ceiling(n / size))
  }

  ## Sanity checks
  stop_if_not(is.integer(npartitions), length(npartitions) == 1L,
              !is.na(npartitions), npartitions > 0L, npartitions <= n)


  ## Partioning
  pidxs <- splitIndices(nx = n, ncl = npartitions)

  ## Target weight per partion
  w_target <- 1 / npartitions

  idxs <- NULL
  ready <- FALSE
  count <- 0L
  while (!ready && count < max_rejections) {
    ## Produce shuffled indices in 1:n
    idxs <- sample.int(n, size = n, replace = FALSE)

    ## Shuffle weights accordingly
    w_idxs <- w[idxs]

    ## Calculate weights per partion
    ws <- lapply(pidxs, FUN = function(partition) sum(w_idxs[partition]))
    ws <- unlist(ws, use.names = FALSE)
    ready <- all(abs(ws - w_target) <= w_tolerance)
    count <- count + 1L
  }

  if (ready) {
    partitions <- lapply(pidxs, FUN = function(partition) idxs[partition])
    attr(partitions, "weights") <- ws
  } else {
    partitions <- NA
  }
  attr(partitions, "w_target") <- w_target
  attr(partitions, "w_tolerance") <- w_tolerance
  attr(partitions, "counts") <- count
  
  partitions
}
