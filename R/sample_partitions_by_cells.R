#' Generate Random, Non-Overlapping Partition of Cells
#'
#' @param reads A data.frame of \eqn{n} reads.
#'
#' @param \dots Argument passed to [sample_partitions_similar_weights].
#'
#' @return A list of random non-overlapping (disjoint) partitions where each
#' element holds indices in \eqn{{1, 2, ..., n}} and where the union of all
#' parts is \eqn{{1, 2, ..., n}}.
#'
#' @references
#' https://github.com/HenrikBengtsson/SegalM_2017-FISH/issues/16
#'
#' @export
sample_partitions_by_cells <- function(reads, ...) {
  stop_if_not("cell_id" %in% colnames(reads))
  cell_weights <- table(reads$cell_id)
  ## Drop non-existing cell_id:s due to empty levels
  cell_weights <- cell_weights[cell_weights > 0]

  ## Partion cells into parts of roughly equal-sized reads
  cell_partitions <- sample_partitions_similar_weights(cell_weights, ...)
  if (length(cell_partitions) == 1L && is.na(cell_partitions)) {
    stop(sprintf("Failed to identify cell partitioning (cell_weights = %g) after %d rejected attempts", cell_weights, attr(cell_partitions, "count", exact = TRUE)))
  }
  
  ## Convert cell partion into read partition
  read_partitions <- lapply(cell_partitions, FUN = function(idxs) {
    which(reads$cell_id %in% names(cell_weights)[idxs])
  })

  ## Sanity check
  idxs <- unlist(read_partitions)
  stop_if_not(length(idxs) == nrow(reads), !any(duplicated(idxs)))
  
  read_partitions
}


#' Generate Random, Non-Overlapping Two-Set Partition of Cells By Half
#'
#' @param reads A data.frame with \eqn{n} reads from \eqn{C} cells.
#'
#' @param fraction A numeric of length two in (0,1/2] specifying the size of
#' each partition relative to `n`, where the first element specifies the
#' "reference" partition and the second the "test" partition.
#' For backward compatible reasons, if of length one, then it's equivalent to
#' specifying `fraction = c(reference = 1/2, test = fraction)`.
#'
#' @param w_tolerance Maximum allowed difference between target weight of
#' each part (e.g. `fraction * 1`) and the actual total weight of the part
#' (i.e. `sum(w[part])`.  There is a one-to-one relationship between
#' weight tolerance and read-count tolerance.  The effective read-count
#' tolerance can be calculated as `fraction * n`.
#' If _all_ parts are within the tolerance, the sampled partition is accepted,
#' otherwise it is rejected.
#'
#' @param max_rejections The maximum number of rejections before giving up.
#'
#' @param warn If `TRUE`, a warning is produced if the parts produced
#' are not of equal size.
#'
#' @param \dots Additional argument passed to
#' [sample_partitions_similar_weights_by_half].
#'
#' @return A list of two disjoint random non-overlapping partitions where each
#' element holds indices in \eqn{{1, 2, ..., n}}.  The first element holds the
#' 'reference' partition and the second the 'test' partition.
#'
#' @references
#' https://github.com/HenrikBengtsson/SegalM_2017-FISH/issues/16
#'
#' @importFrom utils str
#' @export
sample_partitions_by_cells_by_half <- function(reads, fraction, w_tolerance = 0.01, max_rejections = 100L, warn = TRUE, ...) {
  stop_if_not("cell_id" %in% colnames(reads))
  stop_if_not(is.numeric(fraction), length(fraction) %in% 1:2,
              !anyNA(fraction), all(fraction > 0), all(fraction <= 1/2))
  if (length(fraction) == 1L) {
    .Deprecated(msg = "sample_partitions_similar_weights_by_half(): Argument 'fraction' should be of length two")
    fraction <- c(reference = 1/2, test = fraction)
  }
  stopifnot(length(fraction) == 2L)

  cell_weights <- table(reads$cell_id)
  ## Drop non-existing cell_id:s due to empty levels
  cell_weights <- cell_weights[cell_weights > 0]

  ## Partion cells into two partions that have approximately the same number of reads
  cell_sets <- sample_partitions_similar_weights_by_half(cell_weights, fraction = c(reference = 1/2, test = 1/2), w_tolerance = w_tolerance, max_rejections = max_rejections, warn = warn, ...)
  if (length(cell_sets) == 1L && is.na(cell_sets)) {
    stop(sprintf("Failed to identify a cell partition (cell_weights = %g) after %d rejected attempts", cell_weights, attr(cell_sets, "count", exact = TRUE)))
  }
  stop_if_not(length(cell_sets) == 2L)
  
  ## Convert cell partition into read partition
  read_sets <- lapply(cell_sets, FUN = function(cell_idxs) {
    which(reads$cell_id %in% names(cell_weights)[cell_idxs])
  })

  ## Sanity check
  n_sets <- lengths(read_sets)
  n_total <- sum(n_sets)
  props <- n_sets / n_total
  stop_if_not(all(abs(props - 1/2) <= w_tolerance))

  ## Down-sample test set?
  for (ff in seq_along(fraction)) {
    if (fraction[ff] < 1/2) {
      size <- min(floor(fraction[ff] * n_total), n_sets[ff])
      read_idxs <- read_sets[[ff]]
      read_idxs <- read_idxs[sample.int(n_sets[ff], size = size)]
      read_sets[[ff]] <- read_idxs
      read_idxs <- NULL
      
      ## Sanity check
      n_sets <- lengths(read_sets)
      props <- n_sets / n_total
      if (abs(props[ff] - fraction[ff]) > w_tolerance) {
        str(list(read_sets = read_sets, n_sets = n_sets, props = props, fraction_ff = fraction[ff], n_sets = n_sets, n_total = n_total, w_tolerance = w_tolerance))
      }
      stop_if_not(abs(props[ff] - fraction[ff]) <= w_tolerance)
    }
  }

  ## Sanity check
  stop_if_not(length(read_sets) == 2L)
  read_idxs <- unlist(read_sets)
  stop_if_not(all(read_idxs >= 1L), all(read_idxs <= nrow(reads)))
  read_idxs <- NULL
  n_sets <- lengths(read_sets)
  props <- n_sets / n_total
  stop_if_not(all(abs(props - fraction) <= w_tolerance))
  
  read_sets
}
