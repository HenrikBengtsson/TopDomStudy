#' Generate Random, Non-Overlapping Partitions of Cells
#'
#' @param reads A data.frame of reads.
#'
#' @param \dots Argument passed to [sample_partitions_similar_weights].
#'
#' @return A list of random non-overlapping (disjoint) partitions where each
#' element holds indices in \eqn{{1, 2, ..., n}} and where the union of all
#' partitions is \eqn{{1, 2, ..., n}}.
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

  ## Partion cells into partions of roughly equal-sized reads
  cell_partitions <- sample_partitions_similar_weights(cell_weights, ...)
  if (length(cell_partitions) == 1L && is.na(cell_partitions)) {
    stop(sprintf("Failed to identify cell partitioning (cell_weights = %g) after %d rejected attempts", cell_weights, attr(cell_partitions, "count")))
  }
  
  ## Convert into cell partions into read partitions
  read_partitions <- lapply(cell_partitions, FUN = function(idxs) {
    which(reads$cell_id %in% names(cell_weights)[idxs])
  })

  ## Sanity check
  idxs <- unlist(read_partitions)
  stop_if_not(length(idxs) == nrow(reads), !any(duplicated(idxs)))
  
  read_partitions
}


#' Generate Random, Non-Overlapping Partitions of Cells By Half
#'
#' @param reads A data.frame of reads.
#'
#' @param \dots Argument passed to [sample_partitions_similar_weights_by_half].
#'
#' @return A list of random non-overlapping (disjoint) partitions where each
#' element holds indices in \eqn{{1, 2, ..., n}} and where the union of all
#' partitions is \eqn{{1, 2, ..., n}}.
#'
#' @references
#' https://github.com/HenrikBengtsson/SegalM_2017-FISH/issues/16
#'
#' @export
sample_partitions_by_cells_by_half <- function(reads, ...) {
  stop_if_not("cell_id" %in% colnames(reads))
  cell_weights <- table(reads$cell_id)
  ## Drop non-existing cell_id:s due to empty levels
  cell_weights <- cell_weights[cell_weights > 0]

  ## Partion cells into partions of roughly equal-sized reads
  cell_partitions <- sample_partitions_similar_weights_by_half(cell_weights, ...)
  if (length(cell_partitions) == 1L && is.na(cell_partitions)) {
    stop(sprintf("Failed to identify cell partitioning (cell_weights = %g) after %d rejected attempts", cell_weights, attr(cell_partitions, "count")))
  }
  
  ## Convert into cell partions into read partitions
  read_partitions <- lapply(cell_partitions, FUN = function(idxs) {
    which(reads$cell_id %in% names(cell_weights)[idxs])
  })

  ## Sanity check
  idxs <- unlist(read_partitions)
  stop_if_not(all(idxs >= 1L), all(idxs <= nrow(reads)))
  
  read_partitions
}
