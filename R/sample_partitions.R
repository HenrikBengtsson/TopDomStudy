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



#' Generate Random, Non-Overlapping Partitions of Indices with Reference and Sequence of Partitions
#'
#' @param n Number of elements to partition.
#'
#' @param fraction A numeric in (0,1] specifying the size of each partition
#'                 relative to `n`.  If `NULL`, argument `size` is used.
#'
#' @param size An integer in \eqn{{1, 2, ..., n}} specifying the size of each
#'             partition.  If `NULL`, argument `fraction` is used.
#'
#' @param seq A sequence of fractions (if `fraction` is specified) or
#' a sequence of sizes (if `size` is specified)
#'
#' @return A list of random non-overlapping partitions where each
#' element holds indices in \eqn{{1, 2, ..., n}}.
#' The first element holds the `reference` partition \eqn{P_r} with \eqn{|P_r|}
#' indices and the remaining elements holds partitions with indices in
#' \eqn{{1, 2, ..., n} \ P_r} where the size of the partitions corresponds
#' to the sizes specified by `seq`.
#' The reference partition and the "remaining" partitions are always disjoint,
#' but the "remaining" partitions maybe be overlapping.
#'
#' @importFrom parallel splitIndices
#' @export
sample_partitions_ref_vs_seq <- function(n, fraction = NULL, size = NULL, seq) {
  stop_if_not(is.numeric(n), length(n) == 1L, !is.na(n), n > 0)
  stop_if_not(!is.null(fraction) || !is.null(size))
  stop_if_not(is.numeric(seq), length(seq) >= 1L, !anyNA(seq), all(seq > 0))

  ## Number of partitions
  if (!is.null(fraction)) {
    stop_if_not(is.numeric(fraction), length(fraction) == 1L, !is.na(fraction),
                fraction > 0, fraction <= 1)
    size <- round(fraction * n)
    stop_if_not(all(seq <= 1 - fraction))
    seq_labels <- sprintf("fraction=%g", seq)
    seq_sizes <- round(seq * n)
  } else if (!is.null(size)) {
    stop_if_not(is.numeric(size), length(size) == 1L, !is.na(size),
                size >= 1, size <= n)
    size <- round(size)
    stop_if_not(all(seq <= n - size))
    seq_sizes <- round(seq)
    seq_labels <- sprintf("size=%g", seq_sizes)
  }
  
  ## Sanity checks
  stop_if_not(is.numeric(size), length(size) == 1L, !is.na(size),
              size >= 1, size <= n)
  stop_if_not(is.numeric(seq_sizes), length(seq_sizes) >= 1L, !anyNA(seq_sizes),
              all(seq_sizes >= 1), all(seq_sizes <= n - size))

  ## Produce shuffled indices in 1:n
  idxs <- sample.int(n, size = n, replace = FALSE)

  partitions <- vector("list", length = 1L + length(seq_sizes))
  names(partitions) <- c("reference", seq_labels)

  ## Reference partition
  partitions$reference <- idxs[1:size]

  ## Remaining indices to subsample from
  idxs <- idxs[(size+1):n]
  nidxs <- length(idxs)

  for (kk in seq_along(seq_sizes)) {
    idxs_kk <- idxs[sample(nidxs, size = seq_sizes[kk])]
    partitions[[kk + 1L]] <- idxs_kk
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
#' @param seq A sequence of fractions (if `fraction` is specified) or
#' a sequence of sizes (if `size` is specified)
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
sample_partitions_similar_weights_ref_vs_seq <- function(w, fraction = NULL, size = NULL, w_tolerance = 0.01, max_rejections = 100L, seq, warn = TRUE) {
  stop_if_not(is.numeric(w), length(w) > 0, !anyNA(w), all(w > 0))
  stop_if_not(is.numeric(w_tolerance), length(w_tolerance) == 1L,
              !is.na(w_tolerance), w_tolerance >= 0, w_tolerance <= 1)
  stop_if_not(!is.null(fraction) || !is.null(size))
  stop_if_not(is.numeric(max_rejections), length(max_rejections) == 1L,
              !is.na(max_rejections), max_rejections >= 0)
  stop_if_not(is.numeric(seq), length(seq) >= 1L, !anyNA(seq), all(seq > 0))

  ## Normalize weights
  w <- w / sum(w)
  n <- length(w)

  ## Number of partitions
  if (!is.null(fraction)) {
    stop_if_not(is.numeric(fraction), length(fraction) == 1L, !is.na(fraction),
                fraction > 0, fraction <= 1)
    size <- round(fraction * n)
    stop_if_not(all(seq <= 1 - fraction))
    seq_labels <- sprintf("fraction=%g", seq)
    seq_sizes <- round(seq * n)
  } else if (!is.null(size)) {
    stop_if_not(is.numeric(size), length(size) == 1L, !is.na(size),
                size >= 1, size <= n)
    size <- round(size)
    stop_if_not(all(seq <= n - size))
    seq_sizes <- round(seq)
    seq_labels <- sprintf("size=%g", seq_sizes)
  }

  ## Sanity checks
  stop_if_not(is.numeric(size), length(size) == 1L, !is.na(size),
              size >= 1, size <= n)
  stop_if_not(is.numeric(seq_sizes), length(seq_sizes) >= 1L, !anyNA(seq_sizes),
              all(seq_sizes >= 1), all(seq_sizes <= n - size))

  ## Target weights for all partions
  w_targets <- c(size, seq_sizes) / n
  stop_if_not(is.numeric(w_targets), !anyNA(w_targets),
              all(w_targets >= 0), all(w_targets <= 1))
  print(list(w_targets = w_targets))
  
  partitions <- vector("list", length = 1L + length(seq_sizes))
  names(partitions) <- c("reference", seq_labels)
  ws <- rep(NA_real_, times = length(partitions))
  
  idxs <- NULL
  ready <- FALSE
  count <- 0L
  while (!ready && count < max_rejections) {
    ## Produce shuffled indices in 1:n
    idxs <- sample.int(n, size = n, replace = FALSE)

    ## Shuffle weights accordingly
    w_idxs <- w[idxs]

    ## Reference partition
    idxs_ref <- idxs[1:size]
    w_ref <- sum(w_idxs[idxs_ref])
    
    ## Remaining indices to subsample from
    idxs <- idxs[(size+1):n]
    nidxs <- length(idxs)

    ok <- TRUE
    for (kk in seq_along(seq_sizes)) {
      idxs_kk <- idxs[sample(nidxs, size = seq_sizes[kk])]
      w_kk <- sum(w_idxs[idxs_kk])
      w_target <- w_targets[kk + 1L]
      ok <- (abs(w_kk - w_target) <= w_tolerance)
#      str(list(kk = kk, idxs = idxs, w_ref = w_ref, idxs_kk = idxs_kk, w_kk = w_kk, w_target = w_target, ok = ok))
      if (!ok) break
      partitions[[kk + 1L]] <- idxs_kk
      ws[kk + 1L] <- w_kk
    }

    count <- count + 1L
    
    if (ok) {
      partitions[[1L]] <- idxs_ref
      ws[1L] <- w_ref
  
      ## Final assertion that weights are similar
      ready <- all(abs(ws - w_targets) <= w_tolerance)
      if (ready) break
    }    
  }

  if (ready) {
    attr(partitions, "weights") <- ws
  } else {
    partitions <- NA
  }
  attr(partitions, "w_targets") <- w_targets
  attr(partitions, "w_tolerance") <- w_tolerance
  attr(partitions, "counts") <- count

  partitions
}
