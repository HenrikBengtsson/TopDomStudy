#' Summarizes the TopDom Best Scores across TopDomOverlapScores
#'
#' @param fit A [base::list] of K `TopDomOverlapScores` objects
#'
#' @param drop_reference If TRUE, then the reference partition is dropped.
#'
#' @return A [base::data.frame] with (K-1) rows.
#'
#' @importFrom stats mad quantile sd
#' @export
overlap_score_summary <- function(fit, drop_reference = TRUE) {
  stopifnot(is.list(fit))
  lapply(fit, FUN = function(x) stopifnot(inherits(x, "TopDomOverlapScores"), length(x) == 1L))
  n <- length(fit)
  stopifnot(n >= 2L)

  attrs <- attributes(fit)
  stopifnot(is.list(attrs), length(attrs) > 1L)
  
  partition_by <- attrs[["partition_by"]]
  stopifnot(length(partition_by) == 1L, !is.na(partition_by), is.character(partition_by))
  if (partition_by %in% c("reads_by_half", "cells_by_half")) {
    stopifnot(n == 2L)
  }

  chr <- attrs[["chromosome"]]
  stopifnot(length(chr) == 1L, !is.na(chr), is.character(chr))

  bin_size <- attrs[["bin_size"]]
  stopifnot(length(bin_size) == 1L, !is.na(bin_size), is.numeric(bin_size), bin_size >= 1L)

  min_cell_size <- attrs[["min_cell_size"]]
  stopifnot(length(min_cell_size) == 1L, !is.na(min_cell_size), is.numeric(min_cell_size), min_cell_size >= 1L)

  ref <- attrs[["reference_partition"]]
  stopifnot(length(ref) == 1L, !is.na(ref), ref >= 1L)

  if (drop_reference) fit <- fit[-ref]
  
  scores <- lapply(fit, FUN = function(x) {
    if (!inherits(x, "TopDomOverlapScores")) return(double(0L))
    y <- x[[chr]]
    y[["best_scores"]]
  })

  summary <- lapply(scores, FUN = function(x) {
    x <- unlist(x, use.names = FALSE)
    count <- length(x)
    x <- x[!is.na(x)]
    n <- length(x)
    data.frame(chromosome = chr, min_cell_size = min_cell_size, bin_size = bin_size, mean = mean(x), as.list(quantile(x)), sd = sd(x), mad = mad(x), count = count, n = n, check.names = FALSE)
  })
  
  summary <- do.call(rbind, summary)

  attr(summary, "partition_by") <- partition_by
  if (!drop_reference) attr(summary, "reference_partition") <- ref
  attr(summary, "seed") <- attrs[["seed"]]

  summary
}
