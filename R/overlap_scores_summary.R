#' Summarizes the TopDom Best Scores across TopDomOverlapScores
#'
#' @param fit A [base::list] of K `TopDomOverlapScores` objects
#'
#' @param weights A character string specifying how overlap scores across
#' domains should be weighted.
#'
#' @param drop_reference If TRUE, then the reference partition is dropped.
#'
#' @return A [base::data.frame] with (K-1) rows.
#'
#' @importFrom stats mad quantile sd
#' @importFrom matrixStats weightedMean weightedSd weightedMad
#' @importFrom Hmisc wtd.quantile
#' @export
overlap_score_summary <- function(fit, weights = c("uniform", "by_length"), drop_reference = TRUE) {
  stopifnot(is.list(fit))
  lapply(fit, FUN = function(x) stopifnot(inherits(x, "TopDomOverlapScores"), length(x) == 1L))
  weights <- match.arg(weights, choices = c("uniform", "by_length"))
  
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
  
  data <- lapply(fit, FUN = function(x) {
    if (!inherits(x, "TopDomOverlapScores")) return(data.frame(best_score = double(0L), best_length = integer(0L)))
    y <- x[[chr]]
    y[c("best_score", "best_length")]
  })

  summary <- lapply(data, FUN = function(xy) {
    stopifnot(is.data.frame(xy))
    count <- nrow(xy)
    score <- xy[["best_score"]]
    length <- xy[["best_length"]]

    ## Drop missing values
    keep <- !is.na(score) & !is.na(length) & length > 0L
    score <- score[keep]
    length <- length[keep]
    n <- length(score)

    if (weights == "uniform") {
      mean_hat <- mean(score)
      quantile_hat <- quantile(score)
      sd_hat <- sd(score)
      mad_hat <- mad(score)
    } else if (weights == "by_length") {
      mean_hat <- weightedMean(score, w = length)
      quantile_hat <- wtd.quantile(score, weights = length)
      ## Make names consistent with stats::quantile() names
      names(quantile_hat) <- gsub("^[ ]+", "", names(quantile_hat))
      sd_hat <- weightedSd(score, w = length) 
      mad_hat <- weightedMad(score, w = length)
    }
    
    data.frame(
      chromosome = chr,
      min_cell_size = min_cell_size,
      bin_size = bin_size,
      mean = mean_hat,
      as.list(quantile_hat),
      sd = sd_hat,
      mad = mad_hat,
      count = count,
      n = n,
      check.names = FALSE
    )
  })
  
  summary <- do.call(rbind, summary)

  attr(summary, "weights") <- weights
  attr(summary, "partition_by") <- partition_by
  if (!drop_reference) attr(summary, "reference_partition") <- ref
  attr(summary, "seed") <- attrs[["seed"]]

  summary
}
