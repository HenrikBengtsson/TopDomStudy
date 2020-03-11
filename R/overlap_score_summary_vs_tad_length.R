#' @importFrom tibble as_tibble
#' @importFrom utils file_test
#' @importFrom R.cache loadCache saveCache
read_overlap_score_summary_vs_tad_length <- function(dataset, chromosome, bin_sizes, rhos, reference_rhos = rep(1/2, times = length(rhos)), window_size = 5L, nsamples = 50L, weights = c("by_length", "uniform"), domain_length = NULL, path = "overlapScoreSummary", force = FALSE, ..., verbose = FALSE) {
  chromosome <- as.integer(chromosome)
  bin_sizes <- as.integer(bin_sizes)
  stopifnot(is.numeric(rhos), !anyNA(rhos), all(rhos > 0), all(rhos <= 1/2))
  stopifnot(is.numeric(reference_rhos), !anyNA(reference_rhos), all(reference_rhos > 0), all(reference_rhos <= 1/2))
  stopifnot(length(reference_rhos) == length(rhos), all(reference_rhos >= rhos))
  window_size <- as.integer(window_size)
  nsamples <- as.integer(nsamples)
  weights <- match.arg(weights)
  
  if (!file_test("-d", path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    stop_if_not(file_test("-d", path))
  }

  ## Tags
  chromosome_tag <- sprintf("chr=%s", chromosome)
  window_size_tag <- sprintf("window_size=%d", window_size)
  if (!is.null(domain_length)) {
    stop_if_not(is.numeric(domain_length), length(domain_length) == 2L, !anyNA(domain_length), all(domain_length > 0))
    domain_length_tag <- sprintf("domain_length=%.0f-%.0f", domain_length[1], domain_length[2])
  } else {
    domain_length_tag <- NULL
  }
  weights_tag <- sprintf("weights=%s", weights)
  nsamples_tag <- sprintf("nsamples=%d", nsamples)

  if (verbose) {
    message("read_overlap_score_summary_vs_tad_length() ...")
    message("- chromosome: ", chromosome)
    message("- window_size: ", window_size)
    message("- weights: ", weights)
    message("- nsamples: ", nsamples)
  }

  key <- list(dataset = dataset, chromosome = chromosome, bin_sizes = sort(bin_sizes), rhos = rhos, reference_rhos = reference_rhos, window_size = window_size, nsamples = nsamples, weights = weights, domain_length = domain_length)
  dirs <- c("TopDomStudy", dataset)
  if (!force) {
    summary <- loadCache(key = key, dirs = dirs)
    if (!is.null(summary)) {
      if (verbose) {
        message("read_overlap_score_summary_vs_fraction() ... cached")
      }
      return(summary)
    }
  }

  dim <- c(length(bin_sizes), length(rhos))
  summary <- vector("list", length = prod(dim))
  dim(summary) <- dim
  for (bb in seq_along(bin_sizes)) {
    bin_size <- bin_sizes[bb]
    bin_size_tag <- sprintf("bin_size=%.0f", bin_size)
    if (verbose) message(sprintf("Bin size #%d (%s) of %d ...", bb, bin_size_tag, length(bin_sizes)))
    
    for (rr in seq_along(rhos)) {
      rho <- rhos[rr]
      reference_rho <- reference_rhos[rr]
      test_tag <- sprintf("test=%.3f", rho)
      reference_tag <- sprintf("reference=%.3f", reference_rho)
      if (verbose) message(sprintf("Fraction #%d (%s and %s with %s bps on Chr %s) of %d ...", rr, test_tag, reference_tag, bin_size, chromosome, length(rhos)))

      tags <- c(chromosome_tag, "cells_by_half", "avg_score", bin_size_tag, test_tag, reference_tag, window_size_tag, domain_length_tag, weights_tag, nsamples_tag)
  
      fullname <- paste(c(dataset, tags), collapse = ",")
      pathname_summary_kk <- file.path(path, sprintf("%s.rds", fullname))
      if (verbose) message("pathname_summary_kk: ", pathname_summary_kk)
  
      ## Calculate on the fly?
      if (!file_test("-f", pathname_summary_kk)) {
        message("overlap_score_summary_grid() ...")
        res <- overlap_score_summary_grid(dataset = dataset, chromosomes = chromosome, bin_sizes = bin_size, rhos = rho, window_size = window_size, nsamples = nsamples, weights = weights, domain_length = domain_length, verbose = verbose)
        message("overlap_score_summary_grid() ... done")
      }
    
      ## Sanity check
      stop_if_not(file_test("-f", pathname_summary_kk))
      summary[[bb,rr]] <- read_rds(pathname_summary_kk)
            
      if (verbose) message(sprintf("Fraction #%d (%s and %s with %s bps on Chr %s) of %d ... done", rr, test_tag, reference_tag, bin_size, chromosome, length(rhos)))
    } ## for (rr ...)
    
    if (verbose) message(sprintf("Bin size #%d (%s) of %d ... done", bb, bin_size_tag, length(bin_sizes)))
  } ## for (bb ...)

#  str(summary)
#  summary <- unlist(summary, recursive = FALSE, use.names = FALSE)
  summary <- do.call(rbind, summary)
  summary <- as_tibble(summary)
  
  if (verbose) mprint(summary)

  saveCache(summary, key = key, dirs = dirs)

  if (verbose) message("read_overlap_score_summary_vs_tad_length() ... done")

  summary
} ## read_overlap_score_summary_vs_tad_length()
