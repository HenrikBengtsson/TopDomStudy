#' Calculate and Summarize TopDom Overlap Scores Across Chromosomes, Bin Sizes, and Fractions
#'
#' @param dataset (character string) The name of the data set.
#'
#' @param chromosomes (character vector) Chromosomes to process.
#'
#' @param bin_sizes (numeric vector) The set of bin sizes (in bps) to process.
#'
#' @param rhos,reference_rhos (numeric vector) The set of fractions (in (0,0.5]) to process.
#'
#' @param window_size (integer) The TopDom windows size.
#'        Argument passed as `window.size` to [TopDom::TopDom()].
#'
#' @param nsamples (integer) Number of random samples for produce.
#'
#' @param weights (character string) A character string specifying how overlap
#' scores across domains should be weighted.
#' Argument passed as is to [overlap_score_summary()].
#'
#' @param domain_length (optional; character string or numeric vector of length two)
#' If specified, controls how to filter out too short or too long TopDom domains.
#' Argument passed as is to [overlap_score_summary()].
#'
#' @param verbose (logical) If `TRUE`, verbose output is produced.
#'
#' @return A three-dimensional character array of pathname names where the first
#' dimension specify `chromosomes`, the second `bin_sizes`, and the third `rhos` (fractions).
#'
#' @section Parallel processing:
#' The \pkg{future} framework is used to parallelize in three layers:
#'  1. across (chromosome, bin size, fraction)
#'  2. [overlap_scores_partitions()]:
#'     1. across a single chromosome (already subsetted above)
#'     2. across `nsamples` random samples
#'
#' An example of a [future::plan()] setup for parallelization on the
#' local machine is:
#' ```r
#'  plan(list(
#'    chr_bin_rho = sequential,   ## across (chr, bin_size, rho)
#'    mono_chr    = sequential,   ## always a single chromosome
#'    samples     = multiprocess  ## across 1:nsamples
#'  ))
#' ```
#' Another is,
#' ```r
#'  plan(list(
#'    chr_bin_rho = multiprocess, ## across (chr, bin_size, rho)
#'    mono_chr    = sequential,   ## always a single chromosome
#'    samples     = sequential    ## across 1:nsamples
#'  ))
#' ```
#' For parallelization on a HPC cluster via a scheduler,
#' ```r
#'  hpc_scheduler <- tweak(future.batchtools::batchtools_torque,
#'                         resources = list(nodes="1:ppn=8", vmem="32gb"))
#'  plan(list(
#'    chr_bin_rho = hpc_scheduler,
#'    mono_chr    = sequential,
#'    samples     = multiprocess
#'  ))
#' ```
#'
#' @importFrom listenv listenv
#' @importFrom future %<-% plan
#' @importFrom future.apply future_lapply
#' @export
overlap_score_summary_grid <- function(dataset, chromosomes, bin_sizes, rhos, reference_rhos = rep(1/2, times = length(rhos)), window_size = 5L, nsamples = 50L, weights = c("by_length", "uniform"), domain_length = NULL, verbose = FALSE) {
  progressor <- import_progressor()

  chromosomes <- as.character(chromosomes)
  stopifnot(is.character(chromosomes), !anyNA(chromosomes))

  stopifnot(is.numeric(rhos), !anyNA(rhos), all(rhos > 0), all(rhos <= 1/2))
  if (is.character(reference_rhos)) {
    reference_rhos <- switch(reference_rhos,
      "50%" = rep(1/2, times = length(rhos)),
      "same" = rhos,
      stop("Unknown value on 'reference_rhos': ", sQuote(reference_rhos))
    )
  }
  stopifnot(is.numeric(reference_rhos), !anyNA(reference_rhos), all(reference_rhos > 0), all(reference_rhos <= 1/2))
  stopifnot(length(reference_rhos) == length(rhos), all(reference_rhos >= rhos))

  stopifnot(length(window_size) == 1L, is.numeric(window_size), !is.na(window_size), window_size >= 1L)
  window_size <- as.integer(window_size)
  window_size_tag <- sprintf("window_size=%d", window_size)

  stopifnot(length(nsamples) == 1L, is.numeric(nsamples), !is.na(nsamples), nsamples >= 1L)
  nsamples_tag <- sprintf("nsamples=%d", nsamples)

  weights <- match.arg(weights)
  weights_tag <- sprintf("weights=%s", weights)
  
  if (is.numeric(domain_length)) {
    stopifnot(length(domain_length) == 2L, !anyNA(domain_length), all(domain_length > 0))
    domain_length_tag <- sprintf("domain_length=%.0f-%.0f", domain_length[1], domain_length[2])
  } else if (is.character(domain_length)) {
    domain_length <- match.arg(domain_length, choices = c("ref_len_iqr"))
    domain_length_tag <- sprintf("domain_length=%s", domain_length)
  } else {
    domain_length_tag <- NULL
  }

  path <- "overlapScoreSummary"
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  dummy <- listenv()
  dim(dummy) <- c(length(chromosomes), length(bin_sizes), length(rhos))
  dimnames(dummy) <- list(chromosomes, bin_sizes, rhos)
  progress <- progressor(prod(dim(dummy)))

  for (cc in seq_along(chromosomes)) {
    chromosome <- chromosomes[cc]
    chromosome_tag <- sprintf("chr=%s", chromosome)

    if (verbose) message(sprintf("Chromosome #%d (%s) of %d ...", cc, chromosome_tag, length(chromosomes)))

    for (bb in seq_along(bin_sizes)) {
      bin_size <- bin_sizes[bb]
      bin_size_tag <- sprintf("bin_size=%.0f", bin_size)

      if (verbose) message(sprintf("Bin size #%d (%s) of %d ...", bb, bin_size_tag, length(bin_sizes)))

      for (rr in seq_along(rhos)) {
        rho <- rhos[rr]
        reference_rho <- reference_rhos[rr]
        test_tag <- sprintf("test=%.5f", rho)
        reference_tag <- sprintf("reference=%.5f", reference_rho)
        if (verbose) message(sprintf("Fraction #%d (%s and %s with %s bps on Chr %s) of %d ...", rr, test_tag, reference_tag, bin_size, chromosome, length(rhos)))

        if (is.character(domain_length) && domain_length == "ref_len_iqr") {
          limits <- extract_domain_length_limits(
            dataset    = dataset,
            chromosome = chromosome,
            bin_size   = bin_size,
            nsamples   = nsamples,
            weights    = weights,
            verbose    = verbose
          )
          if (verbose) mprint(limits)
          stopifnot(nrow(limits) == 1L)
          domain_length <- c(limits[["ref_len_q0.25"]], limits[["ref_len_q0.75"]])
          if (verbose) message("domain_length:")
          if (verbose) mprint(domain_length)
          domain_length_tag <- sprintf("domain_length=%.0f-%.0f", domain_length[1], domain_length[2])
        }
  
        tags <- c(chromosome_tag, "cells_by_half", "avg_score", bin_size_tag, test_tag, reference_tag, window_size_tag, domain_length_tag, weights_tag, nsamples_tag)
        fullname <- paste(c(dataset, tags), collapse = ",")
        pathname_summary_kk <- file.path(path, sprintf("%s.rds", fullname))
        if (verbose) message("pathname_summary_kk: ", pathname_summary_kk)

        progress(message = paste(c(chromosome_tag, bin_size_tag, test_tag, reference_tag), collapse=", "))

        ## Already processed?
        if (file_test("-f", pathname_summary_kk)) {
          dummy[[cc, bb, rr]] <- pathname_summary_kk
          if (verbose) message(sprintf("Fraction #%d (%s and %s with %s bps on Chr %s) of %d ... already done", rr, test_tag, reference_tag, bin_size, chromosome, length(rhos)))
          next
        }

        dummy[[cc, bb, rr]] %<-% {
          if (verbose) message("Remaining future::plan():")
          if (verbose) mprint(plan("list"))

          filename <- sprintf("%s,unique,chr=%s.rds", dataset, chromosome)
          pathname <- system.file("compiledData", filename, package = "TopDomStudy", mustWork = TRUE)
          if (verbose) message(sprintf("Reads (%s):", pathname))
          reads <- read_rds(pathname)
          if (verbose) mprint(reads)

          if (verbose) message("overlap_scores_partitions() ...")
          res <- overlap_scores_partitions(reads = reads, dataset = sprintf("%s,unique", dataset), bin_size = bin_size,
                                           partition_by = "cells_by_half", min_cell_size = 2L, window_size = window_size, rho = rho, reference_rho = reference_rho,
                                           nsamples = nsamples, chrs = chromosome, seed = 0xBEEF, mainseed = 0xBEEF, verbose = verbose)
          if (verbose) mstr(res)
          if (verbose) message("overlap_scores_partitions() ... done")

          ## Summary of overlap scores and reference domain lengths
          if (verbose) message("Summary of overlap scores and reference domain lengths ...")
          res_chr <- res[[chromosome]]
          summary_kk <- future_lapply(res_chr, FUN = function(pathname) {
            oss <- read_rds(pathname)
            ## Drop failed TopDom fits and possibly skip this sample?
            failed <- unlist(lapply(oss, FUN = inherits, "try-error"))
            if (any(failed)) {
              oss <- oss[!failed]
              if (length(oss) < 2) return(NULL)
            }
            z <- overlap_score_summary(oss, weights = weights, domain_length = domain_length)
            ## Sanity checks
            stopifnot(all(z$fraction == rho), all(z$ref_fraction == reference_rho))
            oss <- failed <- NULL

            ## Locate TopDom fit results
            set <- basename(dirname(pathname))
            path_td <- file.path("topdomData", set)
            stop_if_not(file_test("-d", path_td))
            filename_td <- basename(pathname)
            pathname_td <- file.path(path_td, filename_td)
            stop_if_not(file_test("-f", pathname_td))
            td <- read_rds(pathname_td)

            ref <- grep("^reference", names(td))
            stopifnot(length(ref) == 1L)
            sizes <- td[[ref]]$domain$size
            ## Filter by domain lengths?
            if (!is.null(domain_length)) {
              keep <- (domain_length[1] <= sizes & sizes <= domain_length[2])
              sizes <- sizes[keep]
            }
            probs <- c(0.00, 0.05, 0.25, 0.50, 0.75, 0.95, 1.00)
            qsizes <- quantile(sizes, probs = probs, na.rm = TRUE)
            names(qsizes) <- sprintf("ref_len_q%0.2f", probs)
            z <- cbind(z, as.list(qsizes))
    
            sizes <- td[-ref][[1]]$domain$size
            ## Filter by domain lengths?
            if (!is.null(domain_length)) {
              keep <- (domain_length[1] <= sizes & sizes <= domain_length[2])
              sizes <- sizes[keep]
              }
            probs <- c(0.00, 0.05, 0.25, 0.50, 0.75, 0.95, 1.00)
            qsizes <- quantile(sizes, probs = probs, na.rm = TRUE)
            names(qsizes) <- sprintf("test_len_q%0.2f", probs)
            z <- cbind(z, as.list(qsizes))
   
            z
          })
          summary_kk <- do.call(rbind, summary_kk)
          rownames(summary_kk) <- NULL
          if (verbose) message("Summary of overlap scores and reference domain lengths ... done")

          ## Save intermediate results to file
          save_rds(summary_kk, pathname_summary_kk)
          if (verbose) message("Saved pathname_summary_kk: ", pathname_summary_kk)

          pathname_summary_kk
        } %label% sprintf("osg_%s-%s-%s", bin_size, rho, chromosome)

        if (verbose) message(sprintf("Fraction #%d (%s and %s with %s bps on Chr %s) of %d ... done", rr, test_tag, reference_tag, bin_size, chromosome, length(rhos)))
      } ## for (rr ...)

      if (verbose) message(sprintf("Bin size #%d (%s bps on Chr %s) of %d ... done", bb, bin_size, chromosome, length(bin_sizes)))
    } ## for (bb ...)
    
    if (verbose) message(sprintf("Chromosome #%d (%s) of %d ... done", cc, chromosome_tag, length(chromosomes)))
  } ## for (cc ...)
  
  if (verbose) message("All tasks submitted as futures")
  
  ## Resolve
  dummy <- as.list(dummy)

  ## Coerce to a character array
  pathnames <- unlist(dummy)
  dim(pathnames) <- dim(dummy)
  dimnames(pathnames) <- dimnames(dummy)

  attr(pathnames, "domain_length") <- domain_length

  pathnames
}
