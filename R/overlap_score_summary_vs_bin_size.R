#' Calculate and Summarize TopDom Overlap Scores as Function of Bin Size
#'
#' @param dataset (character string) The name of the data set.
#'
#' @param chromosomes (character vector) Chromosomes to process.
#'
#' @param bin_sizes (numeric vector) The set of bin sizes (in bps) to process.
#'
#' @param rhos (numeric vector) The set of fractions (in (0,0.5]) to process.
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
#' @param figures (logical) If TRUE, PNG images are written to the \file{figures/}
#' folder (created if missing).  This requires the \pkg{ggplot2} packages.
#'
#' @param verbose (logical) If TRUE, verbose output is produced.
#'
#' @return A three-dimensional character array of pathname names where the first
#' dimension specify `chromosomes`, the second `bin_sizes`, and the third 'rhos'.
#'
#' @section Parallel processing:
#' The \pkg{future} framework is used to parallelize [TopDom::TopDom()] in 1+2 layers:
#'  1. across (chromosome, bin size fraction) (arguments `chromosomes`, `bin_sizes`, and `rhos`)
#'  2. [overlap_scores_partitions()] layers:
#'    i. random samples (argument `nsamples`)
#'    ii. partions per sample (typically only two)
#'
#' @importFrom listenv listenv
#' @importFrom future %<-% plan
#' @importFrom future.apply future_lapply
#' @export
overlap_score_summary_vs_bin_size <- function(dataset, chromosomes, bin_sizes, rhos, window_size = 5L, nsamples = 50L, weights = c("by_length", "uniform"), domain_length = NULL, figures = TRUE, verbose = FALSE) {
  stopifnot(length(window_size) == 1L, is.numeric(window_size), !is.na(window_size), window_size >= 1L)
  window_size <- as.integer(window_size)
  window_size_tag <- sprintf("window_size=%d", window_size)

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
  for (cc in seq_along(chromosomes)) {
    chromosome <- chromosomes[cc]
    chromosome_tag <- sprintf("chr=%s", chromosome)

    message(sprintf("Chromosome #%d (%s) of %d ...", cc, chromosome_tag, length(chromosomes)))

    for (bb in seq_along(bin_sizes)) {
      bin_size <- bin_sizes[bb]
      bin_size_tag <- sprintf("bin_size=%.0f", bin_size)

      message(sprintf("Bin size #%d (%s) of %d ...", bb, bin_size_tag, length(bin_sizes)))

      for (rr in seq_along(rhos)) {
        rho <- rhos[rr]
        rho_tag <- sprintf("fraction=%.3f", rho)
        message(sprintf("Fraction #%d (%s with %s bps on Chr %s) of %d ...", rr, rho_tag, bin_size, chromosome, length(rhos)))

        if (is.character(domain_length) && domain_length == "ref_len_iqr") {
          limits <- extract_domain_length_limits(
            dataset    = dataset,
            chromosome = chromosome,
            bin_size   = bin_size,
            nsamples   = nsamples,
            weights    = weights,
            verbose    = verbose
          )
          mprint(limits)
          stopifnot(nrow(limits) == 1L)
          domain_length <- c(limits[["ref_len_q0.25"]], limits[["ref_len_q0.75"]])
          message("domain_length:")
          mprint(domain_length)
          domain_length_tag <- sprintf("domain_length=%.0f-%.0f", domain_length[1], domain_length[2])
        }
  
        tags <- c(chromosome_tag, "cells_by_half", "avg_score", bin_size_tag, rho_tag, window_size_tag, domain_length_tag, weights_tag, nsamples_tag)
        fullname <- paste(c(dataset, tags), collapse = ",")
        pathname_summary_kk <- file.path(path, sprintf("%s.rds", fullname))
        message("pathname_summary_kk: ", pathname_summary_kk)

        ## Already processed?
        if (file_test("-f", pathname_summary_kk)) {
          dummy[[cc, bb, rr]] <- pathname_summary_kk
          message(sprintf("Fraction #%d (%s with %s bps on Chr %s) of %d ... already done", rr, rho_tag, bin_size, chromosome, length(rhos)))
          next
        }

        dummy[[cc, bb, rr]] %<-% {
          message("Remaining future::plan():")
          mprint(plan("list"))

          filename <- sprintf("%s,unique,chr=%s.rds", dataset, chromosome)
          pathname <- system.file("compiledData", filename, package = "TopDomStudy", mustWork = TRUE)
          message(sprintf("Reads (%s):", pathname))
          reads <- read_rds(pathname)
          mprint(reads)

          message("overlap_scores_partitions() ...")
          res <- overlap_scores_partitions(reads = reads, dataset = sprintf("%s,unique", dataset), bin_size = bin_size,
                                           partition_by = "cells_by_half", min_cell_size = 2L, window_size = window_size, rho = rho,
                                           nsamples = nsamples, chrs = chromosome, seed = 0xBEEF, mainseed = 0xBEEF, verbose = verbose)
          mstr(res)
          message("overlap_scores_partitions() ... done")

          ## Summary of overlap scores and reference domain lengths
          message("Summary of overlap scores and reference domain lengths ...")
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
            oss <- failed <- NULL

            pathname_td <- gsub("[.]rds$", ",topdom.rds", pathname)
            td <- read_rds(pathname_td)

            ref <- which(names(td) == "reference")
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
          summary_kk <- cbind(summary_kk, fraction = rho)
          message("Summary of overlap scores and reference domain lengths ... done")

          ## Save intermediate results to file
          saveRDS(summary_kk, file = pathname_summary_kk)
          message("Saved pathname_summary_kk: ", pathname_summary_kk)

          pathname_summary_kk
        } %label% sprintf("%s-%s-%s", chromosome, bin_size, rho)

        message(sprintf("Fraction #%d (%s with %s bps on Chr %s) of %d ...", rr, rho_tag, bin_size, chromosome, length(rhos)))
      } ## for (rr ...)

      message(sprintf("Bin size #%d (%s bps on Chr %s) of %d ... done", bb, bin_size, chromosome, length(bin_sizes)))
    } ## for (bb ...)
    
    message(sprintf("Chromosome #%d (%s) of %d ... done", cc, chromosome_tag, length(chromosomes)))
  } ## for (cc ...)
  
  message("All tasks submitted as futures")
  
  ## Resolve
  dummy <- as.list(dummy)

  ## Coerce to a character array
  pathnames <- unlist(dummy)
  dim(pathnames) <- dim(dummy)
  dimnames(pathnames) <- dimnames(dummy)

  if (figures) {
    message("Plotting ...")
    
    aes <- ggplot2::aes
    aes_string <- ggplot2::aes_string
    geom_boxplot <- ggplot2::geom_boxplot
    geom_jitter <- ggplot2::geom_jitter
    ggplot <- ggplot2::ggplot
    ggsave <- ggplot2::ggsave
    ggtitle <- ggplot2::ggtitle
    stat_summary <- ggplot2::stat_summary
    xlab <- ggplot2::xlab
    ylab <- ggplot2::ylab
    ylim <- ggplot2::ylim
    
    for (cc in seq_along(chromosomes)) {
      chromosome <- chromosomes[cc]
      chromosome_tag <- sprintf("chr=%s", chromosome)

      message(sprintf("Chromosome #%d (%s) of %d ...", cc, chromosome_tag, length(chromosomes)))

      for (rr in seq_along(rhos)) {
        rho <- rhos[rr]
        rho_tag <- sprintf("fraction=%.3f", rho)

        message(sprintf("Fraction #%d (%g on Chr %s) of %d ... done", rr, rho, chromosome, length(rhos)))

        summary <- list()
        
        for (bb in seq_along(bin_sizes)) {
          bin_size <- bin_sizes[bb]
          bin_size_tag <- sprintf("bin_size=%.0f", bin_size)
          message(sprintf("Bin size #%d (%s bps with %g on Chr %s) of %d ...", bb, bin_size, rho, chromosome, length(bin_sizes)))

          tags <- c(chromosome_tag, "cells_by_half", "avg_score", bin_size_tag, rho_tag, window_size_tag, domain_length_tag, weights_tag, nsamples_tag)
          fullname <- paste(c(dataset, tags), collapse = ",")
          pathname_summary_kk <- file.path(path, sprintf("%s.rds", fullname))
          message("pathname_summary_kk: ", pathname_summary_kk)
          message("pathnames[cc,bb,rr]: ", pathnames[cc,bb,rr])
          ## Sanity check
          stop_if_not(identical(pathname_summary_kk, pathnames[cc,bb,rr]))
          stop_if_not(file_test("-f", pathname_summary_kk))
          
          summary[[bb]] <- read_rds(pathname_summary_kk)
            
          message(sprintf("Bin size #%d (%s bps with %g on Chr %s) of %d ... already done", bb, bin_size, rho, chromosome, length(bin_sizes)))
        } ## for (bb ...)

        summary <- do.call(rbind, summary)
        mprint(summary)

        dw <- diff(range(bin_sizes)) / length(bin_sizes)

        length_signals <- c(
          "reference Q25 length"    = "ref_len_q0.25",
          "reference median length" = "ref_len_q0.50",
          "reference Q75 length"    = "ref_len_q0.75",
          "test Q25 length"         = "test_len_q0.25",
          "test median length"      = "test_len_q0.50",
          "test Q75 length"         = "test_len_q0.75"
        )
        signals <- c(mean = "mean", median = "`50%`", length_signals)

        fraction <- NULL; rm(list = "fraction") ## To please R CMD check

        for (signal_label in names(signals)) {
          signal <- signals[[signal_label]]

          gg <- ggplot(summary, aes_string(x = "bin_size", y = signal))

          gg <- gg + geom_boxplot(aes(group = as.factor(bin_size)), width = 0.2*dw)
          gg <- gg + geom_jitter(height = 0, width = 0.05*dw, size = 0.7, colour = "darkgray")

          gg <- gg + stat_summary(aes_string(y = signal, group = 1L),
                                  fun.y = function(x) mean(x, trim = 0.10),
                                  geom = "line", size = 2L, group = 1L)

          params <- c(sprintf("estimator: %s", signal_label),
                      sprintf("weights: %s", weights),
                      sprintf("domains: %.0f-%.0f", domain_length[1], domain_length[2]))
          subtitle <- sprintf("chromosome %s, fraction=%.3f (%d samples) [%s]",
                              chromosome, rho, nsamples, paste(params, collapse = "; "))
  
          gg <- gg + ggtitle(dataset, subtitle = subtitle)
          gg <- gg + xlab("bin size (bps)")
          if (signal_label %in% names(length_signals)) {
            gg <- gg + ylab("domain length (bps)")
            gg <- gg + ylim(0, 2e6)
          } else {
            gg <- gg + ylab("average overlap score")
            gg <- gg + ylim(0, 1)
          }
  
          signal <- gsub("`50%`", "median", signal)
          tags <- sprintf("%s,chr=%s,%s,avg_score-vs-bin_size,fraction=%.3f,window_size=%d,nsamples=%d,signal=%s,weights=%s", dataset, chromosome, "cells_by_half", rho, window_size, nsamples, signal, weights)
          filename <- sprintf("%s.png", paste(c(tags, domain_length_tag), collapse = ","))
          dir.create("figures", recursive = TRUE, showWarnings = FALSE)
          ggsave(gg, filename=file.path("figures", filename))
        } ## for (signal ...)
          
        message(sprintf("Fraction #%d (%g on Chr %s) of %d ... done", rr, rho, chromosome, length(rhos)))
      } ## for (rr ...)
      
      message(sprintf("Chromosome #%d (%s) of %d ... done", cc, chromosome_tag, length(chromosomes)))
    } ## for (cc ...)
  
    message("Plotting ... done")
  } ## if (figures)

  pathnames
}
