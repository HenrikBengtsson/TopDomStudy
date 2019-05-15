#' Calculate and Summarize TopDom Overlap Scores as Function of Sample Fraction
#'
#' @inheritParams overlap_score_summary_grid
#'
#' @return A three-dimensional character array of pathname names where the
#' first dimension specify `chromosomes`, the second `bin_sizes`, and
#' the third 'rhos'.
#'
#' @details
#' PNG images are written to the \file{figures/} folder (created if missing).
#'
#' @seealso
#' Internal, [overlap_score_summary_grid()] is used to calculate overlap
#' scores over (chromosome, bin_size, rho).
#'
#' @importFrom ggplot2 aes aes_string geom_boxplot geom_jitter ggplot ggsave ggtitle stat_summary xlab ylab ylim
#' @importFrom utils file_test
#' @export
overlap_score_summary_vs_fraction <- function(dataset, chromosomes, bin_sizes, rhos, window_size = 5L, nsamples = 50L, weights = c("by_length", "uniform"), domain_length = NULL, verbose = FALSE) {
  weights <- match.arg(weights)
  weights_tag <- sprintf("weights=%s", weights)
  
  pathnames <- overlap_score_summary_grid(dataset, chromosomes = chromosomes, bin_sizes = bin_sizes, rhos = rhos, window_size = window_size, nsamples = nsamples, weights = weights, domain_length = domain_length, verbose = verbose)

  nsamples_tag <- sprintf("nsamples=%d", nsamples)

  window_size <- as.integer(window_size)
  window_size_tag <- sprintf("window_size=%d", window_size)

  domain_length <- attr(pathnames, "domain_length")
  if (!is.null(domain_length)) {
    stop_if_not(is.numeric(domain_length), length(domain_length) == 2L, !anyNA(domain_length), all(domain_length > 0))
    domain_length_tag <- sprintf("domain_length=%.0f-%.0f", domain_length[1], domain_length[2])
  } else {
    domain_length_tag <- NULL
  }

  path <- "overlapScoreSummary"
  stop_if_not(file_test("-d", path))

  message("Plotting ...")
  
  for (cc in seq_along(chromosomes)) {
    chromosome <- chromosomes[cc]
    chromosome_tag <- sprintf("chr=%s", chromosome)

    message(sprintf("Chromosome #%d (%s) of %d ...", cc, chromosome_tag, length(chromosomes)))

    for (bb in seq_along(bin_sizes)) {
      bin_size <- bin_sizes[bb]
      bin_size_tag <- sprintf("bin_size=%.0f", bin_size)
      message(sprintf("Bin size #%d (%s) of %d ...", bb, bin_size_tag, length(bin_sizes)))

      summary <- list()
      for (rr in seq_along(rhos)) {
        rho <- rhos[rr]
        rho_tag <- sprintf("fraction=%.3f", rho)
        message(sprintf("Fraction #%d (%s with %s bps on Chr %s) of %d ...", rr, rho_tag, bin_size, chromosome, length(rhos)))

        tags <- c(chromosome_tag, "cells_by_half", "avg_score", bin_size_tag, rho_tag, window_size_tag, domain_length_tag, weights_tag, nsamples_tag)
        fullname <- paste(c(dataset, tags), collapse = ",")
        pathname_summary_kk <- file.path(path, sprintf("%s.rds", fullname))
        message("pathname_summary_kk: ", pathname_summary_kk)
        message("pathnames[cc,bb,rr]: ", pathnames[cc,bb,rr])
        ## Sanity check
        stop_if_not(identical(pathname_summary_kk, pathnames[cc,bb,rr]))
        stop_if_not(file_test("-f", pathname_summary_kk))

        summary[[rr]] <- read_rds(pathname_summary_kk)
            
        message(sprintf("Fraction #%d (%s with %s bps on Chr %s) of %d ... already done", rr, rho_tag, bin_size, chromosome, length(rhos)))
      } ## for (rr ...)

      summary <- do.call(rbind, summary)
      mprint(summary)
    
      dw <- diff(range(rhos)) / length(rhos)

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

        gg <- ggplot(summary, aes_string(x = "fraction", y = signal))

        gg <- gg + geom_boxplot(aes(group = as.factor(fraction)), width = 0.2*dw)
        gg <- gg + geom_jitter(height = 0, width = 0.05*dw, size = 0.7, colour = "darkgray")

        gg <- gg + stat_summary(aes_string(y = signal, group = 1L),
                                fun.y = function(x) mean(x, trim = 0.10),
                                geom = "line", size = 2L, group = 1L)

        params <- c(sprintf("estimator: %s", signal_label),
                    sprintf("weights: %s", weights),
                    sprintf("domains: %.0f-%.0f", domain_length[1], domain_length[2]))
        subtitle <- sprintf("chromosome %s, bin size=%d, window size=%d (%d samples) [%s]",
                            chromosome, bin_size, window_size, nsamples, paste(params, collapse = "; "))

        gg <- gg + ggtitle(dataset, subtitle = subtitle)
        gg <- gg + xlab("fraction")
        if (signal_label %in% names(length_signals)) {
          gg <- gg + ylab("domain length (bps)")
          gg <- gg + ylim(0, 2e6)
        } else {
          gg <- gg + ylab("average overlap score")
          gg <- gg + ylim(0, 1)
        }

        signal <- gsub("`50%`", "median", signal)
        tags <- sprintf("%s,chr=%s,%s,avg_score-vs-fraction,bin_size=%d,window_size=%d,nsamples=%d,signal=%s,weights=%s", dataset, chromosome, "cells_by_half", bin_size, window_size, nsamples, signal, weights)
        filename <- sprintf("%s.png", paste(c(tags, domain_length_tag), collapse = ","))
        dir.create("figures", recursive = TRUE, showWarnings = FALSE)
        ggsave(gg, filename=file.path("figures", filename))
      } ## for (signal ...)
      
      message(sprintf("Bin size #%d (%s) of %d ... done", bb, bin_size_tag, length(bin_sizes)))
    } ## for (bb ...)

    message(sprintf("Chromosome #%d (%s) of %d ... done", cc, chromosome_tag, length(chromosomes)))
  } ## for (cc ...)

  message("Plotting ... done")

  pathnames
}
