#' Calculate and Summarize TopDom Overlap Scores as Function of Sample Fraction
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
#' The \pkg{future} framework is used to parallelize [TopDom::TopDom()] in 2+3 layers:
#'  1. across (chromosome, bin_size) (arguments `chromosomes` and `bin_sizes`)
#'  2. across all fractions (argument `rhos`)
#'  3. [overlap_scores_partitions()] layers:
#'    a. chromosomes (here a single one)
#'    b. random samples (argument `nsamples`)
#'    c. partions per sample (typically only two)
#'
#' @importFrom listenv listenv
#' @importFrom future %<-% plan
#' @importFrom future.apply future_lapply
#' @export
overlap_score_summary_vs_fraction <- function(dataset, chromosomes, bin_sizes, rhos, window_size = 5L, nsamples = 50L, weights = c("by_length", "uniform"), domain_length = NULL, figures = TRUE, verbose = FALSE) {
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

  if (figures) {
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

      dummy[[cc, bb, 1L]] %<-% {
        ## WORKAROUND: Below is a case of "y %<-% { if (reset) x <- 0; x + 1 }"
	## Help future identify these as globals:
        domain_length_tag
	
        message("Remaining future::plan():")
        mprint(plan("list"))

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

        pathnames_summary <- character(length(rhos))
        summary <- listenv()
        for (rr in seq_along(rhos)) {
          rho <- rhos[rr]
          rho_tag <- sprintf("fraction=%.3f", rho)
          message(sprintf("Fraction #%d (%s with %s bps on Chr %s) of %d ...", rr, rho_tag, bin_size, chromosome, length(rhos)))

          tags <- c(chromosome_tag, "cells_by_half", "avg_score", bin_size_tag, rho_tag, window_size_tag, domain_length_tag, weights_tag, nsamples_tag)
          fullname <- paste(c(dataset, tags), collapse = ",")
          pathname_summary_kk <- file.path(path, sprintf("%s.rds", fullname))
	  pathnames_summary[rr] <- pathname_summary_kk
          message("pathname_summary_kk: ", pathname_summary_kk)

          ## Already processed?
          if (file_test("-f", pathname_summary_kk)) {
            summary[[rr]] <- read_rds(pathname_summary_kk)
            
            message(sprintf("Fraction #%d (%s with %s bps on Chr %s) of %d ... already done", rr, rho_tag, bin_size, chromosome, length(rhos)))
            next
          }

          summary[[rr]] %<-% {
            message("Remaining future::plan():")
            mprint(plan("list"))

            message(sprintf("Fraction #%d (%s with %s bps on Chr %s) of %d ... already done", rr, rho_tag, bin_size, chromosome, length(rhos)))

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
            summary_kk %<-% future_lapply(res_chr, FUN = function(pathname) {
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

            summary_kk
          } %label% sprintf("%s-%s-%s", chromosome, bin_size, rho)

          message(sprintf("Fraction #%d (%s with %s bps on Chr %s) of %d ... done", rr, rho_tag, bin_size, chromosome, length(rhos)))
        } ## for (rr ...)

        ## Resolve futures
        summary <- as.list(summary)

        summary <- do.call(rbind, summary)
        mprint(summary)

        if (figures) {
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
        } ## if (figures)

        pathnames_summary
      } %label% sprintf("%s-%s", chromosome, bin_size)
      
      message(sprintf("Bin size #%d (%s bps on Chr %s) of %d ... done", bb, bin_size, chromosome, length(bin_sizes)))
    } ## for (bb ...)
    
    message(sprintf("Chromosome #%d (%s) of %d ... done", cc, chromosome_tag, length(chromosomes)))
  } ## for (cc ...)
  
  message("All tasks submitted as futures")
  
  ## Resolve
  dummy <- as.list(dummy)

  ## AD HOC: distribute pathnames
  for (aa in seq_along(dim(dummy)[1])) {
    for (bb in seq_along(dim(dummy)[2])) {
      pathnames <- dummy[aa,bb,1]
      dummy[aa,bb,] <- pathnames
    }
  }

  ## Coerce to a character array
  pathnames <- unlist(dummy)
  dim(pathnames) <- dim(dummy)
  dimnames(pathnames) <- dimnames(dummy)

  pathnames
}
