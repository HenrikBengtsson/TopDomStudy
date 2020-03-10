#' Plots TopDom Overlap Scores and Lengths as a Function of Sample Fraction or Bin Size
#'
#' @inheritParams overlap_score_summary_grid
#'
#' @param chromosome (character string) The chromosome to be plotted.
#'
#' @param bin_size,rho (numeric) The bin size (bps) or sample fraction (in (0,0.5]) to be plotted.
#'
#' @param signals (character vector of length two) The two signals to be
#' plotted in the top and bottom panels.
#'
#' @param rho_lim,bin_size_lim The range of the fraction or bin size to be displayed.
#'
#' @param length_lim The range of the domain lenghts to be displayed.
#'
#' @param rel_heights (numeric vector of length two) The relative height
#' of the two figures.  Passed to [cowplot::plot_grid] as is.
#'
#' @param labels (character vector).  Enable or disable some labels.
#'
#' @param line_col The color of the profile curve.
#'
#' @param fig_path If non-NULL, a PNG image is written to this path.
#'
#' @param skip If TRUE and the PNG file already exists, then it is not replotted.
#'
#' @param \ldots Not used.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @example incl/gg_overlap_score_summary_vs_nnn.R
#'
#' @seealso
#' Internal, [overlap_score_summary_grid()] is used to calculate overlap
#' scores over (chromosome, bin_size, rho).
#'
#' @importFrom grid unit
#' @importFrom ggplot2 aes aes_string coord_cartesian element_blank
#'             geom_boxplot geom_jitter geom_text ggplot ggsave ggtitle
#'             element_text scale_x_continuous stat_summary theme xlab ylab ylim
#' @importFrom scales percent_format
#' @importFrom cowplot plot_grid
#' @importFrom utils file_test
#' @export
gg_overlap_score_summary_vs_fraction <- function(dataset, chromosome, bin_size, rhos, window_size = 5L, nsamples = 50L, weights = c("by_length", "uniform"), domain_length = NULL, signals = c("mean", "test_len_q0.50"), rho_lim = c(0, 1/2), length_lim = c(0, 2e3), rel_heights = c(4,1), labels = c("chromosome"), line_col = "black", fig_path = "figures", skip = FALSE, ..., verbose = TRUE) {
  avg_signals <- c(mean = "mean", median = "`50%`", "q25" = "`25%`", "q75" = "`75%`")
  length_signals <- c(
    "reference Q25 length"    = "ref_len_q0.25",
    "reference median length" = "ref_len_q0.50",
    "reference Q75 length"    = "ref_len_q0.75",
    "test Q25 length"         = "test_len_q0.25",
    "test median length"      = "test_len_q0.50",
    "test Q75 length"         = "test_len_q0.75"
  )
  known_signals <- c(avg_signals, length_signals)
  stopifnot(length(signals) == 2L, all(signals %in% known_signals))
  signal_labels <- character(2L)
  signal_labels[1] <- names(known_signals)[signals[1] == known_signals]
  signal_labels[2] <- signals[2]
  stop_if_not(length(signal_labels) == 2L)

  if (verbose) {
    message("gg_overlap_score_summary_vs_fraction() ...")
    message("- chromosome: ", chromosome)
    message("- bin_size: ", bin_size)
    message("- window_size: ", window_size)
    message("- weights: ", weights)
    message("- nsamples: ", nsamples)
    message("- signals: ", paste(sQuote(signals), collapse = ", "))
    message("- signal labels: ", paste(sQuote(signal_labels), collapse = ", "))
  }

  if (!is.null(fig_path) && !file_test("-d", fig_path)) {
    dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)
    stop_if_not(file_test("-d", fig_path))
  }

  ## Tags
  chromosome_tag <- sprintf("chr=%s", chromosome)
  bin_size_tag <- sprintf("bin_size=%.0f", bin_size)
  window_size_tag <- sprintf("window_size=%d", window_size)
  if (!is.null(domain_length)) {
    stop_if_not(is.numeric(domain_length), length(domain_length) == 2L, !anyNA(domain_length), all(domain_length > 0))
    domain_length_tag <- sprintf("domain_length=%.0f-%.0f", domain_length[1], domain_length[2])
  } else {
    domain_length_tag <- NULL
  }
  weights_tag <- sprintf("weights=%s", weights)
  nsamples_tag <- sprintf("nsamples=%d", nsamples)

  fig_pathname <- NULL
  if (!is.null(fig_path)) {
    signals_tag <- sprintf("signals=%s", paste(signal_labels, collapse = "-"))
  
    tags <- c(chromosome_tag, "cells_by_half", "avg_score-vs-fraction", bin_size_tag, window_size_tag, nsamples_tag, signals_tag, weights_tag, domain_length_tag)
    fullname <- paste(c(dataset, tags), collapse = ",")
    filename <- sprintf("%s.png", fullname)
    fig_pathname <- file.path(fig_path, filename)

    ## Nothing to do?
    if (skip && file_test("-f", fig_pathname)) {
      res <- list()
      attr(res, "fig_pathname") <- normalizePath(fig_pathname)
      return(res)
    }
  }

  summary <- read_overlap_score_summary_vs_fraction(dataset, chromosome = chromosome, bin_size = bin_size, rhos = rhos, window_size = window_size, nsamples = nsamples, weights = weights, domain_length = domain_length, ..., verbose = verbose)

  fraction <- NULL; rm(list = "fraction") ## To please R CMD check

  dw <- diff(range(rhos)) / length(rhos)

  signal_label <- names(known_signals)[signals[1] == known_signals]
  params <- c(sprintf("estimator: %s", signal_label),
              sprintf("weights: %s", weights),
              sprintf("domains: %.0f-%.0f", domain_length[1], domain_length[2]))
  subtitle <- sprintf("chromosome %s, bin size=%d, window size=%d (%d samples) [%s]",
                       chromosome, bin_size, window_size, nsamples, paste(params, collapse = "; "))

  ggs <- list()
  for (ss in seq_along(signals)) {
    signal <- signals[ss]
    if (signal %in% length_signals) {
      summary[[signal]] <- summary[[signal]] / 1000
    }
    gg <- ggplot(summary, aes_string(x = "fraction", y = signal))
    gg <- gg + geom_boxplot(aes(group = as.factor(fraction)), width = 0.2*dw)
    gg <- gg + geom_jitter(height = 0, width = 0.05*dw, size = 0.7, colour = "darkgray")
    gg <- gg + stat_summary(aes_string(y = signal, group = 1L),
                            fun.y = function(x) mean(x, trim = 0.10),
                            geom = "line",
                            colour = line_col, size = 2L, group = 1L)
    
    if (ss == 1L) {                              
      gg <- gg + ggtitle(label = NULL, subtitle = subtitle) + theme(plot.subtitle = element_text(size = 8))
      if ("chromosome" %in% labels) gg <- gg + geom_text_top_left(sprintf("Chr %s", chromosome), size=5)
      gg <- gg + geom_text_top_right(sprintf("%.0f kbp", bin_size/1000), size=5)
      gg <- gg + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
      gg <- gg + xlab("")
    } else {
      gg <- gg + scale_x_continuous(labels = percent_format(accuracy = 2L))
      gg <- gg + xlab("sample fraction")
    }
    if (signal %in% length_signals) {
      ylim <- length_lim
      gg <- gg + ylab("domain length (kbp)")
      gg <- gg + theme(plot.margin = unit(c(5, 0, 5, 5), "pt"))
    } else {
      ylim <- c(0, 1)
      gg <- gg + ylab("average overlap score")
      gg <- gg + theme(plot.margin = unit(c(5, 5, -20, 5), "pt"))
    }
    gg <- gg + theme(axis.text = element_text(size = 13))
    gg <- gg + theme(axis.title = element_text(size = 13))

    ## xlim and ylim must be set at the same time
    gg <- gg + coord_cartesian(xlim = rho_lim, ylim = ylim)

    ggs[[signal]] <- gg
  }

  res <- plot_grid(ggs[[1]], ggs[[2]], ncol = 1L, rel_heights = rel_heights, align = "v")

  if (!is.null(fig_pathname)) {
    if (verbose) suppressMessages <- identity
    suppressMessages(ggsave(res, filename = fig_pathname,, scale = getOption("TopDomStudy.ggsave.scale", 0.80), dpi = getOption("TopDomStudy.ggsave.dpi", 360)))
    attr(res, "fig_pathname") <- normalizePath(fig_pathname)
    if (verbose) message(" - Figure written: ", fig_pathname)
  }

  if (verbose) message("gg_overlap_score_summary_vs_fraction() ... done")

  res
} ## gg_overlap_score_summary_vs_fraction()








#' @rdname gg_overlap_score_summary_vs_fraction
#'
#' @importFrom grid unit
#' @importFrom ggplot2 aes aes_string coord_cartesian element_blank
#'             geom_boxplot geom_jitter geom_text ggplot ggsave ggtitle
#'             element_text scale_x_continuous stat_summary theme xlab ylab ylim
#' @importFrom scales percent_format
#' @importFrom cowplot plot_grid
#' @importFrom utils file_test
#'
#' @export
gg_overlap_score_summary_vs_bin_size <- function(dataset, chromosome, bin_sizes, rho, window_size = 5L, nsamples = 50L, weights = c("by_length", "uniform"), domain_length = NULL, signals = c("mean", "test_len_q0.50"), bin_size_lim = c(0, max(bin_sizes)), length_lim = c(0, 2e3), rel_heights = c(4,1), labels = c("chromosome"), line_col = "black", fig_path = "figures", skip = FALSE, ..., verbose = TRUE) {
  avg_signals <- c(mean = "mean", median = "`50%`", "q25" = "`25%`", "q75" = "`75%`")
  length_signals <- c(
    "reference Q25 length"    = "ref_len_q0.25",
    "reference median length" = "ref_len_q0.50",
    "reference Q75 length"    = "ref_len_q0.75",
    "test Q25 length"         = "test_len_q0.25",
    "test median length"      = "test_len_q0.50",
    "test Q75 length"         = "test_len_q0.75"
  )
  known_signals <- c(avg_signals, length_signals)
  stopifnot(length(signals) == 2L, all(signals %in% known_signals))
  signal_labels <- character(2L)
  signal_labels[1] <- names(known_signals)[signals[1] == known_signals]
  signal_labels[2] <- signals[2]
  stop_if_not(length(signal_labels) == 2L)

  if (verbose) {
    message("gg_overlap_score_summary_vs_bin_size() ...")
    message("- chromosome: ", chromosome)
    message("- rho: ", rho)
    message("- window_size: ", window_size)
    message("- weights: ", weights)
    message("- nsamples: ", nsamples)
    message("- signals: ", paste(sQuote(signals), collapse = ", "))
    message("- signal labels: ", paste(sQuote(signal_labels), collapse = ", "))
  }

  if (!is.null(fig_path) && !file_test("-d", fig_path)) {
    dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)
    stop_if_not(file_test("-d", fig_path))
  }

  ## Tags
  chromosome_tag <- sprintf("chr=%s", chromosome)
  rho_tag <- sprintf("test=%.3f", rho)
  window_size_tag <- sprintf("window_size=%d", window_size)
  if (!is.null(domain_length)) {
    stop_if_not(is.numeric(domain_length), length(domain_length) == 2L, !anyNA(domain_length), all(domain_length > 0))
    domain_length_tag <- sprintf("domain_length=%.0f-%.0f", domain_length[1], domain_length[2])
  } else {
    domain_length_tag <- NULL
  }
  weights_tag <- sprintf("weights=%s", weights)
  nsamples_tag <- sprintf("nsamples=%d", nsamples)

  fig_pathname <- NULL
  if (!is.null(fig_path)) {
    signals_tag <- sprintf("signals=%s", paste(signal_labels, collapse = "-"))
  
    tags <- c(chromosome_tag, "cells_by_half", "avg_score-vs-fraction", rho_tag, window_size_tag, nsamples_tag, signals_tag, weights_tag, domain_length_tag)
    fullname <- paste(c(dataset, tags), collapse = ",")
    filename <- sprintf("%s.png", fullname)
    fig_pathname <- file.path(fig_path, filename)

    ## Nothing to do?
    if (skip && file_test("-f", fig_pathname)) {
      res <- list()
      attr(res, "fig_pathname") <- normalizePath(fig_pathname)
      return(res)
    }
  }

  summary <- read_overlap_score_summary_vs_bin_size(dataset, chromosome = chromosome, bin_sizes = bin_sizes, rho = rho, window_size = window_size, nsamples = nsamples, weights = weights, domain_length = domain_length, ..., verbose = verbose)
  bin_size <- NULL; rm(list = "bin_size") ## To please R CMD check

  dw <- diff(range(bin_sizes)) / length(bin_sizes)

  signal_label <- names(known_signals)[signals[1] == known_signals]
  params <- c(sprintf("estimator: %s", signal_label),
              sprintf("weights: %s", weights),
              sprintf("domains: %.0f-%.0f", domain_length[1], domain_length[2]))
  subtitle <- sprintf("chromosome %s, test=%.3f, window size=%d (%d samples) [%s]",
                       chromosome, rho, window_size, nsamples, paste(params, collapse = "; "))

  ggs <- list()
  for (ss in seq_along(signals)) {
    signal <- signals[ss]
    if (signal %in% length_signals) {
      summary[[signal]] <- summary[[signal]] / 1000
    }
    gg <- ggplot(summary, aes_string(x = "bin_size", y = signal))
    gg <- gg + geom_boxplot(aes(group = as.factor(bin_size)), width = 0.2*dw)
    gg <- gg + geom_jitter(height = 0, width = 0.05*dw, size = 0.7, colour = "darkgray")
    gg <- gg + stat_summary(aes_string(y = signal, group = 1L),
                            fun.y = function(x) mean(x, trim = 0.10),
                            geom = "line",
                            colour = line_col, size = 2L, group = 1L)
    
    if (ss == 1L) {                              
      gg <- gg + ggtitle(label = NULL, subtitle = subtitle) + theme(plot.subtitle = element_text(size = 8))
      gg <- gg + xlab("")
      gg <- gg + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
      if ("chromosome" %in% labels) gg <- gg + geom_text_top_left(sprintf("Chr %s", chromosome), size=5)
      gg <- gg + geom_text_top_right(sprintf("%.0f%%", 100*rho), size=5)
    } else {
      gg <- gg + scale_x_continuous(labels = function(x) sprintf("%.0f", x/1000))
      gg <- gg + xlab("bin size (kbp)")
    }
    if (signal %in% length_signals) {
      ylim <- length_lim
      gg <- gg + ylab("domain length (kbp)")
      gg <- gg + theme(plot.margin = unit(c(5, 0, 5, 5), "pt"))
    } else {
      ylim <- c(0, 1)
      gg <- gg + ylab("average overlap score")
      gg <- gg + theme(plot.margin = unit(c(5, 5, -20, 5), "pt"))
    }
    gg <- gg + theme(axis.text = element_text(size = 13))
    gg <- gg + theme(axis.title = element_text(size = 13))

    ## xlim and ylim must be set at the same time
    gg <- gg + coord_cartesian(xlim = bin_size_lim, ylim = ylim)

    ggs[[signal]] <- gg    
  }

  res <- plot_grid(ggs[[1]], ggs[[2]], ncol = 1L, rel_heights = rel_heights, align = "v")

  if (!is.null(fig_pathname)) {
    if (verbose) suppressMessages <- identity
    suppressMessages(ggsave(res, filename = fig_pathname, scale = getOption("TopDomStudy.ggsave.scale", 0.80), dpi = getOption("TopDomStudy.ggsave.dpi", 360)))
    attr(res, "fig_pathname") <- normalizePath(fig_pathname)
    if (verbose) message(" - Figure written: ", fig_pathname)
  }

  if (verbose) message("gg_overlap_score_summary_vs_bin_size() ... done")

  res
} ## gg_overlap_score_summary_vs_bin_size()



#' @importFrom ggplot2 geom_text aes
geom_text_top_left <- function(label, size=5.0, hjust=-0.1, vjust=+1.5) {
  ## To please R CMD check
  xpos <- ypos <- NULL
  
  annotations <- data.frame(
    label = label,
    xpos = -Inf, ypos = +Inf,
    hjust = hjust, vjust = vjust
  )
  geom_text(data=annotations, aes(x=xpos, y=ypos, hjust=hjust, vjust=vjust, label=label), size=size)
}

#' @importFrom ggplot2 geom_text aes
geom_text_top_right <- function(label, size=5.0, hjust=+1.1, vjust=+1.5) {
  ## To please R CMD check
  xpos <- ypos <- NULL
  
  annotations <- data.frame(
    label = label,
    xpos = +Inf, ypos = +Inf,
    hjust = hjust, vjust = vjust
  )
  geom_text(data=annotations, aes(x=xpos, y=ypos, hjust=hjust, vjust=vjust, label=label), size=size)
}



#' @param tad_length_lim The range of the TAD lengths to be displayed.
#'
#' @rdname gg_overlap_score_summary_vs_fraction
#'
#' @importFrom grid unit
#' @importFrom ggplot2 aes aes_string coord_cartesian element_blank
#'             geom_boxplot geom_jitter geom_text ggplot ggsave ggtitle
#'             element_text scale_x_continuous stat_summary theme xlab ylab ylim
#' @importFrom scales percent_format
#' @importFrom utils file_test
#'
#' @export
gg_overlap_score_summary_vs_tad_length <- function(dataset, chromosome, bin_sizes, rhos, window_size = 5L, nsamples = 50L, weights = c("by_length", "uniform"), domain_length = NULL, signals = c("mean", "test_len_q0.50"), tad_length_lim = c(0, 2e6), labels = c("chromosome"), line_col = "black", fig_path = "figures", skip = FALSE, ..., verbose = TRUE) {
  avg_signals <- c(mean = "mean", median = "`50%`", "q25" = "`25%`", "q75" = "`75%`")
  length_signals <- c(
    "reference Q25 length"    = "ref_len_q0.25",
    "reference median length" = "ref_len_q0.50",
    "reference Q75 length"    = "ref_len_q0.75",
    "test Q25 length"         = "test_len_q0.25",
    "test median length"      = "test_len_q0.50",
    "test Q75 length"         = "test_len_q0.75"
  )
  known_signals <- c(avg_signals, length_signals)
  stopifnot(length(signals) == 2L, all(signals %in% known_signals))
  signal_labels <- character(2L)
  signal_labels[1] <- names(known_signals)[signals[1] == known_signals]
  signal_labels[2] <- signals[2]
  stop_if_not(length(signal_labels) == 2L)
  bin_sizes <- sort(bin_sizes)
  stop_if_not(all(bin_sizes > 0))
  rhos <- sort(rhos)
  stop_if_not(all(rhos > 0, rhos <= 1/2))
  
  if (verbose) {
    message("gg_overlap_score_summary_vs_bin_size() ...")
    message("- chromosome: ", chromosome)
    message("- window_size: ", window_size)
    message("- weights: ", weights)
    message("- nsamples: ", nsamples)
    message("- signals: ", paste(sQuote(signals), collapse = ", "))
    message("- signal labels: ", paste(sQuote(signal_labels), collapse = ", "))
  }

  if (!is.null(fig_path) && !file_test("-d", fig_path)) {
    dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)
    stop_if_not(file_test("-d", fig_path))
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
  r <- range(bin_sizes)
  bin_sizes_tag <- sprintf("bin_sizes=%d-%d", r[1], r[2])
  r <- range(rhos)
  rhos_tag <- sprintf("fractions=%.2f-%.2f", r[1], r[2])

  fig_pathname <- NULL
  if (!is.null(fig_path)) {
    signals_tag <- sprintf("signals=%s", paste(signal_labels, collapse = "-"))
  
    tags <- c(chromosome_tag, "cells_by_half", "avg_score-vs-fraction", window_size_tag, nsamples_tag, signals_tag, weights_tag, domain_length_tag, bin_sizes_tag, rhos_tag)
    fullname <- paste(c(dataset, tags), collapse = ",")
    filename <- sprintf("%s.png", fullname)
    fig_pathname <- file.path(fig_path, filename)

    ## Nothing to do?
    if (skip && file_test("-f", fig_pathname)) {
      res <- list()
      attr(res, "fig_pathname") <- normalizePath(fig_pathname)
      return(res)
    }
  }

  summary <- read_overlap_score_summary_vs_tad_length(dataset, chromosome = chromosome, bin_sizes = bin_sizes, rhos = rhos, window_size = window_size, nsamples = nsamples, weights = weights, domain_length = domain_length, ..., verbose = verbose)

  x_signal <- "test_len_q0.50"
#  summary[[x_signal]] <- summary[[x_signal]] / 1000
  x <- summary[["test_len_q0.50"]]
  dw <- diff(range(x)) / length(x)
  signal <- "mean"

  signal_label <- names(known_signals)[signals[1] == known_signals]
  params <- c(sprintf("estimator: %s", signal_label),
              sprintf("weights: %s", weights),
              sprintf("domains: %.0f-%.0f", domain_length[1], domain_length[2]))
  subtitle <- sprintf("chromosome %s, %s, %s, window size=%d, (%d samples) [%s]",
                       chromosome, bin_sizes_tag, rhos_tag, window_size, nsamples, paste(params, collapse = "; "))

  gg <- ggplot(summary, aes_string(x = x_signal, y = signal))
  gg <- gg + geom_boxplot(aes(group = as.factor(x)), width = 0.2*dw)
  gg <- gg + geom_jitter(height = 0, width = 0.05*dw, size = 0.7, colour = "darkgray")
  gg <- gg + stat_summary(aes_string(y = signal, group = 1L),
                          fun.y = function(x) mean(x, trim = 0.10),
                          geom = "line",
                          colour = line_col, size = 2L, group = 1L)
    
  gg <- gg + ggtitle(label = NULL, subtitle = subtitle) + theme(plot.subtitle = element_text(size = 8))
#  gg <- gg + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  if ("chromosome" %in% labels) gg <- gg + geom_text_top_left(sprintf("Chr %s", chromosome), size=5)
  gg <- gg + scale_x_continuous(labels = function(x) sprintf("%.0f", x/1000))
  gg <- gg + xlab("TAD length (kbp)")
  ylim <- c(0, 1)
  gg <- gg + ylab("average overlap score")
#  gg <- gg + theme(plot.margin = unit(c(5, 5, -20, 5), "pt"))
  gg <- gg + theme(axis.text = element_text(size = 13))
  gg <- gg + theme(axis.title = element_text(size = 13))

  ## xlim and ylim must be set at the same time
  gg <- gg + coord_cartesian(xlim = tad_length_lim, ylim = ylim)

  if (!is.null(fig_pathname)) {
    if (verbose) suppressMessages <- identity
    suppressMessages(ggsave(gg, filename = fig_pathname, scale = getOption("TopDomStudy.ggsave.scale", 0.80), dpi = getOption("TopDomStudy.ggsave.dpi", 360)))
    attr(gg, "fig_pathname") <- normalizePath(fig_pathname)
    if (verbose) message(" - Figure written: ", fig_pathname)
  }

  if (verbose) message("gg_overlap_score_summary_vs_bin_size() ... done")

  gg
} ## gg_overlap_score_summary_vs_tad_length()
