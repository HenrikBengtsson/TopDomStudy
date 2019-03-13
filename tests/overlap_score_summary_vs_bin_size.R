library(TopDomStudy)
library(future.apply)
plan(multiprocess, workers = max(1, 1/2 * availableCores()))

dataset <- "human,HAP1"
chromosome <- "22"
nsamples <- 10L
rho <- 0.25

filename <- sprintf("%s,unique,chr=%s.rds", dataset, chromosome)
pathname <- system.file("compiledData", filename, package = "TopDomStudy", mustWork = TRUE)
reads <- read_rds(pathname)
print(reads)

## Overlap scores per partition
bin_sizes <- c(10e3, 50e3, 100e3)

for (weights in c("uniform", "by_length")) {
  summary <- future_lapply(bin_sizes, FUN = function(bin_size) {
    res <- overlap_scores_partitions(reads = reads, dataset = "human,HAP1,unique", bin_size = bin_size, partition_by = "cells_by_half", min_cell_size = 2L, rho = rho, nsamples = nsamples, chrs = chromosome, seed = 0xBEEF, mainseed = 0xBEEF, force = TRUE)
    ## Overlap-score summaries
    summary_kk <- lapply(res[[chromosome]], FUN = function(pathname) {
      oss <- read_rds(pathname)
      ## Drop failed TopDom fits and possibly skip this sample?
      failed <- unlist(lapply(oss, FUN = inherits, "try-error"))
      if (any(failed)) {
        oss <- oss[!failed]
        if (length(oss) < 2) return(NULL)
      }
      z <- overlap_score_summary(oss, weights = weights)
      oss <- failed <- NULL
      
      pathname_td <- gsub("[.]rds$", ",topdom.rds", pathname)
      td <- read_rds(pathname_td)

      ref <- which(names(td) == "reference")
      sizes <- td[[ref]]$domain$size
      probs <- c(0.00, 0.05, 0.25, 0.50, 0.75, 0.95, 1.00)
      qsizes <- quantile(sizes, probs = probs, na.rm = TRUE)
      names(qsizes) <- sprintf("ref_len_q%0.2f", probs)
      z <- cbind(z, as.list(qsizes))

      sizes <- td[-ref][[1]]$domain$size
      probs <- c(0.00, 0.05, 0.25, 0.50, 0.75, 0.95, 1.00)
      qsizes <- quantile(sizes, probs = probs, na.rm = TRUE)
      names(qsizes) <- sprintf("test_len_q%0.2f", probs)
      z <- cbind(z, as.list(qsizes))

      z
    })
    summary_kk <- do.call(rbind, summary_kk)
    rownames(summary_kk) <- NULL
    summary_kk <- cbind(summary_kk, fraction = rho)
    summary_kk
  })
  summary <- do.call(rbind, summary)
  print(summary)

  if (require("ggplot2")) {
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

    for (signal_label in names(signals)) {
      signal <- signals[[signal_label]]
   
      gg <- ggplot(summary, aes_string(x = "bin_size", y = signal))
    
      gg <- gg + geom_boxplot(aes(group = as.factor(bin_size)), width = 0.2*dw)
      gg <- gg + geom_jitter(height = 0, width = 0.05*dw, size = 0.7, colour = "darkgray")
    
      gg <- gg + stat_summary(aes_string(y = signal, group = 1L),
                              fun.y = function(x) mean(x, trim = 0.10),
                              geom = "line", size = 2L, group = 1L)
  
      gg <- gg + ggtitle(dataset,
              subtitle = sprintf("chromosome %s, fraction=%.3f (%d samples) [estimator: %s; weights: %s]",
                                 chromosome, rho, nsamples, signal_label, weights))
      gg <- gg + xlab("bin size (bps)")
      if (signal_label %in% names(length_signals)) {
        gg <- gg + ylab("domain length (bps)")
        gg <- gg + ylim(0, 2e6)
      } else {
        gg <- gg + ylab("average overlap score")
        gg <- gg + ylim(0, 1)
      }

      print(gg)

      signal <- gsub("`50%`", "median", signal)
      filename <- sprintf("%s,chr=%s,%s,avg_score-vs-bin_size,fraction=%.3f,nsamples=%d,signal=%s,weights=%s.png", dataset, chromosome, "cells_by_half", rho, nsamples, signal_label, weights)

      ggsave(gg, filename = filename)
    } ## for (signal ...)
  }

} ## for (weights ...)
