library(TopDomStudy)
library(future.apply)
plan(multiprocess, workers = 3/4 * availableCores())

dataset <- "human,HAP1"
chromosome <- "22"
nsamples <- 30L
bin_size <- 100000

filename <- sprintf("%s,unique,chr=%s.rds", dataset, chromosome)
pathname <- system.file("compiledData", filename, package = "TopDomStudy", mustWork = TRUE)
reads <- read_rds(pathname)
print(reads)

## Overlap scores per partition
summary <- NULL
## FIXME: chromosome = "22", rho = 0.01, bin_size = 10000, nsamples = 1L gives an error
rhos <- c(0.01, 0.02, 0.04, 0.05, 0.06, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50)[-1]
summary <- future_lapply(rhos, FUN = function(rho) {
  res <- overlap_scores_partitions(reads = reads, dataset = "human,HAP1,unique", bin_size = bin_size, partition_by = "cells_by_half", min_cell_size = 2L, rho = rho, nsamples = nsamples, chrs = chromosome, seed = 0xBEEF, mainseed = 0xBEEF)
  
  ## Overlap-score summaries
  summary_kk <- lapply(res[[chromosome]], FUN = function(pathname) {
    oss <- read_rds(pathname)
    ## Drop failed TopDom fits and possibly skip this sample?
    failed <- unlist(lapply(oss, FUN = inherits, "try-error"))
    if (any(failed)) {
      oss <- oss[!failed]
      if (length(oss) < 2) return(NULL)
    }
    overlap_score_summary(oss)
  })
  summary_kk <- do.call(rbind, summary_kk)
  rownames(summary_kk) <- NULL
  summary_kk <- cbind(summary_kk, fraction = rho)
  summary_kk
})
summary <- do.call(rbind, summary)

print(summary)

if (require("ggplot2")) {
  dw <- diff(range(rhos)) / length(rhos)
  gg <- ggplot(summary, aes(x = fraction, y = mean))
  
  gg <- gg + geom_boxplot(aes(group = as.factor(fraction)), width = 0.2*dw)
  gg <- gg + geom_jitter(height = 0, width = 0.05*dw, size = 0.7, colour = "darkgray")
  
  gg <- gg + stat_summary(aes(y = mean, group = 1L),
                          fun.y = function(x) mean(x, trim = 0.10),
                          geom = "line", size = 2L, group = 1L)
			  
  gg <- gg + ggtitle(dataset,
          subtitle = sprintf("chromosome %s, bin size=%d (%d samples)",
                             chromosome, bin_size, nsamples))
  gg <- gg + xlab("fraction") + ylab("average overlap score")
  gg <- gg + ylim(0,1)
  print(gg)

  ggsave(gg, filename=sprintf("%s,chr=%s,%s,avg_score-vs-fraction,bin_size=%d,nsamples=%d.png", dataset, chromosome, "cells_by_half", bin_size, nsamples))
}
