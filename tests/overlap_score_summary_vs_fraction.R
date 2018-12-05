library(TopDomStudy)

dataset <- "human,HAP1"
chromosome <- "22"
nsamples <- 35L

filename <- sprintf("%s,unique,chr=%s.rds", dataset, chromosome)
pathname <- system.file("compiledData", filename, package = "TopDomStudy", mustWork = TRUE)
reads <- read_rds(pathname)
print(reads)

## Overlap scores per partition

summary <- NULL


rhos <- c(0.1, 0.2, 0.3, 0.4, 0.5)
for (kk in seq_along(rhos)) {
  rho <- rhos[kk]
  res <- overlap_scores_partitions(reads = reads, dataset = "human,HAP1,unique", bin_size = 100000, partition_by = "cells_by_half", min_cell_size = 2L, rho = rho, nsamples = nsamples, chrs = chromosome, seed = 0xBEEF, mainseed = 0xBEEF)
  
  ## Overlap-score summaries
  summary_kk <- lapply(res[[chromosome]], FUN = function(pathname) {
    oss <- read_rds(pathname)
    overlap_score_summary(oss)
  })
  summary_kk <- do.call(rbind, summary_kk)
  rownames(summary_kk) <- NULL
  summary_kk <- cbind(summary_kk, fraction = rho)
  summary <- rbind(summary, summary_kk)
}

print(summary)

if (require("ggplot2")) {
  gg <- ggplot(summary, aes(x = as.factor(fraction), y = mean))
  
  gg <- gg + geom_boxplot(width = 0.1)
  gg <- gg + geom_jitter(height = 0, width = 0.1, colour = "darkgray")

  gg <- gg + stat_summary(aes(y = mean, group = 1L),
                          fun.y = function(x) mean(x, trim = 0.10),
                          geom = "line", size = 2L, group = 1L)
			  
  gg <- gg + ggtitle(dataset,
          subtitle = sprintf("chromosome %s (%d partitions)",
                             chromosome, nsamples))
  gg <- gg + xlab("fraction") + ylab("average overlap score")
  gg <- gg + ylim(0,1)
  print(gg)
}
