library(TopDomStudy)

pathname <- system.file("compiledData", "human,HAP1,unique,chr=22.rds", package = "TopDomStudy", mustWork = TRUE)
reads <- read_rds(pathname)
print(reads)

## Overlap scores per partition

summary <- NULL

chromosome <- "22"
rhos <- c(0.1, 0.2, 0.3, 0.4, 0.5)
for (kk in seq_along(rhos)) {
  rho <- rhos[kk]
  res <- overlap_scores_partitions(reads = reads, dataset = "human,HAP1,unique", bin_size = 100000, partition_by = "cells_by_half", min_cell_size = 2L, rho = rho, nsamples = 5L, chrs = chromosome, seed = 0xBEEF, mainseed = 0xBEEF)
  
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
