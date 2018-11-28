library(TopDomStudy)

pathname <- system.file("compiledData", "human,HAP1,unique,chr=22.rds", package = "TopDomStudy", mustWork = TRUE)
reads <- read_rds(pathname)
print(reads)

for (partition_by in c("reads", "cells", "reads_by_half", "cells_by_half")) {
  res <- overlap_scores_partitions(reads = reads, dataset = "human,HAP1,unique", bin_size = 100000, partition_by = partition_by, min_cell_size = 1L, rho = 1/4, nsamples = 2L, chrs = "22", seed = TRUE, mainseed = 0xBEEF, verbose = TRUE)
  print(res)
}
