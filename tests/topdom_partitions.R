source("incl/start.R")

pathname <- system.file("compiledData", sprintf("%s,unique,chr=22.rds", dataset), package = "TopDomStudy", mustWork = TRUE)
reads <- read_rds(pathname)
print(reads)

for (partition_by in c("reads", "cells", "reads_by_half", "cells_by_half")) {
  res <- topdom_partitions(reads = reads, dataset = dataset, bin_size = 100000, partition_by = partition_by, min_cell_size = 2L, rho = 0.20, reference_rho = 1/2, nsamples = 2L, chrs = "22", seed = TRUE, mainseed = 0xBEEF, force = TRUE, verbose = TRUE)
  print(res)

  ## Calling it a second time should skip already existing results (on file)
  res2 <- topdom_partitions(reads = reads, dataset = dataset, bin_size = 100000, partition_by = partition_by, min_cell_size = 2L, rho = 0.20, reference_rho = "50%", nsamples = 2L, chrs = "22", seed = TRUE, mainseed = 0xBEEF, force = FALSE, verbose = TRUE)
  print(res2)
}

source("incl/end.R")
