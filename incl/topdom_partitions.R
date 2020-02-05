dataset <- "human,HAP1,unique,chr=22"
pathname <- system.file("compiledData", sprintf("%s.rds", dataset), package = "TopDomStudy")
reads <- read_rds(pathname)
print(reads)

td <- topdom_partitions(reads, bin_size = 100e3, rho = 0.5, nsamples = 5L,
                        partition_by = "cells_by_half",
                        dataset = dataset, verbose = TRUE)

