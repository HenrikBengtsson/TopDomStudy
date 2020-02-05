dataset <- "human,HAP1,unique"
pathname <- system.file("compiledData", sprintf("%s,chr=22.rds", dataset), package = "TopDomStudy")
reads <- read_rds(pathname)
print(reads)

td <- topdom_partitions(reads, bin_size = 100e3, rho = 0.5, nsamples = 5L,
                        partition_by = "cells_by_half",
                        dataset = dataset, verbose = TRUE)

