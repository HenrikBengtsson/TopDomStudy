library(TopDomStudy)

dataset <- "human,HAP1,unique"
pathname <- system.file("compiledData", sprintf("%s,chr=22.rds", dataset), package = "TopDomStudy", mustWork = TRUE)
reads <- read_rds(pathname)
print(reads)

## Overlap scores per partition
dt <- system.time(res <- overlap_scores_partitions(reads = reads, dataset = dataset, bin_size = 100000, partition_by = "cells_by_half", min_cell_size = 2L, rho = 1/4, nsamples = 10L, chrs = "22", seed = 0xBEEF, mainseed = 0xBEEF, force = TRUE))
print(dt)

## TopDom data
pathname_oss <- res[["22"]][1]
oss <- read_rds(pathname_oss)

## Identify TopDom fit
set <- basename(dirname(pathname_oss))
path_td <- file.path("topdomData", set)
stop_if_not(file_test("-d", path_td))
filename_td <- basename(pathname_oss)
## Ad hoc /HB 2020-02-20
filename_td <- gsub(",reference_type=[^,]+,", ",", filename_td)
pathname_td <- file.path(path_td, filename_td)
stop_if_not(file_test("-f", pathname_td))
tds <- read_rds(pathname_td)
print(tds)

stopifnot(length(tds) == length(oss))

## The reference partition
ref <- attr(oss, "reference_partition", exact = TRUE)
td_ref <- tds[[ref]]

for (kk in seq_along(tds)) {
  td <- tds[[kk]]
  print(td$domain)
  
  os <- oss[[kk]]
  print(os)
  
  scores <- os[["22"]]$best_score
  print(scores)
  
  stopifnot(length(scores) == nrow(td_ref$domain))
}
