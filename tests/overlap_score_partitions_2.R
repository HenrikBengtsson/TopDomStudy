library(TopDomStudy)

pathname <- system.file("compiledData", "human,HAP1,unique,chr=22.rds", package = "TopDomStudy", mustWork = TRUE)
reads <- read_rds(pathname)
print(reads)

dt <- system.time(res <- overlap_scores_partitions(reads = reads, dataset = "human,HAP1,unique", bin_size = 100000, partition_by = "cells_by_half", min_cell_size = 1L, rho = 1/4, nsamples = 100L, chrs = "22", seed = TRUE, mainseed = 0xBEEF))
print(dt)

pathname <- res[["22"]][1]
oss <- read_rds(pathname)
print(oss)

pathname_td <- paste0(tools::file_path_sans_ext(pathname), ",topdom.rds")
tds <- read_rds(pathname_td)
print(tds)

stopifnot(length(tds) == length(oss))

for (kk in seq_along(tds)) {
  td <- tds[[kk]]
  os <- oss[[kk]]
  
  print(os)
  scores <- os[["22"]]$best_scores
  print(scores)
  
  stopifnot(length(scores) == nrow(td$domain))
}

