library(TopDomStudy)

parts <- sample_partitions(100, fraction = 0.5)
str(parts)
idxs <- sort(unlist(parts))
stopifnot(identical(idxs, 1:100))

parts <- sample_partitions_ref_vs_seq(100, fraction = 0.5, seq = seq(from = 0.1, to = 0.5, by = 0.1))
str(parts)
stopifnot(names(parts)[1] == "reference")
for (pp in 2:length(parts)) {
  stopifnot(length(intersect(parts[[pp]], parts$reference)) == 0L)
}
