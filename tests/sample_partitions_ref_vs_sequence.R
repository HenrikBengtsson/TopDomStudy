library(TopDomStudy)

n <- 1000L
parts <- sample_partitions(n, fraction = 0.5)
str(parts)
idxs <- sort(unlist(parts))
stopifnot(identical(idxs, seq_len(n)))

n <- 1000L
seq <- seq(from = 0.1, to = 0.5, by = 0.1)
parts <- sample_partitions_ref_vs_seq(n, fraction = 0.5, seq = seq)
str(parts)
stopifnot(names(parts)[1] == "reference")
for (pp in 2:length(parts)) {
  stopifnot(length(intersect(parts[[pp]], parts$reference)) == 0L)
}

set.seed(0x42)
n <- 1000L
w <- runif(n)
parts <- sample_partitions_similar_weights(w, fraction = 1/3)
str(parts)
idxs <- sort(unlist(parts))
stopifnot(identical(idxs, seq_len(n)))


set.seed(0x42)
n <- 1000L
w <- runif(n)
seq <- seq(from = 0.1, to = 0.5, by = 0.1)
parts <- sample_partitions_similar_weights_ref_vs_seq(w, fraction = 0.5, seq = seq, w_tolerance = 0.01)
str(parts)
idxs <- sort(unlist(parts))
#stopifnot(identical(idxs, seq_len(n)))
