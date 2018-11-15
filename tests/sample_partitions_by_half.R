library(TopDomStudy)

n <- 1000L
parts <- sample_partitions(n, fraction = 0.5)
str(parts)
idxs <- sort(unlist(parts))
stopifnot(identical(idxs, seq_len(n)))

n <- 1000L
parts <- sample_partitions_by_half(n, fraction = 0.2)
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
parts <- sample_partitions_similar_weights_by_half(w, fraction = 0.4, w_tolerance = 0.01)
str(parts)
idxs <- sort(unlist(parts))
#stopifnot(identical(idxs, seq_len(n)))
