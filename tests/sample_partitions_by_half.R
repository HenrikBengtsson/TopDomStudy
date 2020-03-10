library(TopDomStudy)

n <- 1000L

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Equivalence to sample_partitions(fraction = 0.5)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set.seed(0xBEEF)
parts0 <- sample_partitions(n, fraction = 0.5)
str(parts0)
idxs <- sort(unlist(parts0, use.names = FALSE))
stopifnot(
  length(parts0) == 2L,
  length(idxs) == n,
  identical(idxs, seq_len(n))
)

n <- 1000L
set.seed(0xBEEF)
parts <- sample_partitions_by_half(n, fraction = c(reference = 1/2, test = 1/2))
str(parts)
idxs <- sort(unlist(parts, use.names = FALSE))
stopifnot(
  grepl("^reference", names(parts)[1]),
  length(parts) == 2L,
  length(idxs) == n,
  identical(idxs, seq_len(n)),
  all(parts[[1]] == parts0[[1]])
)
for (pp in 2:length(parts)) {
  stopifnot(length(
    intersect(parts[[pp]], parts[[1]])) == 0L,
    length(parts[[pp]]) == length(parts[[1]]),
    length(parts[[pp]]) == length(parts0[[pp]]),
    all(sort(parts[[pp]]) == sort(parts0[[pp]]))  ## note sort()
  )
}


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Equivalence to sample_partitions_similar_weights(fraction = 0.5)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set.seed(0xBEEF)
n <- 1000L
w <- runif(n)

set.seed(0xBEEF)
parts0 <- sample_partitions_similar_weights(w, fraction = 1/2)
str(parts0)
idxs <- sort(unlist(parts0, use.names = FALSE))
stopifnot(
  length(parts0) == 2L,
  length(idxs) == n,
  identical(idxs, seq_len(n))
)

set.seed(0xBEEF)
parts <- sample_partitions_similar_weights_by_half(w, fraction = c(reference = 1/2, test = 1/2), w_tolerance = 0.01)
str(parts)
idxs <- sort(unlist(parts, use.names = FALSE))
stopifnot(
  grepl("^reference", names(parts)[1]),
  length(parts) == 2L,
  length(idxs) == n,
  identical(idxs, seq_len(n)),
  all(parts[[1]] == parts0[[1]])
)
for (pp in 2:length(parts)) {
  stopifnot(length(
    intersect(parts[[pp]], parts[[1]])) == 0L,
    length(parts[[pp]]) == length(parts[[1]]),
    length(parts[[pp]]) == length(parts0[[pp]]),
    all(sort(parts[[pp]]) == sort(parts0[[pp]]))  ## note sort()
  )
}


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Misc.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (fraction in c(0.01, 0.2, 0.5)) {
  set.seed(0xBEEF)
  parts <- sample_partitions_by_half(n, fraction = c(reference = 1/2, test = fraction))
  str(parts)
  idxs <- sort(unlist(parts, use.names = FALSE))
  stopifnot(
    grepl("^reference", names(parts)[1]),
    length(parts) == 2L,
    all(parts[[1]] == parts0[[1]]),
    length(idxs) == (0.5 + fraction) * n
  )
  for (pp in 2:length(parts)) {
    stopifnot(
      length(intersect(parts[[pp]], parts[[1]])) == 0L,
      length(parts[[pp]]) == fraction * n
    )
  }
}


for (fraction in c(0.01, 0.2, 0.5)) {
  set.seed(0xBEEF)
  parts <- sample_partitions_similar_weights_by_half(w, fraction = c(reference = 1/2, test = fraction), w_tolerance = 0.01)
  str(parts)
  idxs <- sort(unlist(parts, use.names = FALSE))
  stopifnot(
    grepl("^reference", names(parts)[1]),
    length(parts) == 2L,
    all(parts[[1]] == parts0[[1]]),
    length(idxs) == (0.5 + fraction) * n
  )
  for (pp in 2:length(parts)) {
    stopifnot(
      length(intersect(parts[[pp]], parts[[1]])) == 0L,
      length(parts[[pp]]) == fraction * n
    )
  }
}

