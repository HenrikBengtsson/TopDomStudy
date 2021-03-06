source("incl/start.R")

message("hic_bin() ...")

message("- reading reads")
pathname <- system.file("compiledData", sprintf("%s,unique,chr=22.rds", dataset), package = "TopDomStudy", mustWork = TRUE)
reads <- read_rds(pathname)
print(reads)

message("- binning")
counts <- hic_bin(reads, bin_size = 50000, intra_only = TRUE)
str(counts)
stopifnot(
  length(counts) == 1L,
  names(counts) == "22",
  is.matrix(counts[[1]]),
  is.integer(counts[[1]]),
  is.list(attr(counts[[1]], "bins", exact = TRUE))
)


progressr::with_progress({
  counts1 <- hic_bin(reads, bin_size = 50000, intra_only = TRUE)
})
str(counts1)
stopifnot(identical(counts1, counts))


message("- Coercing to list of TopDomData objects")
counts <- as_TopDomData(counts)
print(counts)
stopifnot(
  length(counts) == 1L,
  names(counts) == "22",
  inherits(counts[[1]], "TopDomData")
)

message("- Fit TomDom for each TopDomData object")
fits <- lapply(counts, FUN = TopDom, window.size = 5L)
print(fits)
stopifnot(
  length(fits) == 1L,
  names(fits) == "22",
  inherits(fits[[1]], "TopDom")
)

message("hic_bin() ... done")

source("incl/end.R")
