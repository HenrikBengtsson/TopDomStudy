library(TopDomStudy)

pathname <- system.file("compiledData", "human,HAP1,unique,chr=22.rds", package = "TopDomStudy", mustWork = TRUE)
data <- read_rds(pathname)
print(data)

binned <- hic_bin(data, bin_size = 50000, intra_only = TRUE)
str(binned)
stopifnot(
  length(binned) == 1L,
  names(binned) == "22",
  is.matrix(binned[[1]]),
  is.integer(binned[[1]]),
  is.list(attr(binned[[1]], "bins"))
)

