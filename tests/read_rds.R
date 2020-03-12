source("incl/start.R")

message("read_rds() ...")

message("read_rds() - RDS file")

x <- base::letters
tmpfile <- tempfile()
saveRDS(x, file = tmpfile)
y <- read_rds(tmpfile)
stopifnot(identical(y, x))

message("read_rds() - non-existing file")

res <- tryCatch({
  read_rds("non-existing-file.rds")
}, error = identity)
print(res)
stopifnot(inherits(res, "error"))


message("read_rds() - invalid file")

res <- tryCatch({
  read_rds(system.file("DESCRIPTION", package = "TopDomStudy"))
}, error = identity)
print(res)
stopifnot(inherits(res, "error"))

message("read_rds() ... done")

source("incl/end.R")
