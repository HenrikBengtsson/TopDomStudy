library(TopDomStudy)

iris1 <- unique(datasets::iris)
iris2 <- unique_by(datasets::iris, by = "Species")
stopifnot(identical(iris2, iris1))
