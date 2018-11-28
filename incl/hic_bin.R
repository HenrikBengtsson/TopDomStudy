pathname <- system.file("compiledData", "human,HAP1,unique,chr=22.rds", package = "TopDomStudy")
reads <- read_rds(pathname)
print(reads)

binned <- hic_bin(reads, bin_size = 100000, intra_only = TRUE)
str(binned)

image(log2(binned[[1]]), axes = FALSE)
box()

