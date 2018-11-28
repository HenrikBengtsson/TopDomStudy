pathname <- system.file("compiledData", "human,HAP1,unique,chr=22.rds", package = "TopDomStudy")
data <- read_rds(pathname)
print(data)

binned <- hic_bin(data, bin_size = 100000, intra_only = TRUE)
str(binned)

image(log2(binned[[1]]), axes = FALSE)
box()

