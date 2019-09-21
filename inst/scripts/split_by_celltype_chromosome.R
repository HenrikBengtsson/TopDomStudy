library("utils")
library("progressr")     ## https://github.com/HenrikBengtsson/progressr
library("progress")

options(width = 140)
progressr::handlers("progress")

## NOTES: 6 human + 2 mouse cell types takes ~40 minutes to complete

path <- "compiledData"

celltypes <- list(
  human = c("HAP1", "HeLa", "GM12878", "K562", "Asynchronous", "Nocadazole"),
  mouse = c("MEF", "Patski")
)

chromosomes <- 1:23

with_progress({
  p <- progressr::progressor(sum(lengths(celltypes)) * (1 + length(chromosomes)))

for (org in names(celltypes)) {
  message(sprintf("Organism %s ...", sQuote(org)))

  pattern <- sprintf("%s,unique.rds", org)
  pathnames_org <- dir(path = path, pattern = pattern, full.names = TRUE)
  message(sprintf("Data files: [n = %d] %s", length(pathnames_org), paste(sQuote(pathnames_org), collapse = ", ")))
  
  for (celltype in celltypes[[org]]) {
    message(sprintf("Cell type %s ...", sQuote(celltype)))
    p(sprintf("Reading (%s,%s)", sQuote(org), sQuote(celltype)))

    filename <- sprintf("%s,%s,unique.rds", org, celltype)
    pathname_celltype <- file.path(path, filename)

    ## Source not available? Then skip with a warning.
    if (!file_test("-f", pathname_celltype)) {
      warning("Input file not found: ", sQuote(pathname_celltype), immediate. = TRUE)
      next
    }
    
    data <- readRDS(pathname_celltype)
    
    for (chr in chromosomes) {
      message(sprintf("Chromosome %s ...", sQuote(chr)))
      
      filename <- sprintf("%s,%s,unique,chr=%s.rds", org, celltype, chr)
      pathname_chr <- file.path(path, filename)
      p(sprintf("Saving (%s,%s,%s) intra chromosome", sQuote(org), sQuote(celltype), sQuote(chr)))
      
      ## Already processed?
      if (file_test("-f", pathname_chr)) next

      rows <- which(data$chr_a == chr & data$chr_b == chr)
      
      ## Nothing to save?
      if (length(rows) == 0L) next
      
      data_chr <- data[rows, , drop = FALSE]
      data <- data[-rows, , drop = FALSE]
      saveRDS(data_chr, file = pathname_chr)
      data_chr <- NULL
      
      message(sprintf("Chromosome %s ... OK", sQuote(chr)))
    } # for (chr in ...)
    
    message(sprintf("Cell type %s ... OK", sQuote(celltype)))
  } ## for (celltype in ...)
  
  message(sprintf("Organism %s ... OK", sQuote(org)))
} ## for (org in ...)

}) ## with_progress({ ... })
