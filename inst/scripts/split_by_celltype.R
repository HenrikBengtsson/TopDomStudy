library("dplyr")
library("future.apply")
library("progressr")   ## https://github.com/HenrikBengtsson/progressr
library("progress")

options(width = 140)
progressr::handlers("progress")

## NOTES: 6 human + 2 mouse cell types takes ~40 minutes to complete

path <- "compiledData"

celltypes <- list(
  human = c("HAP1", "HeLa", "GM12878", "K562", "Asynchronous", "Nocadazole"),
  mouse = c("MEF", "Patski")
)

with_progress({
  p <- progressr::progressor(4 * sum(lengths(celltypes)))

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

    ## Already processed?
    if (file_test("-f", pathname_celltype)) {
      next
    }

    data <- future_lapply(pathnames_org, FUN = function(pathname) {
      name <- gsub(",.*", "", basename(pathname))
      message(sprintf("Sample: %s (%s)", name, pathname))
      data <- readRDS(pathname)
      type <- celltype
      data <- filter(data, celltype == type)
      
      ## No such celltype in current sample?
      if (nrow(data) == 0) return(NULL)
      
      data$name <- name
      data
    })

    p(sprintf("Merging (%s,%s)", sQuote(org), sQuote(celltype)))
    data <- do.call(rbind, data)
    
    ## Sort by (chr, position); makes the saved file smaller
    p(sprintf("Sorting (%s,%s)", sQuote(org), sQuote(celltype)))
    data <- arrange(data, chr_a, start_a, chr_b, start_b)
    
    p(sprintf("Saving (%s,%s)", sQuote(org), sQuote(celltype)))
    saveRDS(data, file = pathname_celltype)
    
    message(sprintf("Cell type %s ... OK", sQuote(celltype)))
  } ## for (celltype in ...)
  
  message(sprintf("Organism %s ... OK", sQuote(org)))
} ## for (org in ...)

}) ## with_progress({ ... })

