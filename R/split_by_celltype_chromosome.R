#' Split Per-Cell Type Compiled File Data into Per Chromosome
#'
#' @param celltypes A named list of character vectors.
#'
#' @param chromosomes Integer vector of chromosomes to be considered.
#'
#' @param path The folder where input files are located and where output files
#' are written.
#'
#' @return A named list of the same structure as `celltypes` but
#' where the character string are pathnames to files produced.
#'
#' @details
#' This function operates on files produces by [split_by_celltype].
#'
#' @examples
#' \donttest{\dontrun{
#' progressr::with_progress({
#'   files <- TopDomStudy::split_by_celltype_chromosome(
#'                 celltypes=list(human="HAP1"), chromosomes=c(12,16,22),
#'                 path="compiledData")
#' })
#' print(files)
#' # $human
#' # $human$HAP1
#' #                                      chr=12 
#' # "compiledData/human,HAP1,unique,chr=12.rds" 
#' #                                      chr=16 
#' # "compiledData/human,HAP1,unique,chr=16.rds" 
#' #                                      chr=22 
#' # "compiledData/human,HAP1,unique,chr=22.rds"
#' }}
#'
#' @section Progress updates:
#' This function signals [progressr::progression] updates. To visualize,
#' or in other ways render, progress information, wrap the call inside a
#' [progressr::with_progress] call.
#'
#' @importFrom utils file_test
#' @importFrom progressr progressor
#' @export
split_by_celltype_chromosome <- function(celltypes = list(
                                human = c("HAP1", "HeLa", "GM12878", "K562",
				          "Asynchronous", "Nocadazole"),
                                mouse = c("MEF", "Patski")
			     ), chromosomes = 1:23, path = "compiledData") {
  stopifnot(is.list(celltypes), !is.null(names(celltypes)),
            names(celltypes) %in% c("human", "mouse"),
            all(vapply(celltypes, FUN = anyDuplicated, FUN.VALUE = 0L) == 0L))
  stopifnot(is.numeric(chromosomes), length(chromosomes) >= 1L)
  
  if (!file_test("-d", path)) dir.create(path, recursive = TRUE)
  stopifnot(file_test("-d", path))

  p <- progressor(sum(lengths(celltypes)) * (1 + length(chromosomes)))

  res <- lapply(celltypes, FUN = function(x) {
    y <- lapply(x, FUN = function(x) character(0L))
    names(y) <- x
    y
  })

  organisms <- names(celltypes)
  for (org in organisms) {
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
        chr_tag <- sprintf("chr=%s", chr)
        message(sprintf("Chromosome %s (%s) ...", sQuote(chr), sQuote(chr_tag)))
        
        filename <- sprintf("%s,%s,unique,%s.rds", org, celltype, chr_tag)
        pathname_chr <- file.path(path, filename)
        p(sprintf("Saving (%s,%s,%s) intra chromosome", sQuote(org), sQuote(celltype), sQuote(chr_tag)))
        
        ## Already processed?
        if (file_test("-f", pathname_chr)) {
          res[[org]][[celltype]][chr_tag] <- pathname_chr
          next
	}
  
        rows <- which(data$chr_a == chr & data$chr_b == chr)
        
        ## Nothing to save?
        if (length(rows) == 0L) next
        
        data_chr <- data[rows, , drop = FALSE]
        data <- data[-rows, , drop = FALSE]
        saveRDS(data_chr, file = pathname_chr)

        ## Not needed anymore
        data_chr <- NULL

        res[[org]][[celltype]][chr_tag] <- pathname_chr
	
        message(sprintf("Chromosome %s (%s) ... OK", sQuote(chr), sQuote(chr_tag)))
      } # for (chr in ...)

      ## Not needed anymore
      data <- NULL

      message(sprintf("Cell type %s ... OK", sQuote(celltype)))
    } ## for (celltype in ...)
    
    message(sprintf("Organism %s ... OK", sQuote(org)))
  } ## for (org in ...)
  
  res
}
