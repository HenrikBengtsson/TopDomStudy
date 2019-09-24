#' Split Per-Organism Compiled File Data into Subsets of Cell Types
#'
#' @param celltypes A named list of character vectors.
#'
#' @param path The folder where input files are located and where output files
#' are written.
#'
#' @return A named list of the same structure as `celltypes` but
#' where the character string are pathnames to files produced.
#'
#' @details
#' This function operates on files produces by [compile_by_organism].
#'
#' @examples
#' \donttest{\dontrun{
#' progressr::with_progress({
#'   files <- TopDomStudy::split_by_celltype(celltypes=list(human="HAP1"),
#'                                           path="compiledData")
#' })
#' print(files)
#' # $human
#' #                                 HAP1 
#' # "compiledData/human,HAP1,unique.rds"
#' }}
#'
#' @section Progress updates:
#' This function signals [progressr::progression] updates. To visualize,
#' or in other ways render, progress information, wrap the call inside a
#' [progressr::with_progress] call.
#'
#' @importFrom utils file_test
#' @importFrom future.apply future_lapply
#' @importFrom dplyr arrange filter
#' @importFrom progressr progressor
#' @export
split_by_celltype <- function(celltypes = list(
                                human = c("HAP1", "HeLa", "GM12878", "K562",
				          "Asynchronous", "Nocadazole"),
                                mouse = c("MEF", "Patski")
			     ), path = "compiledData") {
  ## Dummy globals to please R CMD check
  celltype <- chr_a <- start_a <- chr_b <- start_b <- NULL
  
  stopifnot(is.list(celltypes), !is.null(names(celltypes)),
            names(celltypes) %in% c("human", "mouse"),
            all(vapply(celltypes, FUN = anyDuplicated, FUN.VALUE = 0L) == 0L))
  if (!file_test("-d", path)) dir.create(path, recursive = TRUE)
  stopifnot(file_test("-d", path))
  
  p <- progressor(4 * sum(lengths(celltypes)))

  res <- lapply(celltypes, FUN = function(x) {
    y <- rep(NA_character_, times = length(x))
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
  
      ## Already processed?
      if (file_test("-f", pathname_celltype)) {
        res[[org]][celltype] <- pathname_celltype
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

      ## Not needed anymore
      data <- NULL
      
      message(sprintf("Cell type %s ... OK", sQuote(celltype)))

      res[[org]][celltype] <- pathname_celltype
    } ## for (celltype in ...)
    
    message(sprintf("Organism %s ... OK", sQuote(org)))
  } ## for (org in ...)

  res
}

