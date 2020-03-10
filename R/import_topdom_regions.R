#' Import all test-set TopDom regions into a single data.frame
#'
#' This script will read all available \file{topdomData/*/*.rds}
#' files, extract the TopDom regions for the test samples (but not the
#' reference) and convert to a data.frame with additional columns on parameter
#' settings and RNG seeds appended.  All these data.frames are stacked into one
#' big data.frame which is saved to a \file{overlapScoreData/*,test-topdom.rds}.
#'
#' @param pattern (character string) A regular expression of the subfolders
#' in folder `path` to import.
#'
#' @param path (character string) The root folder of the overlap score data.
#'
#' @param skip (logical) If TRUE, already existing results are returned, otherwise not.
#'
#' @param save_individual (logical) If TRUE, also the imported data.frames of
#' the subfolders are saved in individual \file{*,test-topdom.rds} file in the
#' `path` folder.
#'
#' @return The pathname of the saved \file{*,test-topdom.rds} file in the
#' `path` folder.
#'
#' @importFrom utils file_test
#' @importFrom future.apply future_lapply
#' @export
import_topdom_regions <- function(pattern = "human,HAP1,unique,bin_size=.*,partition_by=cells_by_half,min_cell_size=2,window_size=.*,test=.*,mainseed=0xBEEF", path = "topdomData", skip = TRUE, save_individual = TRUE) {
  stopifnot(file_test("-d", path))
  stopifnot(is.character(pattern), length(pattern) == 1L,
            !is.na(pattern), nzchar(pattern))
	    
  parts <- strsplit(pattern, split = ",", fixed = TRUE)[[1]]
  parts <- grep("*", parts, fixed = TRUE, invert = TRUE, value = TRUE)
  
  filename <- sprintf("%s,test-topdom.rds", paste(parts, collapse = ","))
  pathname <- file.path(path, filename)
  if (file_test("-f", pathname)) {
    if (skip) return(pathname)
    stop("File already exist: ", sQuote(pathname))
  }

  sets <- dir(path, pattern = pattern, no.. = TRUE)
  sets <- sets[file_test("-d", file.path(path, sets))]
  message("Number of sets: ", length(sets))
  stopifnot(length(sets) > 0L)

  data <- list()
  for (kk in seq_along(sets)) {
    set <- sets[kk]
    pathname_kk <- file.path(path, sprintf("%s,test-topdom.rds", set))
    message(sprintf("Set #%d (%s) of %d ...", kk, set, length(sets)))
    if (file_test("-f", pathname_kk)) {
      data_kk <- read_rds(pathname_kk)
      data[[kk]] <- data_kk
      message(sprintf("Set #%d (%s) of %d ... alread done", kk, set, length(sets)))
      next
    }
    pathnames <- dir(file.path(path, set), pattern = "[.]rds$", recursive = TRUE, full.names = TRUE)
    data_kk <- future_lapply(pathnames, FUN = read_topdom_regions)
    data_kk <- do.call(rbind, data_kk)
    print(data_kk)
    if (save_individual) save_rds(data_kk, pathname_kk)
    data[[kk]] <- data_kk
    message(sprintf("Set #%d (%s) of %d ... saved", kk, set, length(sets)))
  }

  data <- do.call(rbind, data)
  o <- with(data, order(chr, bin_size, fraction, window_size, seed, from.id))
  data <- data[o, ]
  save_rds(data, pathname)
  
  pathname
}
