#' Calculates TopDom Overlap Scores Across Partitions
#'
#' @param reads A [base::data.frame].
#'
#' @param bin_size A positive numeric.
#'
#' @param partition_by A string specifying how to partition;
#'        one of `"reads"`, `"cells"`, `"reads_by_half"`, or `"cells_by_half"`.
#'
#' @param rho A numeric in \eqn{(0,1/2]} specifying the relative size of the partions.
#'
#' @param nsamples Number of random samples.
#'
#' @param seed Random seed for reproducible (parallel) random number generation (RNG).
#'
#' @param min_cell_size (optional, filter) The minimum number of reads for a cell to
#'        be included. Cells with less reads are dropped.
#'
#' @param chrs (optional, filter) Names of chromosomes to iterate over.
#'        Defaults to the chromsomes in `reads$chr_a`.
#'
#' @param cell_ids (optional, filter) ...
#'
#' @param dataset (optional) ...
#'
#' @param path_out The root folder that will contain the `overlapScoreData/` folder
#'        to which RDS files are written.
#'
#' @param save_topdom ...
#'
#' @param mainseed ...
#'
#' @param as Should values or pathnames be returned?
#'
#' @param force If `FALSE`, already processed partitions are skipped, otherwise not.
#'
#' @param verbose If `TRUE`, verbose message are produced, otherwise not.
#'
#' @return A named list of length `length(chrs)` with names as `chrs`.
#' Each list elements contains `nsamples` pathnames of RDS files.
#'
#' @section Parallel processing:
#' Internally, [future.apply::future_lapply] is used to parallelize
#' over [TopDom::TopDom].
#'
#' @importFrom future value %<-% %label% %seed%
#' @importFrom future.apply future_lapply
#' @importFrom listenv listenv
#' @importFrom utils file_test str
#' @importFrom TopDom overlapScores TopDom
#' @export
overlap_scores_partitions <- function(reads, bin_size, partition_by, rho, nsamples = 100L, seed = TRUE,
                                      chrs = NULL, min_cell_size = 1L, dataset, cell_ids = NULL,
				      path_out = ".", save_topdom = TRUE, mainseed = 0xBEEF, force = FALSE,
				      as = c("pathname", "value"),
				      verbose = FALSE) {
  ## To please R CMD check
  cell_id <- chr_a <- NULL; rm(list = c("cell_id", "chr_a"))
  
  stop_if_not(is.data.frame(reads) || is.function(reads))
  stop_if_not(is.numeric(bin_size), length(bin_size) == 1L, !is.na(bin_size), bin_size > 0, is.finite(bin_size))
  partition_by <- match.arg(partition_by, choices = c("reads", "cells", "reads_by_half", "cells_by_half"))
  stopifnot(length(min_cell_size) == 1L, is.numeric(min_cell_size),
            !is.na(min_cell_size), min_cell_size >= 1L)
  stop_if_not(is.numeric(rho), length(rho) == 1L, !is.na(rho), rho > 0.0, rho <= 0.5)
  stop_if_not(is.numeric(nsamples), length(nsamples) == 1L,
              !is.na(nsamples), nsamples >= 1L)
  if (is.null(chrs)) {
    chrs_a <- sort(unique(reads$chr_a))
    chrs_b <- sort(unique(reads$chr_b))  ## don't use this
    chrs <- chrs_a
  }
  stop_if_not(is.character(chrs), length(chrs) >= 1L, !anyNA(chrs))
  chrs <- sort(unique(chrs))

  ## Argument tags
  if (!is.null(cell_ids)) {
    stopifnot(partition_by == "reads", is.character(cell_ids), !anyNA(cell_ids))
    cell_ids_tag <- sprintf("cell_ids=%s", paste(cell_ids, collapse = "_"))
  } else {
    cell_ids_tag <- NULL
  }
  bin_size_tag <- sprintf("bin_size=%g", bin_size)
  partition_by_tag <- sprintf("partition_by=%s", partition_by)
  if (min_cell_size > 1L) {
    min_cell_size_tag <- sprintf("min_cell_size=%d", min_cell_size)
  } else {
    min_cell_size_tag <- NULL
  }
  rho_tag <- sprintf("fraction=%.3f", rho)
  ## Random seeds (use the same for all chromosomes, i.e. invariant to chromosome)
  ## FIXME: Export make_rng_seeds()
  if (is.list(seed)) {
    stop_if_not(length(seeds) == nsamples)
    seeds <- seed
  } else {
    seeds <- future.apply:::make_rng_seeds(nsamples, seed = seed)
  }
  seed_tags <- sprintf("seed=%s", sapply(seeds, FUN = crc32))

  stop_if_not(mainseed == 0xBEEF)
  mainseed_tag <- "mainseed=0xBEEF"

  as <- match.arg(as)

  dataset_out <- paste(c(dataset, cell_ids_tag, bin_size_tag, partition_by_tag, min_cell_size_tag, rho_tag, mainseed_tag), collapse = ",")
  path_out <- file.path("overlapScoreData", dataset_out)
  dir.create(path_out, recursive = TRUE, showWarnings = FALSE)
  stop_if_not(file_test("-d", path_out))

  if (verbose) message("Output path: ", path_out)
  
  res <- listenv()
  length(res) <- length(chrs)
  names(res) <- chrs
  
  ## For each chromosome ...
  for (cc in seq_along(chrs)) {
    chr <- chrs[cc]
    chr_tag <- sprintf("chr=%s", chr)
    if (verbose) mprintf("Chromosome #%d (%s) of %d ...", cc, chr_tag, length(chrs))

    ## Find all input files for this chromosome
    pathnames <- file.path(path_out, sapply(1:nsamples, function(bb) {
      tags <- c(cell_ids_tag, chr_tag, bin_size_tag, partition_by_tag, min_cell_size_tag, rho_tag, seed_tags[bb])
      name <- paste(c(dataset, tags), collapse = ",")
      sprintf("%s.rds", name)
    }))
  
    res[[chr]] <- pathnames
  
    ## Identify samples to be done
    if (force) {
      sample_idxs <- seq_along(pathnames)
    } else {
      done <- file_test("-f", pathnames)
      sample_idxs <- which(!done)
      if (as == "value") {
        res[[chr]][done] <- lapply(pathnames[done], FUN = read_rds)
      }
    }
    
    ## Already done?
    if (length(sample_idxs) == 0L) {
      if (verbose) mprintf("Chromosome #%d (%s) of %d ... ALREADY DONE", cc, chr_tag, length(chrs))
      next
    }
    res[[chr]][sample_idxs] <- NA_character_
  
    if (verbose) mprintf(" - Number of (remaining) samples to process for chromosome %s (%s): %d", chr, chr_tag, length(sample_idxs))

    res[[chr]] %<-% {
      if (is.function(reads)) reads <- reads()

      ## Subset by minimum (whole-genome) cell size?
      if (min_cell_size > 1L) {
        cell_sizes <- table(as.character(reads$cell_id))
        ncells <- length(cell_sizes)
        
        ## Drop non-existing cell_id:s due to empty levels
        cells_keep <- which(cell_sizes >= min_cell_size)
        cell_sizes <- cell_sizes[cells_keep]
        stopifnot(min(cell_sizes) >= min_cell_size)
        
        reads_keep <- which(reads$cell_id %in% names(cell_sizes))
        nreads <- nrow(reads)
        if (verbose) {
	  mprintf("Dropped %d (%.2f%%) cells (out of %d) with less than %d reads each resulting in dropping %d (%.2f%%) reads (out of %d)",
                  ncells - length(cells_keep), 100 * (ncells - length(cells_keep))/ncells, ncells,
                  min_cell_size,
                  nreads - length(reads_keep), 100 * (nreads - length(reads_keep))/nreads, nreads
          )
	}
  	      
        reads <- reads[reads_keep, ]
        stopifnot(length(setdiff(names(cell_sizes), reads$cell_id)) == 0L,
                  length(setdiff(reads$cell_id, names(cell_sizes))) == 0L)
        cell_sizes <- table(as.character(reads$cell_id))
        stopifnot(min(cell_sizes) >= min_cell_size)
        cell_sizes <- cells_keep <- reads_keep <- NULL ## Not needed anymore
      }
  
      ## Subset by cell ids?
      if (!is.null(cell_ids)) {
        reads <- subset(reads, cell_id %in% cell_ids)
      }
  
      ## Subset by chromosome
      reads <- subset(reads, chr_a %in% chr)
      if (verbose) print(reads)
  
      res_kk <- listenv()
  
      ## For each samples ...
      for (kk in seq_along(sample_idxs)) {
        bb <- sample_idxs[kk]
        pathname <- pathnames[bb]
        if (verbose) mprintf("Sample #%d (%s) of %d ...", kk, seed_tags[bb], length(sample_idxs))
    
        ## Already done? (should not happen, but just in case)
        if (!force && file_test("-f", pathname)) {
          if (verbose) mprintf("Sample #%d (%s) of %d ... ALREADY DONE", kk, seed_tags[bb], length(sample_idxs))
          next
        }
        if (verbose) message("- Output pathname: ", pathname)
  
        seed <- seeds[[bb]]
        if (verbose) {
	  message("- Random seed:")
          str(seed)
	}
        
        res_kk[[bb]] %<-% {
          if (verbose) mprintf(" - Random, disjoint partitioning of %s", partition_by)
          if (partition_by == "reads") {
            reads_partitions <- sample_partitions(nrow(reads), fraction = rho)
          } else if (partition_by == "cells") {
            reads_partitions <- sample_partitions_by_cells(reads, fraction = rho)
          } else if (partition_by == "reads_by_half") {
            reads_partitions <- sample_partitions_by_half(nrow(reads), fraction = rho)
          } else if (partition_by == "cells_by_half") {
            reads_partitions <- sample_partitions_by_cells_by_half(reads, fraction = rho)
          }
          reads_partitions <- lapply(reads_partitions, FUN = function(partition) reads[partition, ])
  
          ## TopDom on each partition
          tds <- lapply(reads_partitions, FUN = function(reads_pp) {
            counts <- {
              counts <- hic_bin(reads_pp, intra_only = TRUE, bin_size = bin_size, progress = FALSE)
              stopifnot(is.list(counts), length(counts) == 1L, all(names(counts) == chr))
              counts <- as_TopDomData(counts)
              counts
            }
            stopifnot(is.list(counts), length(counts) == length(chr), all(names(counts) == chr))
  
            ## Used to be Try(TopDom), cf. https://github.com/HenrikBengtsson/TopDom/issues/4 
            tds <- future_lapply(counts, FUN = TopDom, window.size = 5L)
            stopifnot(is.list(tds), length(tds) == 1L, all(names(tds) == chr))
            
            tds_chr <- tds[[chr]]
            if (save_topdom) attr(tds_chr, "counts") <- counts

            counts <- NULL ## Not needed anymore

            tds_chr
          })
  
          read_partitions <- NULL ## Not needed anymore

          if (partition_by %in% c("reads_by_half", "cells_by_half")) {
	    ref <- 1L
	  } else {
            ## Find first TopDom fit that didn't produce an error
            ok <- unlist(lapply(tds, FUN = function(td) !inherits(td, "try-error")))
            stopifnot(length(ok) == length(tds))
            ref <- which(ok)[1]
	  }
  
          td_ref <- tds[[ref]]
  
          if (save_topdom) {
            tt <- tds
            attr(tt, "bin_size") <- bin_size
            attr(tt, "chromosome") <- chr
            attr(tt, "min_cell_size") <- min_cell_size
            attr(tt, "reference_partition") <- ref
            attr(tt, "seed") <- seed
            attr(tt, "reference") <- ref
            attr(tt, "partition_by") <- partition_by
            pathname2 <- sprintf("%s,topdom.rds", tools::file_path_sans_ext(pathname))
            saveRDS(tt, file = pathname2)
            tt <- NULL
          }
 
          overlaps <- lapply(tds, FUN = function(td) Try(overlapScores)(td, td_ref))
          stopifnot(is.list(overlaps), length(overlaps) == length(tds))
          tds <- NULL ## Not needed anymore
          attr(overlaps, "bin_size") <- bin_size
          attr(overlaps, "chromosome") <- chr
          attr(overlaps, "min_cell_size") <- min_cell_size
          attr(overlaps, "partition_by") <- partition_by
          attr(overlaps, "reference_partition") <- ref
          attr(overlaps, "seed") <- seed
          
          saveRDS(overlaps, file = pathname)
          if (verbose) print(overlaps)
  
          pathname
        } %seed% seed %label% paste(c(chr_tag, sprintf("sample=%d", kk)), collapse = "-")
  
        if (verbose) mprintf("Sample #%d (%s) of %d ... DONE", kk, seed_tags[bb], length(sample_idxs))
      } ## for (kk ...)
  
      reads <- NULL  ## Not needed anymore
  
      ## Resolve all samples for current chromosome
      res__kk <- unlist(res_kk)

      if (as == "value") {
        lapply(pathnames, FUN = read_rds)
      } else {
        pathnames
      }	
    } %label% chr_tag
  
    if (verbose) mprintf("Chromosome #%d (%s) of %d ... DONE", cc, chr_tag, length(chrs))
  } ## for (chr ...)
  
  ## Resolve everything
  res <- as.list(res)
  
  if (verbose) print(res)

  res
} ## overlap_scores_partitions()
