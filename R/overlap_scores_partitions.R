#' Calculates TopDom Overlap Scores Across Partitions
#'
#' @param reads ...
#'
#' @param dataset ...
#'
#' @param cell_ids ...
#'
#' @param bin_size ...
#'
#' @param partition_by ...
#'
#' @param min_cell_size ...
#'
#' @param rho ...
#'
#' @param seed_tags ...
#'
#' @param nsamples = 100L ...
#'
#' @param chrs ...
#'
#' @param path_out ...
#'
#' @param save_topdom ...
#'
#' @param seed ...
#'
#' @param verbose ...
#'
#' @return A named list of length `length(chrs)`.
#'
#' @importFrom future value %<-% %label% %seed%
#' @importFrom future.apply future_lapply
#' @importFrom listenv listenv
#' @importFrom utils file_test str
#' @importFrom TopDom overlapScores TopDom
#' @export
overlap_scores_partitions <- function(reads, dataset, cell_ids, bin_size, partition_by, min_cell_size, rho, seed_tags, nsamples = 100L, chrs, path_out = ".", save_topdom = TRUE, seed = TRUE, verbose = FALSE) {
  ## To please R CMD check
  cell_id <- chr_a <- NULL; rm(list = c("cell_id", "chr_a"))
  
  stop_if_not(is.data.frame(reads) || inherits(reads, "Future"))
  stop_if_not(is.numeric(bin_size), length(bin_size) == 1L, !is.na(bin_size), bin_size > 0, is.finite(bin_size))
  stop_if_not(is.numeric(nsamples), length(nsamples) == 1L,
              !is.na(nsamples), nsamples >= 1L)
  stopifnot(length(min_cell_size) == 1L, is.numeric(min_cell_size),
            !is.na(min_cell_size), min_cell_size >= 1L)
  stop_if_not(is.numeric(rho), length(rho) == 1L, !is.na(rho), rho > 0.0, rho <= 0.5)
  stop_if_not(is.character(seed_tags), !anyNA(seed_tags))
  stop_if_not(file_test("-d", path_out))

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
    sample_idxs <- which(!file_test("-f", pathnames))
    
    ## Already done?
    if (length(sample_idxs) == 0L) {
      if (verbose) mprintf("Chromosome #%d (%s) of %d ... ALREADY DONE", cc, chr_tag, length(chrs))
      next
    }
    res[[chr]][sample_idxs] <- NA_character_
  
    if (verbose) mprintf(" - Number of (remaining) samples to process for chromosome %s (%s): %d", chr, chr_tag, length(sample_idxs))

    if (inherits(reads, "Future")) reads <- value(reads)

    res[[chr]] %<-% {      
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
        print(reads)
      }
  
      ## Subset by cell ids?
      if (!is.null(cell_ids)) {
        reads <- subset(reads, cell_id %in% cell_ids)
        print(reads)
      }
  
      ## Subset by chromosome
      reads <- subset(reads, chr_a %in% chr)
      print(reads)
  
      res_kk <- listenv()
  
      ## For each samples ...
      for (kk in seq_along(sample_idxs)) {
        bb <- sample_idxs[kk]
        pathname <- pathnames[bb]
        if (verbose) mprintf("Sample #%d (%s) of %d ...", kk, seed_tags[bb], length(sample_idxs))
    
        ## Already done? (should not happen, but just in case)
        if (file_test("-f", pathname)) {
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
  	}
          reads_partitions <- lapply(reads_partitions, FUN = function(partition) reads[partition, ])
          str(reads_partitions)
  
          ## TopDom on each partition
          tds <- lapply(reads_partitions, FUN = function(reads_pp) {
            counts <- {
              counts <- hic_bin(reads_pp, intra_only = TRUE, bin_size = bin_size, progress = FALSE)
              stopifnot(is.list(counts), length(counts) == 1L, all(names(counts) == chr))
              counts <- as_TopDomData(counts)
              counts
            }
            str(counts, nchar.max = 60L)
            stopifnot(is.list(counts), length(counts) == length(chr), all(names(counts) == chr))
  
            ## Used to be Try(TopDom), cf. https://github.com/HenrikBengtsson/TopDom/issues/4 
            tds <- future_lapply(counts, FUN = TopDom, window.size = 5L)
            stopifnot(is.list(tds), length(tds) == 1L, all(names(tds) == chr))
            counts <- NULL ## Not needed anymore
            
            tds[[chr]]
          })
  
          read_partitions <- NULL ## Not needed anymore
  
          ## Find first TopDom fit that didn't produce an error
          ok <- unlist(lapply(tds, FUN = function(td) !inherits(td, "try-error")))
          stopifnot(length(ok) == length(tds))
  
          ref <- which(ok)[1]
          td_ref <- tds[[ref]]
  
          if (save_topdom) {
  	  tt <- tds
            attr(tt, "bin_size") <- bin_size
            attr(tt, "chromosome") <- chr
            attr(tt, "min_cell_size") <- min_cell_size
            attr(tt, "reference_partition") <- ref
            attr(tt, "seed") <- seed
  	  attr(tt, "reference") <- ref
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
          attr(overlaps, "reference_partition") <- ref
          attr(overlaps, "seed") <- seed
          
          saveRDS(overlaps, file = pathname)
          print(overlaps)
  
          pathname
        } %seed% seed %label% paste(c(chr_tag, sprintf("sample=%d", kk)), collapse = "-")
  
        if (verbose) mprintf("Sample #%d (%s) of %d ... DONE", kk, seed_tags[bb], length(sample_idxs))
      } ## for (kk ...)
  
      reads <- NULL  ## Not needed anymore
  
      ## Resolve all samples for current chromosome
      res__kk <- unlist(res_kk)
  
      pathnames
    } %label% chr_tag
  
    if (verbose) mprintf("Chromosome #%d (%s) of %d ... DONE", cc, chr_tag, length(chrs))
  } ## for (chr ...)
  
  ## Resolve everything
  res <- as.list(res)
  
  if (verbose) print(res)
  
  res
} ## overlap_scores_partitions()
