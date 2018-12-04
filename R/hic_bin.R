#' Bins HiC Reads in a Data Frame into Chromosome-Pair Count Matrices
#' 
#' @param data An N-by-K data.frame with
#' characters columns `chr_a`, `chr_b`, and
#' numeric columns `start_a`, `end_a`, `start_b`, and `end_b`.
#' 
#' @param bin_size (integer) The bin size (in basepairs).
#'
#' @param intra_only (logical) If `TRUE`, only intra-chromosome pairs are
#' binned.
#' 
#' @param known_chrs (character vector) The names of the L chromosomes
#' to bin over.
#' 
#' @param progress (logical) If `TRUE`, progress is displayed.
#'
#' @return
#' If `intra_only = FALSE`, an L-by-L list matrix where each element \eqn{(i, j)} in turn
#' holds an integer matrix of binned counts for chromosome pair \eqn{(i, j)}.
#' If `intra_only = TRUE`, an L list where each element \eqn{i} in turn
#' holds an integer matrix of binned counts for intra-chromosome pair \eqn{(i, i)}.
#'
#' @example incl/hic_bin.R
#'
#' @author Henrik Bengtsson adopted from a script by Adam Olshen.
#'
#' @importFrom listenv listenv
#' @importFrom future %<-% %label%
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
hic_bin <- function(data, bin_size, intra_only = FALSE, known_chrs = c(1:22, "X", "Y", "M"), progress = TRUE) {
  stop_if_not(
    all(c("chr_a", "start_a", "end_a", "chr_b", "start_b",  "end_b")
        %in% colnames(data))
  )
  
  bin_size <- as.integer(bin_size)
  stop_if_not(bin_size >= 1L)

  ## Focus only on intra-chromosome pairs?
  if (intra_only) {
    data <- subset(data, chr_a == chr_b)
  }
  
  known_chrs <- intersect(known_chrs, data$chr_a)

  if (progress) {
    max_progress <- if (intra_only) length(known_chrs) else length(known_chrs)^2
    max_progress <- max_progress + length(known_chrs)
    pb <- txtProgressBar(max = max_progress, file = stderr())
    iter <- 0L
  }

  ## Identify bins for each chromosome
  chr_bins <- lapply(known_chrs, FUN = function(chr) {
    chr_len <- max(data$end_a[data$chr_a == chr], data$end_b[data$chr_b == chr])
    seq(from = 0L, to = chr_len + bin_size, by = bin_size)
  })
  names(chr_bins) <- known_chrs

  ## Bin chr_a and chr_b positions independently into bin indices
  data$chr_a_bin <- data$chr_b_bin <- integer(nrow(data))
  for (chr in known_chrs) {
    if (progress) setTxtProgressBar(pb, value = (iter <- iter + 1))
    bins <- chr_bins[[chr]]

    ## chr_a
    rows <- which(data$chr_a == chr)
    if (length(rows) > 0) {
      pos <- with(data[rows, ], (start_a + end_a) / 2)
      data$chr_a_bin[rows] <- findInterval(pos, vec = bins)
    }

    ## chr_b
    rows <- which(data$chr_b == chr)
    if (length(rows) > 0) {
      pos <- with(data[rows, ], (start_b + end_b) / 2)
      data$chr_b_bin[rows] <- findInterval(pos, vec = bins)
    }

    ## Sanity check
    nbins <- length(bins) - 1L
    rows <- which(data$chr_a == chr)
    if (length(rows) > 0) {
      idxs <- data[rows, ]$chr_a_bin
      stop_if_not(all(idxs >= 1), all(idxs <= nbins))
    }

    rows <- which(data$chr_b == chr)
    if (length(rows) > 0) {
      idxs <- data[rows, ]$chr_b_bin
      stop_if_not(all(idxs >= 1), all(idxs <= nbins))
    }
  }

  ## Sanity checks
  stop_if_not(all(data$chr_a_bin >= 1), all(data$chr_b_bin >= 1))


  ## Count reads per (chr_a_bin, chr_b_bin) pair
  counts <- listenv()
  dim(counts) <- c(length(known_chrs), ncol = length(known_chrs))
  dimnames(counts) <- list(known_chrs, known_chrs)

  for (chr_a in known_chrs) {
    data_a <- data[which(data$chr_a == chr_a), ]
    bins <- list(chr_a = data_a$chr_a_bin, chr_b = NULL)
    chr_bins_a <- chr_bins[[chr_a]]
    nbins_a <- length(chr_bins_a) - 1L

    ## Sanity checks
    stop_if_not(all(data_a$chr_a == chr_a))
    stop_if_not(all(data_a$chr_a_bin >= 1), all(data_a$chr_a_bin <= nbins_a))

    chrs_b <- if (intra_only) chr_a else known_chrs
    for (chr_b in chrs_b) {
      if (progress) setTxtProgressBar(pb, value = (iter <- iter + 1))
      
      if (intra_only) {
        data_ab <- data_a  ## already done at the very top
      } else {
        data_ab <- data_a[which(data_a$chr_b == chr_b), ]
      }

      if (nrow(data_ab) == 0) {
        counts_t <- matrix(0L,
                           nrow = length(chr_bins_a) - 1L,
                           ncol = length(chr_bins_b) - 1L)
        attr(counts_t, "bins") <- list(chr_a = chr_bins_a, chr_b = chr_bins_b)
        counts[[chr_a, chr_b]] <- counts_t

        counts_t <- NULL ## Not needed anymore      
	next
      }	

      bins$chr_b <- data_ab$chr_b_bin
      chr_bins_b <- chr_bins[[chr_b]]
      nbins_b <- length(chr_bins_b) - 1L

      ## Sanity checks
      stop_if_not(all(data_ab$chr_a == chr_a), all(data_ab$chr_b == chr_b))
      stop_if_not(all(data_ab$chr_a_bin >= 1), all(data_ab$chr_a_bin <= nbins_a))
      stop_if_not(all(data_ab$chr_b_bin >= 1), all(data_ab$chr_b_bin <= nbins_b))

      counts[[chr_a, chr_b]] %<-% {
        counts_pair <- table(bins)
        idxs <- lapply(dimnames(counts_pair), FUN = as.integer)

        counts_t <- matrix(0L, nrow = nbins_a, ncol = nbins_b)
        counts_t[idxs$chr_a, idxs$chr_b] <- counts_pair
        attr(counts_t, "bins") <- list(chr_a = chr_bins_a, chr_b = chr_bins_b)
        counts_pair <- NULL ## Not needed anymore      
	counts_t
      }	%label% sprintf("%s-%s", chr_a, chr_b)
    }
  }

  ## Resolve
  counts <- as.list(counts)
  
  if (progress) {
    setTxtProgressBar(pb, value = nrow(counts))
    close(pb)
  }

  if (intra_only) counts <- diag(counts)
  
  counts
}
