library("TopDomStudy") ## https://github.com/HenrikBengtsson/TopDomStudy
library("ramani")      ## https://github.com/HenrikBengtsson/ramani
library("progressr")   ## https://github.com/HenrikBengtsson/progressr
library("progress")
library("dplyr")
library("future.apply")
library("listenv")

options(width = 140)
progressr::handlers("progress")


## NOTES: GSM2438426_ML4,human takes ~14 GiB RAM and ~25 mins to finish

## ---------------------------------------------------------------------
## Setup data
## ---------------------------------------------------------------------
names <- c(
  "GSM2254215_ML1", "GSM2254216_ML2", "GSM2254217_ML3",
  "GSM2254218_PL1", "GSM2254219_PL2", "GSM2438426_ML4"
)
organisms <- c("human", "mouse")

res <- listenv()
dim(res) <- c(length(names), length(organisms))
dimnames(res) <- list(names, organisms)

with_progress({
  p <- progressr::progressor(length(names) + (2+7)*length(res))
  
for (name in names) {
  message(sprintf("Sample %s ...", sQuote(name)))
  p(sprintf("Reading 'percentages' for %s", sQuote(name)))

  file <- sprintf("hicData/GSE84920/%s.percentages.txt.gz", name)
  per <- read_percentages(file) ## Takes 1-2 sec to read (270 kB)
  
  ## ---------------------------------------------------------------------
  ## Subset of barcodes to focus on
  ## ---------------------------------------------------------------------
  ## The 'celltype' field in "percentages" is "Undetermined" unless either
  ## hg19_frac or mm10_frac >= 0.95.  See ?ramani::read_percentages.
  per <- filter(per, celltype != "Undetermined")
  print(per)
  
  print(nrow(per))
  ## [1] 13987
  
  print(table(per$celltype, useNA = "always"))
  ##   HAP1   HeLa    MEF Patski   <NA>
  ##   3413   3686   3842   3046      0
  
  
  
  ## ---------------------------------------------------------------------
  ## Split by organism and annotate with celltypes
  ## ---------------------------------------------------------------------
  path_dest <- "compiledData"
  dir.create(path_dest, recursive = TRUE, showWarnings = FALSE)

  for (org in c("human", "mouse")) {
    message("Organism: ", org)
  
    filename <- sprintf("%s,%s,unique.rds", name, org)
    pathname <- file.path(path_dest, filename)

    ## Already done?
    if (file_test("-f", pathname)) {
      res[[name, org]] <- pathname
      p(sprintf("Already done (%s,%s)", name, org), amount=9L)
      next
    }

    p(sprintf("Processing (%s,%s)", name, org))

    ## Keep only celltypes for organism of interest
    per_org <- switch(org,
      human = filter(per, hg19_frac >= 0.95),
      mouse = filter(per, mm10_frac >= 0.95)
    )

    print(c(table(per_org$celltype), total = nrow(per_org)))
    ## human:
    ##  HAP1  HeLa total
    ##  3413  3686  7099
    ##
    ## mouse:
    ##   MEF Patski  total
    ##  3842   3046   6888
      
    ## Assert that (inner_barcode, outer_barcode) -> celltype is unique
    t <- select(per_org, c(inner_barcode, outer_barcode, celltype))
    stopifnot(all(unique(t)[,1:2] == unique(t[,1:2])))
    
    celltypes <- unique(select(per_org, c(inner_barcode, outer_barcode, celltype)))
    print(celltypes)
    ## # A tibble: 3,272 x 3
    ##    inner_barcode outer_barcode celltype
    ##    <chr>         <chr>         <chr>
    ##  1 ACCACCAC      TCAGATGC      HAP1
    ##  2 CATAGCGC      ACTTGATA      HAP1
    ##  3 GGCCGTTC      GCCATTAA      HAP1
    ##  4 GTCGGCAT      TTGACCAT      HeLa
    ##  5 GTTCCACC      CGTTACTT      HAP1
    ##  6 GAGCTCGA      CACTAGCT      HeLa
    ##  7 GACTGCCA      TCACTGAG      HAP1
    ##  8 CTCGCCGA      GATTCTTA      HAP1
    ##  9 CGCGAGTA      CTGATTCT      HeLa
    ## 10 CGATGCTC      ATATCAGA      HAP1
    ## # ... with 3,262 more rows

    p(sprintf("Processing (%s,%s) for %d cells", name, org, nrow(celltypes)))

    res[[name, org]] %<-% {
      ## Valid pairs for this organism
      path <- file.path("hicData", "GSE84920")
      filename <- sprintf("%s.validPairs.txt.gz", name)
      pathname <- file.path(path, filename)
      p(sprintf("Reading 'validPairs' for %s: %s", name, filename))
      vp <- read_validpairs(pathname, columns = c("chr_a", "start_a", "end_a", "chr_b", "start_b", "end_b", "inner_barcode", "outer_barcode"))
      
      ## Subset for organism
      p(sprintf("Subsetting %s for organism %s", name, org))
      vp <- vp[grepl(org, vp$chr_a, fixed = TRUE) & grepl(org, vp$chr_b, fixed = TRUE), ]
      vp <- as_tibble(vp)
      print(vp)
      ## # A tibble: 18,823,611 x 8
      ##    chr_a         start_a     end_a chr_b         start_b     end_b inner_barcode outer_barcode
      ##    <chr>           <int>     <int> <chr>           <int>     <int> <chr>         <chr>
      ##  1 human_chr2   88637036  88637116 human_chr2   89109729  89109866 GAGGAGCA      CGATGACA
      ##  2 human_chr4  126924974 126925025 human_chr4  127334997 127335124 GCTACGGT      AGTCGTAT
      ##  3 human_chr15  42130039  42130103 human_chr15  42677209  42677290 AGGTGCGA      ATACATGT
      ##  4 human_chr1  209393848 209393905 human_chr1  232468122 232468159 GCCTCGAA      GAGTACGT
      ##  5 human_chr5    5565796   5565876 human_chr5  149955845 149955889 GCTCGCTA      CTAGTGAA
      ##  6 human_chr1  181632778 181632810 human_chr1  181731968 181732207 CAGGCTTG      GATATAAC
      ##  7 human_chr13  73644136  73644240 human_chr13  74667914  74667979 TCGGATCG      GATATAAC
      ##  8 human_chr2   80491806  80491879 human_chr2   80575049  80575117 CTCGCCGA      TACCGGAA
      ##  9 human_chr1  168466976 168467108 human_chr1  168467668 168467791 GAGCTCGA      TCCGTCTA
      ## 10 human_chr18  70871038  70871080 human_chr18  71498300  71498395 TCCGTGAG      ATATCAGA
      ## # ... with 18,823,601 more rows
      
      
      p(sprintf("Joining subsetted 'validPairs' with 'celltypes'"))
      data <- left_join(vp, celltypes, by = c("inner_barcode", "outer_barcode"))
      pattern <- switch(org, human = "human_chr", mouse = "mouse_ch")
      data <- mutate(data, chr_a = gsub(pattern, "", chr_a, fixed = TRUE), chr_b = gsub(pattern, "", chr_b, fixed = TRUE))
      print(data)
      ## # A tibble: 18,823,611 x 9
      ##    chr_a   start_a     end_a chr_b   start_b     end_b inner_barcode outer_barcode celltype
      ##    <chr>     <int>     <int> <chr>     <int>     <int> <chr>         <chr>         <chr>   
      ##  1 2      88637036  88637116 2      89109729  89109866 GAGGAGCA      CGATGACA      HeLa    
      ##  2 4     126924974 126925025 4     127334997 127335124 GCTACGGT      AGTCGTAT      HAP1    
      ##  3 15     42130039  42130103 15     42677209  42677290 AGGTGCGA      ATACATGT      HeLa    
      ##  4 1     209393848 209393905 1     232468122 232468159 GCCTCGAA      GAGTACGT      HeLa    
      ##  5 5       5565796   5565876 5     149955845 149955889 GCTCGCTA      CTAGTGAA      HeLa    
      ##  6 1     181632778 181632810 1     181731968 181732207 CAGGCTTG      GATATAAC      HeLa    
      ##  7 13     73644136  73644240 13     74667914  74667979 TCGGATCG      GATATAAC      NA      
      ##  8 2      80491806  80491879 2      80575049  80575117 CTCGCCGA      TACCGGAA      HAP1    
      ##  9 1     168466976 168467108 1     168467668 168467791 GAGCTCGA      TCCGTCTA      NA      
      ## 10 18     70871038  70871080 18     71498300  71498395 TCCGTGAG      ATATCAGA      HeLa    
      ## # ... with 18,823,601 more rows
    
      ## Assert that no duplicates were introduced
      stopifnot(nrow(data) == nrow(vp))
      vp <- NULL  ## Not needed anymore
      
      ## Arrange by chromosome and position
      ## (mostly for neatness, but also because RDS file shrink to 40%!)
      p(sprintf("Arrange by chromosome and position"))
      data <- arrange(data, chr_a, start_a, chr_b, start_b)
    
      print(c(table(data$celltype, useNA = "always"), total = nrow(data)))
      ## human:
      ##     HAP1     HeLa     <NA>     total
      ##  6071095 11001283  1751233  18823601
      ##
      ## mouse:
      ##      MEF   Patski     <NA>     total
      ##   552055  7466099   745860   8764004
    
      ## Clean up: drop "funny" chromosomes
      good_chrs <- c(1:23, "X", "Y", "M")
      data <- filter(data, chr_a %in% good_chrs & chr_b %in% good_chrs)
    
      ## Clean up: map (inner_barcode, outer_barcode) <-> cell_index
      barcodes <- data[, c("inner_barcode", "outer_barcode")]
      barcodes <- with(barcodes, paste(inner_barcode, outer_barcode, sep = "-"))
      data$cell_id <- as.factor(barcodes)
      data$inner_barcode <- data$outer_barcode <- NULL
      barcodes <- NULL  ## Not needed anymore
      
      filename <- sprintf("%s,%s.rds", name, org)
      pathname <- file.path(path_dest, filename)
      p(sprintf("Saving to file: %s", filename))
      saveRDS(data, file = pathname)
    
      ## Keep only unique read pairs
      p(sprintf("Dropping duplicated read pairs"))
      ## NOTE: unique() on a data.frame can be very memory hungry.
      ## Because of this, we run unique() per chromosome.
      uchrs <- unique(data$chr_a)
      udata <- vector("list", length = length(uchrs))
      for (ii in seq_along(uchrs)) {
        chr <- uchrs[ii]
	rows <- which(data$chr_a == chr)
	data_ii <- data[rows, , drop = FALSE]
	data <- data[-rows, , drop = FALSE]
	udata[[ii]] <- unique(data_ii)
	data_ii <- NULL
      }
      data <- NULL
      data <- do.call(rbind, udata)
      udata <- NULL
      filename <- sprintf("%s,%s,unique.rds", name, org)
      pathname <- file.path(path_dest, filename)
      p(sprintf("Saving to file: %s", filename))
      saveRDS(data, file = pathname)

      data <- NULL  ## Not needed anymore
      
      pathname
    } %label% paste(name, org, sep = "-") ## %<-%    
  } ## for (org ...)
  
  message(sprintf("Sample %s ... OK", sQuote(name)))
} ## for (name ...)
print(res)
res <- as.list(res)
}) ## with_progress({ ... })

print(res)
