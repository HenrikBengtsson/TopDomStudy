%\VignetteIndexEntry{TopDomStudy: Data Preprocessing}
%\VignetteAuthor{Henrik Bengtsson}
%\VignetteEngine{TopDomStudy::self}
%\VignetteEncoding{UTF-8}
%\VignetteKeyword{R}
%\VignetteKeyword{package}
%\VignetteKeyword{vignette}

In this work, we study the human HAP1 cell types part of the Ramani et al. (2017) data set. In addition to provide a set of R functions to help reproduce the results in our study, this package also provides a set of pre-processed data files to simply reproducibility in the `system.file("compiledData", package="TomDopStudy")` folder, which are used for this study, were generated from the Ramani et al. (2017) data set.
This document describes how these data files where produced.

## Overview

<!--
path <- system.file("compiledData", package="TomDopStudy")
path <- "inst/compiledData"
filenames <- dir(path = path, pattern = ",chr=.+[.]rds$")
chrs <- sub(".*,chr=(.+)[.]rds$", "\\1", filenames)
filenames <- filenames[gtools::mixedorder(chrs)]
chrs <- sub(".*,chr=(.+)[.]rds$", "\\1", filenames)
pathnames <- file.path(path, filenames)
sizes <- sapply(pathnames, FUN = file.size)
data <- lapply(pathnames, FUN = readRDS)
names(data) <- sprintf("Chr %s", chrs)
bps <- sapply(data, FUN = function(df) max(df$end_b) - min(df$start_a))
nreads <- sapply(data, FUN = nrow)
ncells <- sapply(data, FUN = function(df) length(unique(df$cell_id)))
files <- data.frame(chromosome = chrs, bps = bps, cells = ncells, reads = nreads, size = sizes, filename = filenames, stringsAsFactors = FALSE)
files <- files[gtools::mixedorder(files$chromosome),]
rownames(files) <- NULL

tbl <- files[,c("chromosome", "bps", "cells", "reads", "size")]
total <- list("Total", NA_real_, NA_integer_, sum(tbl$reads), sum(tbl$size))
tbl <- rbind(tbl, total)
tbl$chromosome <- sprintf("Chr %s", tbl$chromosome)
tbl$bps <- round(tbl$bps/1e6, digits=1L)
tbl$size <- round(tbl$size/1e6, digits=1L)
colnames(tbl) <- c("Chromosome", "Length (Mbps)", "Unique Cells", "Unique Read Pairs", "File size (MB)")
options(knitr.kable.NA = "")
knitr::kable(tbl, format.args = list(big.mark = ","))
-->

|Chromosome | Length (Mbps)| Unique Cells| Unique Read Pairs| File size (MB)|
|:----------|-------------:|------------:|-----------------:|--------------:|
|Chr 1      |         249.1|        1,872|           898,587|           10.9|
|Chr 2      |         243.2|        1,882|           972,893|           11.7|
|Chr 3      |         197.9|        1,854|           812,916|            9.8|
|Chr 4      |         191.0|        1,837|           744,741|            9.0|
|Chr 5      |         180.8|        1,833|           696,648|            8.4|
|Chr 6      |         170.8|        1,861|           662,910|            8.0|
|Chr 7      |         159.1|        1,847|           591,839|            7.2|
|Chr 8      |         146.1|        1,815|           569,888|            6.9|
|Chr 9      |         141.1|        1,797|           417,634|            5.1|
|Chr 10     |         135.4|        1,812|           508,691|            6.2|
|Chr 11     |         134.8|        1,818|           524,712|            6.4|
|Chr 12     |         133.8|        1,851|           537,659|            6.5|
|Chr 13     |          96.1|        1,781|           372,950|            4.5|
|Chr 14     |          88.3|        1,794|           345,941|            4.2|
|Chr 15     |          82.5|        1,804|           400,169|            4.6|
|Chr 16     |          90.2|        1,785|           265,621|            3.3|
|Chr 17     |          81.2|        1,783|           274,831|            3.4|
|Chr 18     |          78.0|        1,756|           290,220|            3.5|
|Chr 19     |          59.0|        1,756|           186,474|            2.3|
|Chr 20     |          62.9|        1,757|           232,448|            2.8|
|Chr 21     |          38.7|        1,745|           118,787|            1.5|
|Chr 22     |          34.9|        1,718|           112,704|            1.4|
|Chr 23     |         153.6|        1,827|           557,273|            6.8|
|Chr 24     |           0.1|            8|                 8|            0.0|
|Chr 25     |           0.0|           24|                83|            0.0|
|Total      |              |             |        11,096,627|          134.5|

_Table S1: Summary of HiC read pair data across chromosomes._


## Distribution of cell sizes

<!--
cell_sizes <- lapply(data, FUN = function(x) {
  t <- table(x$cell_id)
  t <- t[t > 0]
  t <- sort(t, decreasing=TRUE)
  t
})

distr <- lapply(cell_sizes, FUN = function(sizes) {
  summary(as.integer(sizes))
})
distr <- do.call(rbind, distr)
distr <- as.data.frame(distr)
distr <- lapply(distr, FUN = round, digits = 1L)

tbl <- files[,c("chromosome", "cells")]
tbl <- cbind(tbl, distr)
tbl$chromosome <- sprintf("Chr %s", tbl$chromosome)
colnames(tbl)[1:2] <- c("Chromosome", "Unique Cells")
options(knitr.kable.NA = "")
knitr::kable(tbl, format.args = list(big.mark = ","))
-->

|Chromosome | Unique Cells| Min.| 1st Qu.| Median|  Mean| 3rd Qu.|   Max.|
|:----------|------------:|----:|-------:|------:|-----:|-------:|------:|
|Chr 1      |        1,872|    1|    98.0|  252.0| 480.0|   580.0| 15,326|
|Chr 2      |        1,882|    1|   105.0|  279.0| 516.9|   629.0| 16,057|
|Chr 3      |        1,854|    1|    96.2|  235.0| 438.5|   539.5| 13,581|
|Chr 4      |        1,837|    1|    90.0|  225.0| 405.4|   485.0| 11,867|
|Chr 5      |        1,833|    1|    84.0|  210.0| 380.1|   463.0| 11,455|
|Chr 6      |        1,861|    1|    76.0|  194.0| 356.2|   434.0| 10,870|
|Chr 7      |        1,847|    1|    69.0|  173.0| 320.4|   394.0|  9,905|
|Chr 8      |        1,815|    1|    69.0|  173.0| 314.0|   377.0|  8,993|
|Chr 9      |        1,797|    1|    53.0|  124.0| 232.4|   278.0|  7,068|
|Chr 10     |        1,812|    1|    61.8|  150.0| 280.7|   340.0|  8,494|
|Chr 11     |        1,818|    1|    64.0|  154.0| 288.6|   353.0|  8,915|
|Chr 12     |        1,851|    1|    62.0|  156.0| 290.5|   350.0|  9,048|
|Chr 13     |        1,781|    1|    51.0|  120.0| 209.4|   254.0|  5,850|
|Chr 14     |        1,794|    1|    44.0|  105.0| 192.8|   234.8|  6,081|
|Chr 15     |        1,804|    1|    47.0|  116.5| 221.8|   273.2|  6,552|
|Chr 16     |        1,785|    1|    32.0|   79.0| 148.8|   181.0|  4,757|
|Chr 17     |        1,783|    1|    32.0|   78.0| 154.1|   187.0|  5,054|
|Chr 18     |        1,756|    1|    39.0|   94.0| 165.3|   198.2|  4,679|
|Chr 19     |        1,756|    1|    21.0|   51.0| 106.2|   125.0|  3,641|
|Chr 20     |        1,757|    1|    30.0|   71.0| 132.3|   161.0|  4,105|
|Chr 21     |        1,745|    1|    17.0|   38.0|  68.1|    82.0|  1,982|
|Chr 22     |        1,718|    1|    13.0|   32.0|  65.6|    78.0|  2,108|
|Chr 23     |        1,827|    1|    67.5|  166.0| 305.0|   369.0|  8,851|
|Chr 24     |            8|    1|     1.0|    1.0|   1.0|     1.0|      1|
|Chr 25     |           24|    1|     1.0|    2.0|   3.5|     3.5|     17|

_Table S2: Summary of number of read pairs per unique cell across chromosomes._


![](/home/hb/sf/TopDomStudy/figures/human,HAP1,Chr1,cell_size,histogram.png)

_Figure S1: Histogram of human HAP1 cell sizes (number of read-pairs per cell) on Chr 1._


<!--
single_read_cell_stats <- lapply(cell_sizes, FUN = function(n) { is_single <- (n == 1L); c(n_cells=length(n), single_cells=sum(is_single), prop_single_cells=sum(is_single)/length(n), n_reads=sum(n), prop_single_reads=sum(n[is_single])/sum(n)) })
tbl <- do.call(rbind, single_read_cell_stats)
tbl <- cbind(data.frame(chromosome = rownames(tbl), stringsAsFactors = FALSE), tbl)
rownames(tbl) <- NULL
colnames(tbl)[1:6] <- c("Chromosome", "Unique Cells", "Single-Read Cells", "Fraction Single-Read Cells", "Read Pairs", "Fraction Single Reads")
cols <- c("Fraction Single-Read Cells", "Fraction Single Reads")
tbl[cols] <- lapply(tbl[cols], FUN = signif, digits = 2L)
options(knitr.kable.NA = "")
knitr::kable(tbl, format.args = list(big.mark = ","))
-->

|Chromosome | Unique Cells| Single-Read Cells| Fraction Single-Read Cells| Read Pairs| Fraction Single Reads|
|:----------|------------:|-----------------:|--------------------------:|----------:|---------------------:|
|Chr 1      |        1,872|               126|                      0.067|    898,587|               0.00014|
|Chr 2      |        1,882|               135|                      0.072|    972,893|               0.00014|
|Chr 3      |        1,854|               126|                      0.068|    812,916|               0.00015|
|Chr 4      |        1,837|                97|                      0.053|    744,741|               0.00013|
|Chr 5      |        1,833|               102|                      0.056|    696,648|               0.00015|
|Chr 6      |        1,861|               115|                      0.062|    662,910|               0.00017|
|Chr 7      |        1,847|               117|                      0.063|    591,839|               0.00020|
|Chr 8      |        1,815|                90|                      0.050|    569,888|               0.00016|
|Chr 9      |        1,797|                75|                      0.042|    417,634|               0.00018|
|Chr 10     |        1,812|                87|                      0.048|    508,691|               0.00017|
|Chr 11     |        1,818|                95|                      0.052|    524,712|               0.00018|
|Chr 12     |        1,851|               115|                      0.062|    537,659|               0.00021|
|Chr 13     |        1,781|                59|                      0.033|    372,950|               0.00016|
|Chr 14     |        1,794|                73|                      0.041|    345,941|               0.00021|
|Chr 15     |        1,804|                78|                      0.043|    400,169|               0.00019|
|Chr 16     |        1,785|                62|                      0.035|    265,621|               0.00023|
|Chr 17     |        1,783|                66|                      0.037|    274,831|               0.00024|
|Chr 18     |        1,756|                44|                      0.025|    290,220|               0.00015|
|Chr 19     |        1,756|                60|                      0.034|    186,474|               0.00032|
|Chr 20     |        1,757|                51|                      0.029|    232,448|               0.00022|
|Chr 21     |        1,745|                46|                      0.026|    118,787|               0.00039|
|Chr 22     |        1,718|                46|                      0.027|    112,704|               0.00041|
|Chr 23     |        1,827|                96|                      0.053|    557,273|               0.00017|
|Chr 24     |            8|                 8|                      1.000|          8|               1.00000|
|Chr 25     |           24|                11|                      0.460|         83|               0.13000|

_Table S3: Summary of cells with a single read pair._





The Ramani data set is published on NCBI's Gene Expression Omnibus (GEO) in the GEO series [GSE84920] \(titled 'Massively multiplex single-cell Hi-C'\), which contains:

| GEO Sample      | GEO Title                        | Cell Types                                             |
| --------------- | -------------------------------- | ------------------------------------------------------ |
| [GSM2254215]    | Combinatorial scHi-C Library ML1 | human ('HAP1', 'HeLa'), mouse ('MEF', 'Patski')        |
| [GSM2254216]    | Combinatorial scHi-C Library ML2 | human ('HAP1', 'HeLa'), mouse ('MEF', 'Patski')        |
| [GSM2254217]    | Combinatorial scHi-C Library ML3 | human ('GM12878', 'K562'), mouse ('MEF', 'Patski')     |
| [GSM2254218]    | Combinatorial scHi-C Library PL1 | human ('HAP1', 'HeLa'), mouse ('MEF', 'Patski')        |
| [GSM2254219]    | Combinatorial scHi-C Library PL2 | human ('HAP1', 'HeLa'), mouse ('MEF', 'Patski')        |
| [GSM2438426](*) | Combinatorial scHi-C Library ML4 | human ('Asynchronous', 'Nocadazole'), mouse ('Patski') |

_Table S3: Overview of the content in the six GEO samples part of GEO series GSE84920._

## Step 1. Downloading Raw Data

In this study, we focus on the human HAP1 cell types which data is available in four out of the above six data sets.  For each of the four data sets, there are three files we need to download.  We download all of the 12 (=4*3) files using the `system.file("scripts", "download.sh", package="TopDomStudy")` script:

```sh
#!/usr/bin/env bash

url_path="https://www.ncbi.nlm.nih.gov/geo/download"
samples=(GSM2254215_ML1 GSM2254216_ML2 GSM2254218_PL1 GSM2254219_PL2)
types=(percentages validPairs assignments)
dest_path=hicData/GSE84920
mkdir -p "$dest_path"
for sample in "${samples[@]}"; do
  for type in "${types[@]}"; do
    file=$sample.$type.txt.gz
    echo "File: $file"
    if [[ ! -f "$dest_path/$file" ]]; then
      url="$url_path/?acc=${sample//_*}&format=file&file=$file"
      curl "$url" -o "$dest_path/$file"
    fi
  done
done
```

Running this Bash script, e.g.

```sh
$ path=$(Rscript -e "cat(system.file('scripts', package='TopDomStudy'))")
$ $path/download.sh
...
```

will download the twelve `*.txt.gz` files to local folder `hicData/GSE84920/`:

```sh
$ ls -l hicData/GSE84920/
total 5533616
-rw-r--r-- 1 alice alice        362 Oct  8 13:49 GSM2254215_ML1.assignments.txt.gz
-rw-r--r-- 1 alice alice     273357 Oct  8 13:46 GSM2254215_ML1.percentages.txt.gz
-rw-r--r-- 1 alice alice 1224864620 Oct  8 13:49 GSM2254215_ML1.validPairs.txt.gz
-rw-r--r-- 1 alice alice        362 Oct  8 13:57 GSM2254216_ML2.assignments.txt.gz
-rw-r--r-- 1 alice alice     199468 Oct  8 13:54 GSM2254216_ML2.percentages.txt.gz
-rw-r--r-- 1 alice alice 1192510493 Oct  8 13:57 GSM2254216_ML2.validPairs.txt.gz
-rw-r--r-- 1 alice alice        362 Oct  8 13:59 GSM2254218_PL1.assignments.txt.gz
-rw-r--r-- 1 alice alice     466679 Oct  8 13:57 GSM2254218_PL1.percentages.txt.gz
-rw-r--r-- 1 alice alice 1278669926 Oct  8 13:59 GSM2254218_PL1.validPairs.txt.gz
-rw-r--r-- 1 alice alice        362 Oct  8 14:02 GSM2254219_PL2.assignments.txt.gz
-rw-r--r-- 1 alice alice     524826 Oct  8 13:59 GSM2254219_PL2.percentages.txt.gz
-rw-r--r-- 1 alice alice 1968865275 Oct  8 14:02 GSM2254219_PL2.validPairs.txt.gz
```


## Step 2. Extract Human Data

In R, call:

```r
progressr::with_progress({
  files <- TopDomStudy::compile_by_organism(
             samples=c("GSM2254215_ML1", "GSM2254219_PL2",
                       "GSM2254216_ML2", "GSM2254218_PL1"),
             organisms="human",
             path="hicData/GSE84920", path_dest="compiledData"
           )
})
print(files)
#                human                                         
# GSM2254215_ML1 "compiledData/GSM2254215_ML1,human,unique.rds"
# GSM2254219_PL2 "compiledData/GSM2254219_PL2,human,unique.rds"
# GSM2254216_ML2 "compiledData/GSM2254216_ML2,human,unique.rds"
# GSM2254218_PL1 "compiledData/GSM2254218_PL1,human,unique.rds"
```

_Comment_: This step takes a few hours to complete.


## Step 3. Split Human Data by Cell Type

In R, call:

```r
files <- TopDomStudy::split_by_celltype(
           celltypes=list(human="HAP1"),
           path="compiledData")
print(files)
# $human
#                                 HAP1 
# "compiledData/human,HAP1,unique.rds"
```

_Comment_: This step takes approximately a minute to complete.



## Step 4. Split Human Data by Cell Type and Chromosomes

In R, call:

```r
files <- TopDomStudy::split_by_celltype_chromosome(
           celltypes=list(human="HAP1"),
	   chromosomes=1:25,
           path="compiledData")
str(files)
## List of 1
##  $ human:List of 1
##   ..$ HAP1: Named chr [1:25] "compiledData/human,HAP1,unique,chr=1.rds" "compiledData/human,HAP1,unique,chr=2.rds" "compiledData/human,HAP1,unique,chr=3.rds" "compiledData/human,HAP1,unique,chr=4.rds" ...
##   .. ..- attr(*, "names")= chr [1:25] "chr=1" "chr=2" "chr=3" "chr=4" ...
```

_Comment_: This step takes less than a minute to complete.

The `compiles/human,HAP1,unique,chr=*.rds` files correspond to the RDS files that are installed with this package in folder `system.file("compiledData", package="TomDopStudy")`.  The content of these files look like:
```r
> data <- readRDS("compiledData/human,HAP1,unique,chr=22.rds")
> tibble::as_tibble(data)
# A tibble: 112,704 x 9
   chr_a  start_a    end_a chr_b  start_b    end_b celltype cell_id           name          
   <chr>    <int>    <int> <chr>    <int>    <int> <chr>    <fct>             <chr>         
 1 22    16304723 16304853 22    16368550 16368588 HAP1     GGTCAGTG-TGTCTGCA GSM2254215_ML1
 2 22    16344591 16344666 22    17082891 17082926 HAP1     AAGCCGGT-CTACTAGG GSM2254215_ML1
 3 22    16357581 16357715 22    17723422 17723517 HAP1     TCGACTGC-TTAATCGA GSM2254219_PL2
 4 22    16433346 16433395 22    17060321 17060372 HAP1     ACCACCAC-TGTAATCG GSM2254216_ML2
 5 22    16433811 16433879 22    17137580 17137702 HAP1     TTGTGCCG-CGTTACTT GSM2254215_ML1
 6 22    16499667 16499748 22    17462757 17462829 HAP1     CGCGCAAT-CTTAGAAG GSM2254215_ML1
 7 22    16551301 16551348 22    21911741 21911808 HAP1     CGACATGG-CAGCATAT GSM2254215_ML1
 8 22    16554345 16554502 22    17900200 17900292 HAP1     AACGGTCG-TGCAGTGA GSM2254219_PL2
 9 22    16848715 16848855 22    16872125 16872172 HAP1     AAGCCGGT-AACGCGTA GSM2254216_ML2
10 22    16852229 16852468 22    16856842 16857008 HAP1     GCTGAGAC-CCTTATAG GSM2254215_ML1
# ... with 112,694 more rows
```


## References

1. Ramani, V., Deng, X., Qiu, R., Gunderson, K. L., Steemers, F. J., Disteche, C. M., … Shendure, J. (2017). Massively multiplex single-cell Hi-C. Nature methods, 14(3), 263–266. doi:10.1038/nmeth.4155, [PMC5330809](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5330809/)

[R]: https://www.r-project.org/
[GSE84920]: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84920
[GSM2254215]: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2254215
[GSM2254216]: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2254216
[GSM2254217]: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2254217
[GSM2254218]: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2254218
[GSM2254219]: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2254219
[GSM2438426]: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2438426
