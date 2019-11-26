%\VignetteIndexEntry{TopDomStudy: Data Preprocessing}
%\VignetteAuthor{Henrik Bengtsson}
%\VignetteEngine{TopDomStudy::self}
%\VignetteEncoding{UTF-8}
%\VignetteKeyword{R}
%\VignetteKeyword{package}
%\VignetteKeyword{vignette}

In this work, we study the human HAP1 cell types part of the Ramani et al. (2017) data set. In addition to provide a set of R functions to help reproduce the results in our study, this package also provides a set of pre-processed data files to simply reproducibility:

<!--
path <- system.file("compiledData", package="TomDopStudy")
path <- "inst/compiledData"
filenames <- dir(path = path, pattern = ",chr=.+[.]rds$")
pathnames <- file.path(path, filenames)
sizes <- sapply(pathnames, FUN = file.size)
nreads <- sapply(pathnames, FUN = function(f) nrow(readRDS(f)))
chrs <- sub(".*,chr=(.+)[.]rds$", "\\1", filenames)
files <- data.frame(chromosome = chrs, size = sizes, reads = nreads, stringsAsFactors = FALSE, filename = filenames)
files <- files[gtools::mixedorder(files$chromosome),]
rownames(files) <- NULL
-->

<!--
tbl <- files[,c("chromosome", "size", "reads")]
tbl$chromosome <- sprintf("Chr %s", tbl$chromosome)
tbl <- rbind(tbl, list("total", sum(tbl$size), sum(tbl$reads)))
colnames(tbl) <- c("Chromosome", "File size (bytes)", "Number of Read Pairs")
knitr::kable(tbl, format.args = list(big.mark = ','))
-->

|Chromosome | File size (bytes)| Number of Read Pairs|
|:----------|-----------------:|--------------------:|
|Chr 1      |        10,897,300|              898,587|
|Chr 2      |        11,749,566|              972,893|
|Chr 3      |         9,790,836|              812,916|
|Chr 4      |         8,987,704|              744,741|
|Chr 5      |         8,446,368|              696,648|
|Chr 6      |         8,005,164|              662,910|
|Chr 7      |         7,183,925|              591,839|
|Chr 8      |         6,909,628|              569,888|
|Chr 9      |         5,056,207|              417,634|
|Chr 10     |         6,207,292|              508,691|
|Chr 11     |         6,366,477|              524,712|
|Chr 12     |         6,498,021|              537,659|
|Chr 13     |         4,516,108|              372,950|
|Chr 14     |         4,206,679|              345,941|
|Chr 15     |         4,599,528|              400,169|
|Chr 16     |         3,302,761|              265,621|
|Chr 17     |         3,390,448|              274,831|
|Chr 18     |         3,526,220|              290,220|
|Chr 19     |         2,342,704|              186,474|
|Chr 20     |         2,846,595|              232,448|
|Chr 21     |         1,466,148|              118,787|
|Chr 22     |         1,397,216|              112,704|
|total      |       127,692,895|           10,539,263|

This document describes how these data files where produced.

in the `system.file("compiledData", package="TomDopStudy")` folder, which are used for this study, were generated from the Ramani et al. (2017) data set.

The Ramani data set is published on NCBI's Gene Expression Omnibus (GEO) in the GEO series [GSE84920] \(titled 'Massively multiplex single-cell Hi-C'\), which contains:

| GEO Sample      | GEO Title                        | Cell Types
| --------------- | -------------------------------- | ------------------------------------------------------
| [GSM2254215]    | Combinatorial scHi-C Library ML1 | human ('HAP1', 'HeLa'), mouse ('MEF', 'Patski')
| [GSM2254216]    | Combinatorial scHi-C Library ML2 | human ('HAP1', 'HeLa'), mouse ('MEF', 'Patski')
| [GSM2254217]    | Combinatorial scHi-C Library ML3 | human ('GM12878', 'K562'), mouse ('MEF', 'Patski')
| [GSM2254218]    | Combinatorial scHi-C Library PL1 | human ('HAP1', 'HeLa'), mouse ('MEF', 'Patski')
| [GSM2254219]    | Combinatorial scHi-C Library PL2 | human ('HAP1', 'HeLa'), mouse ('MEF', 'Patski')
| [GSM2438426](*) | Combinatorial scHi-C Library ML4 | human ('Asynchronous', 'Nocadazole'), mouse ('Patski')



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
	   chromosomes=1:22,
           path="compiledData")
str(files)
## List of 1
##  $ human:List of 1
##   ..$ HAP1: Named chr [1:22] "compiledData/human,HAP1,unique,chr=1.rds" "compiledData/human,HAP1,unique,chr=2.rds" "compiledData/human,HAP1,unique,chr=3.rds" "compiledData/human,HAP1,unique,chr=4.rds" ...
##   .. ..- attr(*, "names")= chr [1:22] "chr=1" "chr=2" "chr=3" "chr=4" ...
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
