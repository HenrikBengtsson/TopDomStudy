%\VignetteIndexEntry{TopDomStudy: Compiling Raw Data}
%\VignetteAuthor{Henrik Bengtsson}
%\VignetteEngine{TopDomStudy::self}
%\VignetteEncoding{UTF-8}
%\VignetteKeyword{R}
%\VignetteKeyword{package}
%\VignetteKeyword{vignette}

In this work, we study the human HAP1 cell types part of the Ramani et al. (2017) data set. For simplicity and memory performance, we focus on Chromosomes 12, 16, and 22.  The purpose of this document is to show how the three files:

 1. human,HAP1,unique,chr=12.rds
 2. human,HAP1,unique,chr=16.rds
 3. human,HAP1,unique,chr=22.rds

in the `system.file("compiledData", package="TomDopStudy")` folder, which are used for this study, were generated from the Ramani et al. (2017) data set.

The Ramani data set is published on NCBI's Gene Expression Omnibus (GEO) in the GEO series [GSE84920] \(titled 'Massively multiplex single-cell Hi-C'\), which contains:

| GEO Sample          | GEO Title                         | Cell Types
| ------------------- | --------------------------------- | ------------
| [GSM2254215]        | Combinatorial scHi-C Library ML1  | human ('HAP1', 'HeLa'), mouse ('MEF', 'Patski')
| [GSM2254216]        | Combinatorial scHi-C Library ML2  | human ('HAP1', 'HeLa'), mouse ('MEF', 'Patski')
| [GSM2254217]        | Combinatorial scHi-C Library ML3  | human ('GM12878', 'K562'), mouse ('MEF', 'Patski')
| [GSM2254218]        | Combinatorial scHi-C Library PL1  | human ('HAP1', 'HeLa'), mouse ('MEF', 'Patski')
| [GSM2254219]        | Combinatorial scHi-C Library PL2  | human ('HAP1', 'HeLa'), mouse ('MEF', 'Patski')
| [GSM2438426]**      | Combinatorial scHi-C Library ML4  | human ('Asynchronous', 'Nocadazole'), mouse ('Patski')



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
```sh
> stopifnot(file_test("-d", "hicData/GSE84920/")
> source(system.file("scripts", "compile_by_organism.R", package="TopDomStudy"))
```

## Step 2. Split Human Data by Cell Type

In R, call:
```sh
> stopifnot(file_test("-d", "hicData/GSE84920/")
> source(system.file("scripts", "split_by_celltype.R", package="TopDomStudy"))
```


## Step 3. Split Human Data by Cell Type and Chromosomes

In R, call:
```sh
> stopifnot(file_test("-d", "hicData/GSE84920/")
> source(system.file("scripts", "split_by_celltype_chromosome.R", package="TopDomStudy"))
```

This latter step will produce the three `human,HAP1,unique,chr=*.rds` files (located under hicData/GSE84920/).


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
