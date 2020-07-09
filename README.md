# TopDomStudy: TopDom Study by Segal et al.


## Requirements

This [R] package depends on several R packages.  Most of them are available on [CRAN] but not all.  Currently, the following packages are available only on GitHub:

* **[TopDom]** - An efficient and Deterministic Method for identifying Topological Domains in Genomes
* **[GSE84920.parser]** - Parsers for the GSE84920 Hi-C Data Sets of Ramani et al. (2017)

However, you do _not_ need to install those manually as long as you follow the installation instruction below; they will be install automatically.

Moreover, this packages produces and makes use of pathnames that are longer than 255 characters.  Unfortunately, this means that this package _does not run on MS Windows_ (https://github.com/HenrikBengtsson/TopDomStudy/issues/3).  The package tests have been validated on Linux and macOS.


## References

1. Ramani, V., Deng, X., Qiu, R., Gunderson, K. L., Steemers, F. J., Disteche,
   C. M., ..., Shendure, J. (2017). Massively multiplex single-cell Hi-C.
   Nature methods, 14(3), 263–266.
   doi:10.1038/nmeth.4155,
   [PMC5330809](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5330809/).

2. Shin, et al., TopDom: an efficient and deterministic method for
   identifying topological domains in genomes,
   Nucleic Acids Res. 2016 Apr 20; 44(7): e70., 2016.
   doi: 10.1093/nar/gkv1505,
   PMCID: [PMC4838359](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4838359/),
   PMID: 26704975.


## Installation

R package **TopDomStudy** is only available via [GitHub](https://github.com/HenrikBengtsson/TopDomStudy) and can be installed in R as:
```r
remotes::install_github("HenrikBengtsson/TopDomStudy", ref="master")
```

### Pre-release version

To install the pre-release version that is available in Git branch `develop` on GitHub, use:
```r
remotes::install_github("HenrikBengtsson/TopDomStudy", ref="develop")
```
This will install the package from source.  



## Contributions

This Git repository uses the [Git Flow](http://nvie.com/posts/a-successful-git-branching-model/) branching model (the [`git flow`](https://github.com/petervanderdoes/gitflow-avh) extension is useful for this).  The [`develop`](https://github.com/HenrikBengtsson/TopDomStudy/tree/develop) branch contains the latest contributions and other code that will appear in the next release, and the [`master`](https://github.com/HenrikBengtsson/TopDomStudy) branch contains the code of the latest release.

Contributing to this package is easy.  Just send a [pull request](https://help.github.com/articles/using-pull-requests/).  When you send your PR, make sure `develop` is the destination branch on the [TopDomStudy repository](https://github.com/HenrikBengtsson/TopDomStudy).  Your PR should pass `R CMD check --as-cran`, which will also be checked by <a href="https://travis-ci.org/HenrikBengtsson/TopDomStudy">Travis CI</a> and <a href="https://ci.appveyor.com/project/HenrikBengtsson/TopDomStudy">AppVeyor CI</a> when the PR is submitted.


## Software status

| Resource:     | GitHub        | Travis CI       | AppVeyor         |
| ------------- | ------------------- | --------------- | ---------------- |
| _Platforms:_  | _Multiple_          | _Linux & macOS_ | _Windows_        |
| R CMD check   |  | <a href="https://travis-ci.org/HenrikBengtsson/TopDomStudy"><img src="https://travis-ci.org/HenrikBengtsson/TopDomStudy.svg" alt="Build status"></a>   | <a href="https://ci.appveyor.com/project/HenrikBengtsson/TopDomStudy"><img src="https://ci.appveyor.com/api/projects/status/github/HenrikBengtsson/TopDomStudy?svg=true" alt="Build status"></a> |
| Test coverage |                     | <a href="https://codecov.io/gh/HenrikBengtsson/TopDomStudy"><img src="https://codecov.io/gh/HenrikBengtsson/TopDomStudy/branch/develop/graph/badge.svg" alt="Coverage Status"/></a>     |                  |



[R]: https://www.r-project.org
[CRAN]: https://cran.r-project.org
[Bioconductor]: https://www.bioconductor.org
[TopDom]: https://github.com/HenrikBengtsson/TopDom
[GSE84920.parser]: https://github.com/HenrikBengtsson/GSE84920.parser
[progressr]: https://github.com/HenrikBengtsson/progressr

