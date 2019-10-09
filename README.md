# TopDomStudy: TopDom Study by Segal et al.


## Installation
R package TopDomStudy is only available via [GitHub](https://github.com/HenrikBengtsson/TopDomStudy) and can be installed in R as:
```r
remotes::install_github("HenrikBengtsson/TopDomStudy")
```

### Pre-release version

To install the pre-release version that is available in Git branch `develop` on GitHub, use:
```r
remotes::install_github("HenrikBengtsson/TopDomStudy@develop")
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
