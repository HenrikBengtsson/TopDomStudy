PKGNAME=TopDomStudy

build_fast:
	(cd ..; R CMD build --no-resave-data "$(PKGNAME)")

build:
	(cd ..; R CMD build "$(PKGNAME)")

install:
	(cd ..; R CMD INSTALL "$(PKGNAME)"_*.tar.gz)

check:
	(cd ..; R CMD check "$(PKGNAME)"_*.tar.gz)

spelling:
	Rscript -e "spelling::spell_check_files('NEWS', ignore=readLines('inst/WORDLIST'))"
	Rscript -e "spelling::spell_check_package()"
