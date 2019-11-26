PKGNAME=TopDomStudy

build_fast:
	(cd ..; R CMD build --no-resave-data "$(PKGNAME)")

build:
	(cd ..; R CMD build "$(PKGNAME)")

install: build
	(cd ..; R CMD INSTALL "$(PKGNAME)_*.tar.gz")

install_fast: build_fast
	(cd ..; R CMD INSTALL "$(PKGNAME)_*.tar.gz")
