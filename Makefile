REXE = R --vanilla
RCMD = $(REXE) CMD
RCMD_ALT = R --no-save --no-restore CMD
RSCRIPT = Rscript --vanilla

PDFLATEX = pdflatex
BIBTEX = bibtex
MAKEIDX = makeindex

RM = rm -f
CP = cp
TOUCH = touch
INSTALL = install

PKG = $(shell perl -ne 'print $$1 if /Package:\s+((\w+[-\.]?)+)/;' DESCRIPTION)
VERSION = $(shell perl -ne 'print $$1 if /Version:\s+((\d+[-\.]?)+)/;' DESCRIPTION)
PKGVERS = $(PKG)_$(VERSION)
SOURCE=$(shell ls R/*R src/*.c src/*.h data/*)
CSOURCE=$(shell ls src/*.c)
TESTS=$(shell ls tests/*R)
REVDEPS=spaero epimdr CollocInfer

default:
	@echo $(PKGVERS)

.PHONY: clean win wind tests check

dist manual vignettes: export R_QPDF=qpdf
headers: export LC_COLLATE=C
roxy headers dist manual vignettes: export R_HOME=$(shell $(REXE) RHOME)
check xcheck xxcheck: export FULL_TESTS=yes
dist revdeps session tests check xcheck xxcheck: export R_KEEP_PKG_SOURCE=yes
revdeps xcheck tests: export R_PROFILE_USER=$(CURDIR)/.Rprofile
revdeps session xxcheck htmldocs vignettes data tests manual: export R_LIBS=$(CURDIR)/library
session: export R_DEFAULT_PACKAGES=datasets,utils,grDevices,graphics,stats,methods,pomp,tidyverse

includes: inst/include/spatPomp_defines.h

inst/include/%.h: src/%.h
	$(CP) $^ $@

NEWS: inst/NEWS

inst/NEWS: inst/NEWS.Rd
	$(RCMD) Rdconv -t txt $^ -o $@

dist: NEWS $(PKGVERS).tar.gz

$(PKGVERS).tar.gz: $(SOURCE) $(TESTS) includes
	$(RCMD) build --force --no-manual --resave-data --compact-vignettes=both --md5 .

win: dist
	curl -T $(PKGVERS).tar.gz ftp://win-builder.r-project.org/R-release/

wind: dist
	curl -T $(PKGVERS).tar.gz ftp://win-builder.r-project.org/R-devel/

xcheck: dist
	mkdir -p check library
	$(RCMD_ALT) check --no-stop-on-test-error --as-cran --library=library -o check $(PKGVERS).tar.gz
