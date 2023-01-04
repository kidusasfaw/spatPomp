
#### spatPomp 0.30.2.0 

spatPomp 0.30.2.0 is a routine update of 0.30.0.1 with various minor updates and bug fixes. 

#### R CMD CHECK --as-cran

one NOTE:
* checking package dependencies ... NOTE
Suggests orphaned package: ‘doRNG’

Version 1.8.3 of doRNG was posted to CRAN on 2022-12-19 so perhaps this package is about to be restored to good standing?

## Reverse dependencies

none found by devtools::revdep("spatPomp")

#### Quality control

CHECK  rhub::check_with_valgrind() found no memory leaks

covr::package_coverage(type="tests") shows 86.7% coverage, up from 83.1% in 0.30.0.1. Further increases in coverage are a future development goal.

#### Test environments:

## devtools::check_win_devel() for win-builder (win-builder.r-project.org)

## rhub::check_for_cran()

1.  Platform:   Windows Server 2022, R-devel, 64 bit

OK. No notes.

2.   Platform:   Fedora Linux, R-devel, clang, gfortran

In addition to the orphaned doRNG discussed above:

NOTE:
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found
Skipping checking math rendering: package 'V8' unavailable

Please let me know if this needs attention. The same note was
present for 0.30.0.1

3.   Platform:   Ubuntu Linux 20.04.1 LTS, R-release, GCC

NOTE:
  Found the following (possibly) invalid DOIs:
    DOI: 10.1029/94JC00572
      From: DESCRIPTION
      Status: Service Unavailable
      Message: 503

This DOIs checks out as okay. The same note arose for 0.30.0.1

4.   Platform:   Debian Linux, R-devel, GCC ASAN/UBSAN

OK. No notes.



