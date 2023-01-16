
#### spatPomp 0.31.0.0 

spatPomp 0.31.0.0 is an update of 0.30.2.0 necessitated by changes in pomp 4.6 as well as including 
various minor updates and bug fixes. 

#### R CMD check --as-cran

one NOTE:
* checking package dependencies ... NOTE
Suggests orphaned package: ‘doRNG’

Version 1.8.3 of doRNG was posted to CRAN on 2022-12-19 so perhaps this package is about to be restored 
to good standing? 60 packages show up with devtools::revdep("doRNG") so this is an issue than extends 
beyond spatPomp.

## Reverse dependencies

none found by devtools::revdep("spatPomp")

#### Quality control

rhub::check_with_valgrind() found a small memory leak which turned out to be due to png() rather than any spatPomp functionality. This was reported to r-devel@r-project.org. The problem can be reproduced by

R -d "valgrind --tool=memcheck --track-origins=yes --leak-check=full" --vanilla -e "png(filename='p.png'); plot(1:10); dev.off()"
## HAS LEAK
==1021711== LEAK SUMMARY:
==1021711==    definitely lost: 9,216 bytes in 30 blocks
==1021711==    indirectly lost: 19,370 bytes in 838 blocks
==1021711==      possibly lost: 3,868 bytes in 8 blocks

R -d "valgrind --tool=memcheck --track-origins=yes --leak-check=full" --vanilla -e "pdf(file='p.pdf'); plot(1:10); dev.off()"
## NO LEAK
==1031300== LEAK SUMMARY:
==1031300==    definitely lost: 0 bytes in 0 blocks
==1031300==    indirectly lost: 0 bytes in 0 blocks
==1031300==      possibly lost: 0 bytes in 0 blocks

covr::package_coverage(type="tests") shows 88.9% coverage, up from 86.7% in 0.30.2.0. Further increases in coverage are a future development goal.

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
present for 0.30.0.1 and 0.30.2.0

3.   Platform:   Ubuntu Linux 20.04.1 LTS, R-release, GCC

NOTE:
  Found the following (possibly) invalid DOIs:
    DOI: 10.1029/94JC00572
      From: DESCRIPTION
      Status: Service Unavailable
      Message: 503

This DOIs checks out as okay. The same note arose for 0.30.0.1 and 0.30.2.0

4.   Platform:   Debian Linux, R-devel, GCC ASAN/UBSAN

OK. No notes.



