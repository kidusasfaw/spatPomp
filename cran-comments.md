
#### spatPomp 0.34.0

spatPomp 0.34.0 is an update of 0.33.0 with some minor updates and
bug fixes. It also has updates essential for compatibility with
pomp 5.5.0.0.

#### Test environments:

## devtools::check_win_devel()
Status: OK

## devtools::check_win_release()
Status: OK

## devtools::check_mac_release()

The following WARNING and the corresponding NOTEs did not seem to
require any action. The spatPomp package is not supposed to include
the source for the packages that it requires, and installing them
from CRAN is the intended action.

Building package dependency tree..
(from /Volumes/PkgBuild/work/1701616455-16fc4133d5f711ca/packages/CRAN/src/contrib)
./dtree --desc /Volumes/PkgBuild/work/1701616455-16fc4133d5f711ca/packages/CRAN/meta/src/contrib /Volumes/PkgBuild/work/1701616455-16fc4133d5f711ca/packages/CRAN/src/contrib > /Volumes/PkgBuild/work/1701616455-16fc4133d5f711ca/packages/CRAN/dep.list
NOTE: there is no source for package pomp [pomp], installed CRAN version will be used.
NOTE: there is no source for package testthat [testthat], installed CRAN version will be used.
NOTE: there is no source for package doParallel [doParallel], installed CRAN version will be used.
NOTE: there is no source for package parallel [parallel], installed CRAN version will be used.
NOTE: there is no source for package doRNG [doRNG], installed CRAN version will be used.
NOTE: there is no source for package foreach [foreach], installed CRAN version will be used.
NOTE: there is no source for package dplyr [dplyr], installed CRAN version will be used.
NOTE: there is no source for package tidyr [tidyr], installed CRAN version will be used.
NOTE: there is no source for package stringr [stringr], installed CRAN version will be used.
NOTE: there is no source for package abind [abind], installed CRAN version will be used.
NOTE: there is no source for package rlang [rlang], installed CRAN version will be used.
NOTE: there is no source for package magrittr [magrittr], installed CRAN version will be used.
NOTE: there is no source for package ggplot2 [ggplot2], installed CRAN version will be used.
NOTE: there is no source for package pomp [pomp], installed CRAN version will be used.

Status: 1 WARNING
  - installing from sources

## Rhub builder Debian Linux, R-release, GCC
Status: OK

## Reverse dependencies

none found by devtools::revdep("spatPomp")

#### Quality control

## Memory leak

rhub::check_with_valgrind() found a small memory leak which turned out to be due to png() rather than any spatPomp functionality. This was reported to r-devel@r-project.org, and still remains an issue in R4.3.2. The problem can be reproduced on Ubuntu 22.04.3 by

R -d "valgrind --tool=memcheck --track-origins=yes --leak-check=full" --vanilla -e "png(filename='p.png'); plot(1:10); dev.off()"
## HAS LEAK
==540949== LEAK SUMMARY:
==540949==    definitely lost: 8,448 bytes in 28 blocks
==540949==    indirectly lost: 12,424 bytes in 532 blocks
==540949==      possibly lost: 11,866 bytes in 328 blocks

R -d "valgrind --tool=memcheck --track-origins=yes --leak-check=full" --vanilla -e "pdf(file='p.pdf'); plot(1:10); dev.off()"
## NO LEAK
==541137== LEAK SUMMARY:
==541137==    definitely lost: 0 bytes in 0 blocks
==541137==    indirectly lost: 0 bytes in 0 blocks
==541137==      possibly lost: 0 bytes in 0 blocks

## Unit tests

Unit test code coverage is fairly constant in recent releases, at about 90%
covr::package_coverage(type="tests") was run in R4.3.2 on Ubuntu 22.04.3
0.34.0   coverage: 89.76%
0.33.0   coverage: 89.64%
0.32.0   coverage: 89.43%
0.31.0.0 coverage: 88.9% 

