
#### spatPomp 0.34.2

spatPomp 0.34.2 is an update of 0.34.0 with some minor updates and
bug fixes. 

#### Test environments:

## devtools::check_win_devel()
Status: OK

## devtools::check_win_release()
Status: OK

## devtools::check_mac_release()

The following WARNING and the corresponding NOTEs did not seem to
require any action. The spatPomp package is not supposed to include
the source for the packages that it requires, and installing them
from CRAN is the intended action. The same WARNING was reported
in the 0.34.0 submission.

Building package dependency tree..
(from /Volumes/PkgBuild/work/1708265385-8ca2acc1b1297bab/packages/CRAN/src/contrib)
./dtree --desc /Volumes/PkgBuild/work/1708265385-8ca2acc1b1297bab/packages/CRAN/meta/src/contrib /Volumes/PkgBuild/work/1708265385-8ca2acc1b1297bab/packages/CRAN/src/contrib > /Volumes/PkgBuild/work/1708265385-8ca2acc1b1297bab/packages/CRAN/dep.list
NOTE: there is no source for package pomp [pomp], installed CRAN version will be used.
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
Checking all packages..
1708265404:20240219:031004:spatPomp:
  /Volumes/PkgBuild/work/1708265385-8ca2acc1b1297bab/packages/CRAN/src/contrib/spatPomp_0.34.2.tar.gz
 - looking for /Volumes/PkgBuild/work/1708265385-8ca2acc1b1297bab/packages/CRAN/meta/src/contrib/spatPomp_0.34.2.DESCRIPTION
 - using cached DESCRIPTION file for detection
  - unpacking
  - checking spatPomp 
* checking whether package ‘spatPomp’ can be installed ... [5s/6s] WARNING

Status: 1 WARNING
  - installing from sources

## R4.3.2 on Ubuntu 22.04.3

One WARNING from
R --no-save --no-restore CMD check --no-stop-on-test-error --as-cran --library=library -o check spatPomp_0.34.2.tar.gz

Warning: invalid uid value replaced by that for user 'nobody'

This is apparently a known issue and not a serious problem:
https://stackoverflow.com/questions/30599326/warning-message-during-building-an-r-package-invalid-uid-value-replaced-by-that

## Reverse dependencies

none found by devtools::revdep("spatPomp")

#### Quality control

## Memory leak

rhub::check_with_valgrind() found a small memory leak which turned out to be due to png() rather than any spatPomp functionality. This was reported to r-devel@r-project.org on 1/16/23, and still remains an issue in R4.3.2. The problem can be reproduced on Ubuntu 22.04.3 by

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
0.34.2   coverate: 90.17%
0.34.0   coverage: 89.76%
0.33.0   coverage: 89.64%
0.32.0   coverage: 89.43%
0.31.0.0 coverage: 88.9% 

