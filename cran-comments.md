
#### spatPomp 0.33.0

spatPomp 0.33.0 is an update of 0.33.0 with some minor
updates and bug fixes, and some new methods for ibpf and bpfilter objects.

This release also fixes an issue with M1 mac which led to spatPomp being
temporarily pulled from CRAN. The issue involved inconsistent results from
stats::arima on M1 compared to other platforms. Specifically, running
arima on the same data on other platforms gave a log-likelihood of
206.2569593745 compared to 196.3173327754. This issue has been avoided
in 0.33.0 by omitting a unit test for spatPomp::arima_benchmark.

#### Test environments:

## devtools::check_win_devel()

One note:
Maintainer: 'Edward Ionides <ionides@umich.edu>'
New submission
Package was archived on CRAN

## devtools::check_win_release()

One note:
Maintainer: 'Edward Ionides <ionides@umich.edu>'
New submission
Package was archived on CRAN

## devtools::check_mac_release()

Status: OK

Various warnings were produced, which did not result in NOTES.
The reason for these warnings is diagnosed as an issue with the
pomp package, resolved by the pomp development version,
commit c42cef25b2013ecfddfaf6838bee78a00dcc5ccc at
https://github.com/kingaa/pomp, which concerns an M1 mac
treatment of null objects returned by .Call. 

## R CMD check --as-cran in R4.3.1 on Ubuntu 22.04.2

One note:
Maintainer: ‘Edward Ionides <ionides@umich.edu>’
New submission
Package was archived on CRAN
CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2023-06-09 as issues were not corrected
    in time.

## Reverse dependencies

none found by devtools::revdep("spatPomp")

#### Quality control

rhub::check_with_valgrind() found a small memory leak which turned out to be due to png() rather than any spatPomp functionality. This was reported to r-devel@r-project.org, and still remains an issue in R4.3.1. The problem can be reproduced on Ubuntu 22.04.2 by

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

Unit test code coverage is fairly constant in recent releases, at about 90%
covr::package_coverage(type="tests") was run in R4.3.1 on Ubuntu 22.04.2
0.33.0   coverage: 89.64%
0.32.0   coverage: 89.43%
0.31.0.0 coverage: 88.9% 

