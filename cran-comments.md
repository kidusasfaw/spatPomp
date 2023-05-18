
#### spatPomp 0.32.0

spatPomp 0.32.0 is an update of 0.31.0.0 with various minor
updates and bug fixes. 

#### Test environments:

## devtools::check_win_devel() 
check: OK

## devtools::check_win_release() 
check: OK

## devtools::check_mac_release()
one test discrepancy:

Comparing ‘gbm.Rout’ to ‘gbm.Rout.save’ ...46c46
< [1] "ARMA benchmark: 206.2569593745"
---
> [1] "ARMA benchmark: 196.3173327754"

on the ARMA likelihood evaluation:

arma_benchmark(gbm_model,order=c(1,0,0))
paste("ARMA benchmark:", round(a1$total,10))

using R version 4.3.0 beta (2023-04-11 r84222)
* using platform: aarch64-apple-darwin20 (64-bit)
* R was compiled by
    Apple clang version 14.0.0 (clang-1400.0.29.202)
    GNU Fortran (GCC) 12.2.0
* running under: macOS Ventura 13.3.1

Since this problem did not recur in other tests, including r_devel and a
subsequent R 4.3.0, this problem is being ignored for now.

## R CMD check --as-cran in R4.3.0 on Ubuntu 22.04.2
check: OK

## R CMD check --as-cram in R4.2.3 on macOS Ventura 13.3.1(a)
check: OK

## rhub::check_for_cran()

This stalled for some unknown reason, running R4.2.3 on MacOS Ventura 13.3.1(a).
After re-registing with rhub, there was nothing else easy to try so testing was 
done with devtools test environments (for Windows and Mac) as well as local Mac
and Linux machines.

## Reverse dependencies

none found by devtools::revdep("spatPomp")

#### Quality control

rhub::check_with_valgrind() found a small memory leak which turned out to be due to png() rather than any spatPomp functionality. This was reported to r-devel@r-project.org, and still remains an issue in R4.3.0. The problem can be reproduced by

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


covr::package_coverage(type="tests") was run in R4.3.0 on Ubuntu 22.04.2
shows 89.43% coverage, a small increase from 88.9% in 0.31.0.0 (and 86.7% in 0.30.2.0). Further increases in coverage are a future development goal.



