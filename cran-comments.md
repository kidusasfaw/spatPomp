
#### spatPomp 0.35.0

spatPomp 0.35.0 has various technical updates without significant changes for the user. Test coverage has increased to 100%.

#### Test environments:

## devtools::check_win_devel()
Status: OK

## devtools::check_win_release()
Status: OK

## devtools::check_mac_release()
This gave an error message:
  Internal Server Error (HTTP 500). Failed to Uploading package.
Instead, spatPomp was checked with the current version of R on
an M1 and x86 Mac

Mac Sonoma 14.5 on Apple M1 Max with R4.4.0
R CMD check --as-cran spatPomp_0.35.0.tar.gz
Status: OK

Mac Sonoma 14.4.1 on x86 with R4.4.0
R CMD check --as-cran spatPomp_0.35.0.tar.gz
Status: OK

## R4.4.0 on Ubuntu 22.04.4
Status: OK

## Reverse dependencies

none found by devtools::revdep("spatPomp")

#### Quality control

## Memory leak
R4.4.0 on Ubuntu 22.04.4

R --vanilla -d "valgrind --tool=memcheck --track-origins=yes --leak-check=full" < tests/bm.R 2>&1 | tee valgrind-bm.Rout

This identifies a memory leak. As identified in previous CRAN uploads of spatPomp, the problem seems to be in png.

R -d "valgrind --tool=memcheck --track-origins=yes --leak-check=full" --vanilla -e "png(filename='p.png'); plot(1:10); dev.off()"
## HAS LEAK
==4041830== LEAK SUMMARY:
==4041830==    definitely lost: 9,984 bytes in 31 blocks
==4041830==    indirectly lost: 20,273 bytes in 864 blocks
==4041830==      possibly lost: 4,813 bytes in 51 blocks

R -d "valgrind --tool=memcheck --track-origins=yes --leak-check=full" --vanilla -e "pdf(file='p.pdf'); plot(1:10); dev.off()"

## NO LEAK
==4041953== LEAK SUMMARY:
==4041953==    definitely lost: 0 bytes in 0 blocks
==4041953==    indirectly lost: 0 bytes in 0 blocks
==4041953==      possibly lost: 0 bytes in 0 blocks

## Unit tests

Unit test code coverage is fairly constant in recent releases, at about 90%
covr::package_coverage(type="tests") was run in R4.4.0 on intel Mac (Sonoma 14.4.1)
coverage: 100%
