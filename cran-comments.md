## R CMD CHECK --as-cran results

one NOTE:
Maintainer: ‘Edward Ionides <ionides@umich.edu>’

New maintainer:
  Edward Ionides <ionides@umich.edu>
Old maintainer(s):
  Kidus Asfaw <kasfaw@umich.edu>


Kidus will email from kidusasfaw1990@gmail.com to confirm the transfer of the maintainer role.

## Reverse dependencies

none

## Quality control

rhub::check_with_valgrind() found no memory leaks

covr::package_coverage(type="tests") shows 83.1% coverage. The main functionality is tested, but not yet all edge cases. TODO.md identifies further increases in coverage as a future development goal.

## Test environments: rhub::check_for_cran()

Unresolved notes other than change of maintainer:

1.  Platform:   Windows Server 2022, R-devel, 64 bit

NOTE:
  Found the following files/directories:
    'lastMiKTeXException'

* No explanation. Let me know if this needs attention *

2.   Platform:   Fedora Linux, R-devel, clang, gfortran

NOTE:
  Skipping checking HTML validation: no command 'tidy' found
  Skipping checking math rendering: package 'V8' unavailable


3.   Platform:   Ubuntu Linux 20.04.1 LTS, R-release, GCC

NOTE:
  Found the following (possibly) invalid DOIs:
    DOI: 10.1029/94JC00572
      From: DESCRIPTION
      Status: Service Unavailable
      Message: 503
    DOI: 10.1214/14-AAP1061
      From: DESCRIPTION
      Status: Internal Server Error
      Message: 500

* These DOIs check out as okay *


4.   Platform:   Debian Linux, R-devel, GCC ASAN/UBSAN

no notes



