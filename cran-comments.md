
#### Resubmission comments

spatPomp 0.30.0.1 is a resubmission of 0.30.0.0 ammended to fix an error
arising on win-builder. The previous submission was tested only on
rhub::check_for_cran() and apparently there is a relevant difference
in treatment of foreach between
rhub::check(path=".",platform="windows-x86_64-devel") and
devtools::check_win_devel(). Running foreach in series by omitting the
registration of a parallel backend is okay for both. Running the nominally
parallel code by registering a 1-core parallel backend led to an error on
win-builder.

#### R CMD CHECK --as-cran results

one NOTE:
Maintainer: ‘Edward Ionides <ionides@umich.edu>’

New maintainer:
  Edward Ionides <ionides@umich.edu>
Old maintainer(s):
  Kidus Asfaw <kasfaw@umich.edu>


Kidus emailed to confirm the transfer of the maintainer role.

## Reverse dependencies

none

#### Quality control

rhub::check_with_valgrind() found no memory leaks

covr::package_coverage(type="tests") shows 83.1% coverage. The main functionality is tested, but not yet all edge cases. TODO.md identifies further increases in coverage as a future development goal.

#### Test environments:

## devtools::check_win_devel()

No unresolved notes other than change of maintainer

## rhub::check_for_cran()

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



