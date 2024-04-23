options(digits=3)

library(spatPomp)

try(spatPomp:::conc())
try(spatPomp:::conc("a","b"))

sp1 <- bm2(U = 2, N = 20, shared_names = "rho")
sp2 <- bm2(U = 2, N = 20, shared_names = "rho")
sp3 <- bm2(U = 3, N = 10, shared_names = "rho")

(t_spatPompList <- is(c(sp1, sp2), "spatPompList"))

class( c(c(sp1,sp2),sp3))

(t_bpfilterList <- is(
  c(
    bpfilter(sp1, Np = 10, block_size = 1),
    bpfilter(sp2, Np = 10, block_size = 1)
  ),
  "bpfilterList"
))

## ibpfilterList is tested in He10

stopifnot(all(t_spatPompList, t_bpfilterList))

