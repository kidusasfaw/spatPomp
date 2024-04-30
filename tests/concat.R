options(digits=3)

library(spatPomp)

try(spatPomp:::conc())
try(spatPomp:::conc("a","b"))

sp1 <- bm(U = 2, N = 4)
sp2 <- sp1
sp3 <- bm(U = 3, N = 3)

(t_spatPompList <- is(c(sp1, sp2), "spatPompList"))

class( c(c(sp1,sp2),sp3))

(t_bpfilterList <- is(
  c(
    bpfilter(sp1, Np = 5, block_size = 1),
    bpfilter(sp2, Np = 5, block_size = 1)
  ),
  "bpfilterList"
))

## ibpfilterList is tested in He10

stopifnot(all(t_spatPompList, t_bpfilterList))

