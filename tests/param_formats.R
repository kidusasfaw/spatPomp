
library(spatPomp)

p_expanded <- c(a1=0,b1=0,b2=1,b3=2,c1=4,c2=4,c3=4)

p_contracted <- contract_params(p_expanded,expandedParNames="c",U=3)
  
p_contracted

p_expanded2 <- expand_params(p_contracted,expandedParNames="c",U=3)

p_expanded2

if(any(p_expanded[names(p_expanded2)]!=p_expanded2)) stop(
  "failed inverse for expand_params() and contract_params()"
)

mean_by_unit(p_expanded,expandedParNames=c("b","c"),U=3)

mean_by_unit(p_expanded,expandedParNames=c("c"),U=3)

