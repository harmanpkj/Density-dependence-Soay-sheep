# Convert the transition matrix, g, output of roedeer matrices list() into G as a list.
# Values for Roe deer only, but easily modified to any age+stage IPM.


g_to_gmat <- function(matrices){
  # there are 12 age group c(1:11, >11), but 4 age based transition
  # Four age classes: yearlings, 2-7 years old, 8-11 years old and >11 years
  G <- list(matrices[[1]]$G, #yearlings
            matrices[[2]]$G, matrices[[2]]$G, matrices[[2]]$G, matrices[[2]]$G, matrices[[2]]$G, matrices[[2]]$G,#2-7 years old
            matrices[[3]]$G, matrices[[3]]$G, matrices[[3]]$G, matrices[[3]]$G, #8-11 years old
            matrices[[4]]$G) # >11 years
  return(G)
}
