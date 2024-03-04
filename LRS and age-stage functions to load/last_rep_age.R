# Obtain/set the last reproductive age.

last_rep_age <- function(kappamat){
  last_reproductive_age <- max(which(kappamat[1,] != 1))
  return(last_reproductive_age)
}
