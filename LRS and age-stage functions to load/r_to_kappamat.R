# Convert the reproduction, r, output of roedeer matrices list() into kappamat as a list.
# Uses Bionomial reproduction. Values for Roe deer only, but easily modified to any age+stage IPM.


r_to_kappamat <- function(matrices){
  # there are 12 age group c(1:11, >11), but 4 age based reproduction
  # Four age classes: yearlings, 2-7 years old, 8-11 years old and >11 years
  kappamat <- list(matrices[[1]]$R, #yearlings
                   matrices[[2]]$R, matrices[[2]]$R, matrices[[2]]$R, matrices[[2]]$R, matrices[[2]]$R, matrices[[2]]$R,#2-7 years old
                   matrices[[3]]$R, matrices[[3]]$R, matrices[[3]]$R, matrices[[3]]$R, #8-11 years old
                   matrices[[4]]$R) # >11 years
  return(kappamat)
}
