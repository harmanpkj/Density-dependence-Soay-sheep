# Convert the survival, s, output of roedeer matrices list() into pvec as a matrix.
# Values for Roe deer only, but easily modified to any age+stage IPM.

s_to_pvec <- function(matrices){
  # there are 12 age group c(1:11, >11), but 4 age based survival
  # Four age classes: yearlings, 2-7 years old, 8-11 years old and >11 years
  pvec <- cbind(diag(matrices[[1]]$S), #yearlings
                diag(matrices[[2]]$S), diag(matrices[[2]]$S), diag(matrices[[2]]$S), diag(matrices[[2]]$S), diag(matrices[[2]]$S), diag(matrices[[2]]$S),#2-7 years old
                diag(matrices[[3]]$S), diag(matrices[[3]]$S), diag(matrices[[3]]$S), diag(matrices[[3]]$S), #8-11 years old
                diag(matrices[[4]]$S)) # >11 years
  return(pvec)
}
