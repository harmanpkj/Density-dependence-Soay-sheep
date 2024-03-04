# Convert fertility to number of offspring distribution by binomial distribution.
# It can be replaced if any other number of offspring distribution is used.

fbern_to_kappamat <- function(Fmatrix, max.offspring = 12){
  # Since we only consider one state for offspring, kappamat is a matrix here
  # otherwise it should be an array
  n.stages <- dim(Fmatrix)[2]
  kappamat <- matrix(0, nrow = max.offspring + 1, ncol = n.stages)
  kappamat[1, ] <- 1 - Fmatrix[1, ]
  kappamat[2, ] <- Fmatrix[1, ]
  return(kappamat)
}
