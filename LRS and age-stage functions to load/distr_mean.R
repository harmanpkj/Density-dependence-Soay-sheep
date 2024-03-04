# Calculate mean of an arbitrary (discrete) probability distribution

distr_mean <- function(distr){
  ##distr is a data frame which has density and mids
  output <- sum(distr$density * distr$mids)
  return(output)
}
