# Calculate variance of an arbitrary (discrete) probability distribution.

distr_var <- function(distr){
  ##distr is a data frame which has density and mids
  output <- sum(distr$density*(distr$mids - distr_mean(distr))^2)
  return(output)
}
