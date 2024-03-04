# Calculate the lifetime reproductive success distribution for age only models.


age_only_repro_distr <- function(Fmatrix, Umatrix,
                                 max.kids = 1,
                                 distri_type = "bernoulli"){
  #### Obtain the lrs distribution for age specific fertilities only
  ##======End at the age of last reproduction
  if(distri_type == "poisson"){
    kappamat <- fpois_to_kappamat(Fmatrix, max.kids)# for age only input max.kids(/birth) but stages using max.offspring(/lifetime)
  }else if(distri_type == "bernoulli"){
    kappamat <- fbern_to_kappamat(Fmatrix, max.kids)# for age only input max.kids(/birth) but stages using max.offspring(/lifetime)
  }
  pvec <- c(diag(Umatrix[-1,])) # survival at age i
  omega <- dim(kappamat)[2] # maximum age
  max.offspring <- omega*max.kids # maximum kids per life time
  lvec <- c(1, cumprod(pvec))
  lvec <- lvec[-length(lvec)] # survival to age i
  phi <- c(lvec * (1 - pvec), 0) # death distribution
  last_reproductive_age <- last_rep_age(kappamat)
  gamma <- matrix(0, nrow=max.offspring, ncol=last_reproductive_age)
  for(jj in 1:last_reproductive_age){
    if(jj==1){
      z <- kappamat[,jj]
      gamma[1:length(z),jj] <- kappamat[,jj] * phi[jj]
    }else{
      if(jj < last_reproductive_age){
        z <- convolve(kappamat[,jj],rev(z), type = "o")
        gamma[1:length(z),jj] <- z * phi[jj]
      }else{
        z <- convolve(kappamat[,jj],rev(z), type = "o")
        gamma[1:length(z),jj] <- z * (phi[jj] + 1 - sum(phi[1:jj]))
      }
    }
  }
  return(rowSums(gamma)) # Get the LRS distribution
}
