##################################################
# Thanks to Wenyun Zuo for the LRS functions
##################################################

# Calculate the lifetime reproductive success distribution for age only models.

s_to_pvec_death <- function(matrices, age.at.death){
  # as many age classes as age.at.death
  pvec <- vector()
  for(i in 1:(age.at.death))
  {
    pvec1 <- diag(matrices[[i]]$S)
    pvec <- cbind(pvec, pvec1)
  }
  vec0 <- rep(0, ncol = length(pvec1))
  pvec <- cbind(pvec, vec0)
  return(pvec)
}

r_to_kappamat1_death <- function(matrices, age.at.death){
  # as many age classes as age.at.death
  kappamat <- list()
  kappamat1 <- list()
  # i=1
  # i=2
  #
  for(i in 1:(age.at.death))
  {
    kappamat1 <- list(matrices[[i]]$R_kappa)
    kappamat <- c(kappamat, kappamat1)
  }
  return(kappamat)
}

g_to_gmat_death <- function(matrices, age.at.death){
  G <- list()
  G1 <- list()
  for(i in 1:(age.at.death))
  {
    G1 <- list(matrices[[i]]$G)
    G <- c(G, G1)
  }
  return(G)
}


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

# Using fast Fourier transform (FFT) methods to obtain the lifetime reproductive success distribution
# for stage only cases.


lrs_fft_stage_only <- function(Fmatrix, Umatrix,
                               max.offspring,
                               initial.stage = 1,
                               end.stage = 1,
                               distri_type){# = "poisson"){
  ## For reproduction, there are two type of distribution considered, poisson and bernoulli
  ####===== FFT method to calculate the LRS distribution
  U <- Umatrix
  if(distri_type == "poisson"){
    kappamat <- fpois_to_kappamat(Fmatrix, max.offspring - 1)
  }else if(distri_type == "bernoulli"){
    kappamat <- fbern_to_kappamat(Fmatrix, max.offspring - 1)
  }
  kappamat[is.na(kappamat)] <- 0
  n.kids <- dim(kappamat)[1]
  n.stages <- dim(kappamat)[2]
  # p <- fftw::planFFT(n.kids)
  p <- max.offspring
  alphahat.mat <- betahat.mat <- matrix(0, nrow = n.kids, ncol = n.stages)
  for (i in 1:n.stages) {
    alphahat.mat[,i] <- fft(kappamat[,i])
  }
  I <- diag(1, n.stages)
  e <- rep(1, n.stages)
  for (i in 1:n.kids) {
    W <- diag(alphahat.mat[i,])
    betahat.mat[i,] <- W %*% solve(I - t(U) %*% W) %*% (I - t(U)) %*% e
  }
  # gamma <- Re(fftw::IFFT(betahat.mat[, initial.stage], plan = p))
  #
  gamma <- Re(fft(betahat.mat[, initial.stage], inverse = T))/p
  if(end.stage < initial.stage){
    end.stage <- initial.stage
  }
  if(end.stage > initial.stage){
    for (i in (initial.stage + 1):end.stage) {
      gamma <- cbind(gamma, Re(fft(betahat.mat[,i],inverse = T))/p)
    }
    colnames(gamma) <- paste("stage", initial.stage:end.stage, sep = "")
  }
  return(gamma)
}


# Use the hybrid method to calculate the lifetime reproductive success distribution for age + stage cases.

age_stage_hybrid <- function(population.name, kappamat, G, pvec,
                             initial.stage = 1,
                             max.offspring){
  omega <- length(kappamat)
  n.stages <- dim(pvec)[1]
  ##============== Age + stage methods: FFT method
  #### age methods for first omega - 1 years
  gamma_part1 <- dynamic_programming(kappamat[-omega], G[-omega], pvec[, -omega], max.offspring = max.offspring, initial.stage)
  #### stage methods for years from omega, max.offspring need to be set to different species, 2^4 is for Roe deer
  Fmat <- matrix(rep(kappamat[[omega]][2,], n.stages), nrow = n.stages, ncol = n.stages, byrow = T)
  Umat <- G[[omega]] %*% diag(pvec[, omega])
  Mhat <- lrs_fft_stage_only(Fmat, Umat, max.offspring,
                             initial.stage = initial.stage, end.stage = n.stages,
                             distri_type = "bernoulli") # distri_tuype = c("bernoulli", "poisson")
  gamma_part2_semi <- matrix(0, nrow = (length(Mhat[, 1])+ length(gamma_part1[, 1])-1), ncol = n.stages)
  for (i in seq_len(n.stages)) {
    gamma_part2_semi[, i] <- convolve(Mhat[, i], rev(gamma_part1[, i]), type = "o")
  }
  gamma_part2 <- rowSums(gamma_part2_semi)
  #### add two part together
  final.gamma <- c(gamma_part1[, n.stages + 1], rep(0, length(gamma_part2) - length(gamma_part1[, n.stages + 1]))) +
    gamma_part2
  gamma.distr <- data.frame(density = final.gamma,
                            mids = seq_along(final.gamma) - 1)
  ##====plot the LRS distribution and the theoretic distribution
  plot_lrs_distr_allpop(population.name, gamma.distr)
  mtext(population.name)
  # ##====Gini index
  gini_plot(population.name, gamma.distr)
  ##====Return the LRS distribution
  return(gamma.distr)
}


# Create big matrix to use in fast Fourier transform (FFT) method for age+stage cases.

bigmat_agestage <- function(kappamat, G, pvec, max.offspring){
  omega <- length(kappamat) # kappamat is a list here. Each matrix in the list is corresponsing to an age.
  n.stages <- dim(G[[1]])[1]

  ### convert list of kappamat, G to single big matrices
  kappamat.single <- c(100,100)
  U <- matrix(0, nrow = omega*n.stages, ncol = omega*n.stages)
  for(i in seq_len(omega)){
    if(i < omega){
      U[i*n.stages + (1:n.stages), (i-1)*n.stages + (1:n.stages)] <- G[[i]] %*% diag(pvec[,i])
    }else{
      U[(i-1)*n.stages + (1:n.stages), (i-1)*n.stages + (1:n.stages)] <- G[[i]] %*% diag(pvec[,i])
    }
    kappamat.single <- cbind(kappamat.single, kappamat[[i]])
  }
  kappamat.fft <- matrix(0, nrow = max.offspring, ncol = omega*n.stages)
  #kappamat.fft[1:2,] <- kappamat.single[,-1]
  kappamat.fft[seq_len(dim(kappamat[[1]])[1]),] <- kappamat.single[,-1]
  return(list(kappamat = kappamat.fft,
              Umatrix = U))
}



# An age+stage model has a unique age+stage combination is written (a,s),
# and there are A x S such combinations.
# In some cases, the general method of the preceding section may be computationally lengthy
# and the block method described below is faster.


block_matrices_solve <- function(Klist, Glist, Pmatrix,
                                 max.offspring = 2^4,
                                 initial.age = 1,
                                 end.age = 1,
                                 stage_per_age){
  n.ages <- length(Klist)
  n.stages <- dim(Klist[[1]])[2]
  p <- max.offspring#fftw::planFFT(max.offspring)
  pad.mat <- matrix(0, nrow = max.offspring - dim(Klist[[1]])[1], ncol = n.stages)
  Klist.pad <- Klist # set kappahat.pad as a list as Klist
  for(i in seq_len(n.ages)){
    Klist.pad[[i]] <- rbind(Klist[[i]], pad.mat) # pad the fertility matrix to desired size
  }
  kappahat.mat <- Klist.pad # set kappahat.mat as a list as Klist.pad
  for(i in seq_len(n.ages)){
    Klist.pad[[i]][is.na(Klist.pad[[i]])] <- 0
    for (j in seq_len(n.stages)) {
      #kappahat.mat[[i]][,j] <- fftw::FFT(Klist.pad[[i]][,j], plan=p)
      kappahat.mat[[i]][,j] <- fft(Klist.pad[[i]][,j])
    }
  }
  I <- diag(n.stages)
  B <- list()
  d <- 1 - Pmatrix # colunm is age, row is stage
  d.hat <- d # set up d.hat and w.mat as a matrix with same dimension
  betahat.mat <-matrix(0, nrow = max.offspring, ncol = n.stages*n.ages)
  for (k in 1:max.offspring) {
    wvec <- c()
    dhatvec <- c()
    for (i in seq_len(n.ages - 1)) {
      wvec <- c(wvec, kappahat.mat[[i]][k,])
      W <- diag(kappahat.mat[[i + 1]][k,])
      Q <- t(Glist[[i]] %*% diag(Pmatrix[,i]))
      B[[i]] <- -Q %*% W # B[[i]] the last matrix in the list is not useful. The last stage will use C
    }
    wvec <- c(wvec, kappahat.mat[[n.ages]][k,])
    Q <- t(Glist[[n.ages]] %*% diag(Pmatrix[,n.ages]))
    C <- I - Q %*% W
    d.hat[, n.ages] <- solve(C) %*% d[, n.ages]
    dhatvec <- d.hat[, n.ages]

    for (i in seq_len(n.ages - 1)) { # for n.ages - 1 colunms
      d.hat[, n.ages - i] <- d[, n.ages - i] - B[[n.ages - i]] %*% d.hat[, n.ages - i + 1]
      dhatvec <- c(d.hat[, n.ages - i], dhatvec)
    }
    #betahat.mat[k,] <- diag(wvec) %*% dhatvec # element products, as.vector(d.hat) read matrix from colunm to colunm
    betahat.mat[k,] <- wvec * dhatvec # element products, as.vector(d.hat) read matrix from colunm to colunm
  }

  #gamma <- Re(fftw::IFFT(betahat.mat[, (initial.age - 1) * stage_per_age + 1], plan = p))
  gamma <- Re(fft(betahat.mat[, (initial.age - 1) * stage_per_age + 1], inverse = T))/p
  if(end.age < initial.age){
    end.age <- initial.age
  }
  for (i in 2:(stage_per_age * end.age)) {
    gamma <- cbind(gamma, Re(fft(betahat.mat[,i], inverse = T))/p)
  }
  #write.table(gamma, paste("age", (end.age- initial.age + 1), "X", stage_per_age, "stages.roedeer.block_fft", max.offspring,".txt", sep = ""), quote = F, row.names = F, col.names = F)
  return(gamma)
}

# Calculate mean of an arbitrary (discrete) probability distribution

distr_mean <- function(distr){
  ##distr is a data frame which has density and mids
  output <- sum(distr$density * distr$mids)
  return(output)
}


# Calculate variance of an arbitrary (discrete) probability distribution.

distr_var <- function(distr){
  ##distr is a data frame which has density and mids
  output <- sum(distr$density*(distr$mids - distr_mean(distr))^2)
  return(output)
}


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


g_to_gmat <- function(matrices){
  G <- list()
  G1 <- list()
  for(i in 1:(age.at.death))
  {
    G1 <- list(matrices[[i]]$G)
    G <- c(G, G1)
  }
  return(G)
}

# Obtain/set the last reproductive age.

last_rep_age <- function(kappamat){
  last_reproductive_age <- max(which(kappamat[1,] != 1))
  return(last_reproductive_age)
}


# Calculate the lifetime reproductive success distribution for age + stage cases
# by treating each age stage combination as a stage, and then use the stage only method.


lrs_fft_agestage <- function(kappamat, U, stage_per_age,
                             initial.age = 1,
                             end.age = 1){
  ####===== FFT method to calculate the LRS distribution
  kappamat[is.na(kappamat)] <- 0
  max.offspring <- dim(kappamat)[1]
  n.stages <- dim(kappamat)[2]
  #set.plan <- fftw::planFFT(max.offspring)
  set.plan <- max.offspring
  alphahat.mat <- betahat.mat <-matrix(0, nrow = max.offspring, ncol = n.stages)
  for (i in 1:n.stages) {
    alphahat.mat[,i] <- fft(kappamat[,i])
  }
  I <- diag(1, n.stages)
  e <- rep(1, n.stages)
  d <- (I - t(U)) %*% e
  for (i in 1:max.offspring) {
    W <- diag(alphahat.mat[i,])
    B1 <- solve(I - t(U) %*% W) # take about 39 seconds
    betahat.mat[i,] <- W %*% B1 %*% d
  }
  gamma <- Re(fft(betahat.mat[, (initial.age - 1) * stage_per_age + 1], inverse = T))/set.plan
  if(end.age < initial.age){
    end.age <- initial.age
  }
  for (i in 2:(stage_per_age * end.age)) {
    gamma <- cbind(gamma, Re(fft(betahat.mat[,i], inverse = T))/set.plan)
  }

  write.table(gamma, paste((end.age - initial.age + 1), "agesX", stage_per_age, "stages.roedeer.fft", max.offspring,".txt", sep = ""), quote = F, row.names = F, col.names = F)
  return(gamma) # gamma is a matrix including all initial stages
}

lrs_fft_agestage <- function(kappamat, U, stage_per_age,
                             initial.age = 1,
                             end.age = 1){
  ####===== FFT method to calculate the LRS distribution
  kappamat[is.na(kappamat)] <- 0
  max.offspring <- dim(kappamat)[1]
  n.stages <- dim(kappamat)[2]
  #set.plan <- fftw::planFFT(max.offspring)
  set.plan <- max.offspring
  alphahat.mat <- betahat.mat <-matrix(0, nrow = max.offspring, ncol = n.stages)
  for (i in 1:n.stages) {
    alphahat.mat[,i] <- fft(kappamat[,i])
  }
  I <- diag(1, n.stages)
  e <- rep(1, n.stages)
  d <- (I - t(U)) %*% e
  for (i in 1:max.offspring) {
    W <- diag(alphahat.mat[i,])
    B1 <- solve(I - t(U) %*% W) # take about 39 seconds
    betahat.mat[i,] <- W %*% B1 %*% d
  }
  gamma <- Re(fft(betahat.mat[, (initial.age - 1) * stage_per_age + 1], inverse = T))/set.plan
  if(end.age < initial.age){
    end.age <- initial.age
  }
  for (i in 2:(stage_per_age * end.age)) {
    gamma <- cbind(gamma, Re(fft(betahat.mat[,i], inverse = T))/set.plan)
  }

  write.table(gamma, paste((end.age - initial.age + 1), "agesX", stage_per_age, "stages.roedeer.fft", max.offspring,".txt", sep = ""), quote = F, row.names = F, col.names = F)
  return(gamma) # gamma is a matrix including all initial stages
}

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

r_to_kappamat1 <- function(matrices){
  # as many age classes as age.at.death
  kappamat <- list()
  kappamat1 <- list()
  # i=1
  # i=2
  #
  for(i in 1:(age.at.death))
  {


    # matF[1,] <- diag(matrices[[i]]$R)
    # add.zero <- matrix(0, ncol = ncol(M$R), nrow=ncol(M$R)-1)
    # matF[2:ncol(M$R),] <- add.zero
    # matrices[[1]]$
    kappamat1 <- list(matrices[[i]]$R_kappa)
    # kappamat1 <- list(matF)
    kappamat <- c(kappamat, kappamat1)
  }
  return(kappamat)
}


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

# convert to lists for age stage
s_to_pvec <- function(matrices){
  # as many age classes as age.at.death
  pvec <- vector()
  for(i in 1:(age.at.death))
  {
    pvec1 <- diag(matrices[[i]]$S)
    pvec <- cbind(pvec, pvec1)
  }
  vec0 <- rep(0, ncol = length(pvec1))
  pvec <- cbind(pvec, vec0)
  return(pvec)
}

