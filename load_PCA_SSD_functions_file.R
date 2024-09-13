# Harman Jaggi
# Functions for running PCA and SSD results in the manuscript
# Code for Ecology Letters manuscript on Sep 2024: Density dependence

# Set the body size and life history parameters for the species. Soay sheep in our case.
minsize <- 0.1  ### make this smaller than observed value
maxsize <- 38   ### make this larger than observed value
nn <- n <- 50
max.age <- 16

## Define the IPM functions

# Survival function S(z,t)
# Sets up survival function based on formula in Coulson 2012

S.fun <- function(z,intercept,z.slope,n.slope,year.eff,id.eff,N) {
  u<-exp(intercept+z.slope*z+n.slope*N+year.eff+id.eff)
  return((u/(1+u)))
}

# Recruitment function. Note this is identical in this case to the survival function.  Needs altering for a different mating system
R.fun <- function(z,intercept,z.slope,n.slope,year.eff,id.eff,N) {
  u <- exp(intercept+z.slope*z+n.slope*N+year.eff+id.eff)
  return((u/(1+u)))
}

G.fun <- function(z, zz, intercept.mu,
                  z.slope.mu, n.slope.mu, year.eff.mu,
                  id.eff.mu, intercept.va, z.slope.va,
                  n.slope.va, year.eff.va, id.eff.va, N) {
  
  mu.z <- intercept.mu + z.slope.mu*z + n.slope.mu*N + year.eff.mu + id.eff.mu # Eg(z)
  sigma.z2 <- intercept.va+z.slope.va*z+n.slope.va*N+year.eff.va+id.eff.va
  sigma.z2 <- ifelse(sigma.z2<0,0.0001,sigma.z2)
  sigma.z <- sqrt(sigma.z2)
  temp1 <- sqrt(2*pi)*sigma.z
  temp2 <- ((zz-mu.z)^2)/(2*sigma.z2)
  return(exp(-temp2)/temp1)
}

# expected growth increment
G.inc.fun <- function(z, intercept.mu,
                      z.slope.mu, n.slope.mu, year.eff.mu,
                      id.eff.mu, intercept.va, z.slope.va,
                      n.slope.va, year.eff.va, id.eff.va, N) {
  
  junk <- data.frame()
  mu.z <- c()
  g.inc <- c()
  diffg <- c()
  for (i in 1:length(z)) {
    mu.z[i] <- intercept.mu + z.slope.mu*z[i] + n.slope.mu*N + year.eff.mu + id.eff.mu # Eg(z)
    diffg[i] <- mu.z[i] - z[i]
    
  }
  return(diffg)
}

# expected growth increment ratio
G.prop.ratio.fun <- function(z, intercept.mu,
                             z.slope.mu, n.slope.mu, year.eff.mu,
                             id.eff.mu, intercept.va, z.slope.va,
                             n.slope.va, year.eff.va, id.eff.va, N) {
  
  junk <- data.frame()
  mu.z <- c()
  g.inc <- c()
  ratiog <- c()
  for (i in 1:length(z)) {
    mu.z[i] <- intercept.mu + z.slope.mu*z[i] + n.slope.mu*N + year.eff.mu + id.eff.mu # Eg(z)
    ratiog[i] <- (mu.z[i]-z[i])/(z[i])
  }
  return(list(mu.z, z, ratiog))
}

# Inheritance function D(z'|z)  zz = z' for R code.  Note same structure as G in this case -- could have different probability distribution if needed
D.fun <- function(z,zz, intercept.mu, z.slope.mu,
                  n.slope.mu,year.eff.mu,id.eff.mu,intercept.va,
                  z.slope.va,n.slope.va,year.eff.va,id.eff.va,N) {
  mu.z <- intercept.mu+z.slope.mu*z+n.slope.mu*N+year.eff.mu+id.eff.mu
  sigma.z2 <- intercept.va+z.slope.va*z+n.slope.va*N+year.eff.va+id.eff.va
  sigma.z2 <- ifelse(sigma.z2<0,0.0001,sigma.z2)
  sigma.z <- sqrt(sigma.z2)
  temp1 <- sqrt(2*pi)*sigma.z
  temp2 <- ((zz-mu.z)^2)/(2*sigma.z2)
  return(exp(-temp2)/temp1)
}

# The function estimates matrix outputs by running the IPM functions
bigmatrix_agestage_eqm <-function(n, s.params,r.params,g.params,d.params,N) {
  
  # boundary points b and mesh points y
  b <- minsize+c(0:n)*(maxsize-minsize)/n
  y <- 0.5*(b[1:n]+b[2:(n+1)])
  # y <- y*sad
  
  # create S, R, G and D matrices and f0 and f1
  
  S <- (diag(S.fun(y,s.params[1],s.params[2],s.params[3],s.params[4], s.params[5],N)))
  
  # Reproduction distribution or kappamat for binomial case cuz you have 0 or 1 kids
  R_kappa <- matrix(c(1-R.fun(y, r.params[1],r.params[2],r.params[3],r.params[4], r.params[5],N),
                      R.fun(y,r.params[1],r.params[2],r.params[3],r.params[4], r.params[5],N)),
                    nrow=2, byrow=T)
  
  R <- (diag(R.fun(y,r.params[1],r.params[2],r.params[3],r.params[4], r.params[5],N)))
  G <- (t(outer(y,y,G.fun,g.params[1],g.params[2],g.params[3], g.params[4],g.params[5],g.params[6],g.params[7],g.params[8],g.params[9], g.params[10],N)))
  D <- (t(outer(y,y,D.fun,d.params[1],d.params[2],d.params[3], d.params[4],d.params[5],d.params[6],d.params[7],d.params[8],d.params[9], d.params[10],N)))
  
  svec <- S.fun(y,s.params[1],s.params[2],s.params[3],s.params[4], s.params[5],N)
  rvec <- R.fun(y,r.params[1],r.params[2],r.params[3],r.params[4], r.params[5],N)
  g.inc.vec <- G.inc.fun(y, g.params[1], g.params[2], g.params[3], g.params[4], g.params[5], g.params[6], g.params[7], g.params[8], g.params[9], g.params[10], N)
  # g.ratio.vec <- G.ratio.fun(y, g.params[1], g.params[2], g.params[3], g.params[4], g.params[5], g.params[6], g.params[7], g.params[8], g.params[9], g.params[10], N)
  # g.prop.ratio.vec <- G.prop.ratio.fun(y, g.params[1], g.params[2], g.params[3], g.params[4], g.params[5], g.params[6], g.params[7], g.params[8], g.params[9], g.params[10], N)
  g.prop.ratio.check1 <- G.prop.ratio.fun(y, g.params[1], g.params[2], g.params[3], g.params[4], g.params[5], g.params[6], g.params[7], g.params[8], g.params[9], g.params[10], N)[[1]]
  g.prop.ratio.check2 <- G.prop.ratio.fun(y, g.params[1], g.params[2], g.params[3], g.params[4], g.params[5], g.params[6], g.params[7], g.params[8], g.params[9], g.params[10], N)[[2]]
  g.prop.ratio.check3 <- G.prop.ratio.fun(y, g.params[1], g.params[2], g.params[3], g.params[4], g.params[5], g.params[6], g.params[7], g.params[8], g.params[9], g.params[10], N)[[3]]
  # d.ratio.vec <- D.ratio.fun(y, d.params[1], d.params[2], d.params[3], d.params[4], d.params[5], d.params[6], d.params[7], d.params[8], d.params[9], d.params[10], N)
  
  
  
  # scale D and G so columns sum to 1
  G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  D <- D/matrix(as.vector(apply(D,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  return(list(S=S, R=R, R_kappa=R_kappa, G=G, D=D, svec=svec, rvec=rvec,
              g.inc=g.inc.vec,
              # d.ratio=d.ratio.vec,
              g.prop.ratio.check1=g.prop.ratio.check1, #Ez
              g.prop.ratio.check2=g.prop.ratio.check2, #z
              g.prop.ratio.check3=g.prop.ratio.check3, #ratio
              meshpts=y))
}

# The function estimates matrix outputs for LRS analysis by running the IPM functions
bigmatrix_agestage_lrs <-function(n, s.params,r.params,g.params,d.params,N) {
  
  # boundary points b and mesh points y
  b <- minsize+c(0:n)*(maxsize-minsize)/n
  y <- 0.5*(b[1:n]+b[2:(n+1)])
  # y <- y*sad
  # S.fun(y,s.params[1],s.params[2],s.params[3],s.params[4], s.params[5],N)
  
  S <- (diag(S.fun(y,s.params[1],s.params[2],s.params[3],s.params[4], s.params[5],N)))
  
  ## Reproduction distribution or kappamat for binomial case cuz you have 0 or 1 kids
  R_kappa <- matrix(c(1-R.fun(y, r.params[1],r.params[2],r.params[3],r.params[4], r.params[5],N),
                      R.fun(y,r.params[1],r.params[2],r.params[3],r.params[4], r.params[5],N)),
                    nrow=2, byrow=T)
  R <- (diag(R.fun(y,r.params[1],r.params[2],r.params[3],r.params[4], r.params[5],N)))
  G <- (t(outer(y,y,G.fun,g.params[1],g.params[2],g.params[3], g.params[4],g.params[5],g.params[6],g.params[7],g.params[8],g.params[9], g.params[10],N)))
  D <- (t(outer(y,y,D.fun,d.params[1],d.params[2],d.params[3], d.params[4],d.params[5],d.params[6],d.params[7],d.params[8],d.params[9], d.params[10],N)))
  
  svec <- S.fun(y,s.params[1],s.params[2],s.params[3],s.params[4], s.params[5],N)
  rvec <- R.fun(y,r.params[1],r.params[2],r.params[3],r.params[4], r.params[5],N)
  gvec <- G.fun(y,y, g.params[1],g.params[2],g.params[3], g.params[4],g.params[5],g.params[6],g.params[7],g.params[8],g.params[9], g.params[10],N)
  dvec <- D.fun(y,y,d.params[1],d.params[2],d.params[3], d.params[4],d.params[5],d.params[6],d.params[7],d.params[8],d.params[9], d.params[10],N)
  g.inc.vec <- G.inc.fun(y, g.params[1], g.params[2], g.params[3], g.params[4], g.params[5], g.params[6], g.params[7], g.params[8], g.params[9], g.params[10], N)
  # g.ratio.vec <- G.ratio.fun(y, g.params[1], g.params[2], g.params[3], g.params[4], g.params[5], g.params[6], g.params[7], g.params[8], g.params[9], g.params[10], N)
  g.prop.ratio.vec <- G.prop.ratio.fun(y, g.params[1], g.params[2], g.params[3], g.params[4], g.params[5], g.params[6], g.params[7], g.params[8], g.params[9], g.params[10], N)
  
  # scale D and G so columns sum to 1
  G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  D <- D/matrix(as.vector(apply(D,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  return(list(S=S, R=R, R_kappa=R_kappa, G=G, D=D, svec=svec, rvec=rvec, g.inc=g.inc.vec,
              # g.ratio=g.ratio.vec,
              g.prop.ratio=g.prop.ratio.vec, meshpts=y))
}

# The function estimates equilibrium carrying capacity and returns K
# for each parameter set
eqm_func  <- function(toprms){
  
  s.params <- c(toprms[1], toprms[2],toprms[3],0,0)
  r.params <- c(toprms[4],toprms[5],toprms[6],0,0)
  g.params <- c(toprms[7],toprms[8],toprms[9],0,0,toprms[10],toprms[11],0,0,0)
  d.params <- c(toprms[12],toprms[13],toprms[14],0,0,toprms[15],toprms[16],0,0,0)
  
  # check the distribution of K for different starting population
  
  N <- 80
  diff <- 1
  diff_new <- 1
  junk1 <- data.frame()
  tol <- 0.0000001
  
  while(diff>0){
    # N <- sum(nt)
    M <- bigmatrix_agestage_lrs(nn,s.params,r.params,g.params,d.params,N)
    
    # matF <- matrix(NA, ncol=ncol(M$R), nrow=ncol(M$R))
    
    # mat <- M$G %*% M$S + M$D %*% matF
    mat <- M$G %*% M$S + M$D %*% M$R
    mat <- as.matrix(mat)
    
    # # using block matrix to compare R0 etc
    # matF <- M$D %*% M$R
    # matU <- M$G %*% M$S
    #
    # block_matrix_func(matF, matU, age.at.death = 16)
    
    lam <- Re(eigen(mat)$values[1])
    # matU <-  M$G %*% M$S
    # rec.mat <- M$D %*% M$R
    # print(lambda)
    diff <- (lam - 1) #abs(sum(nt.new)-sum(nt))
    junk <- data.frame(diff=diff, lam=lam, old=N, new=N+10)
    print(junk)
    N <- N+10
    # junk <- data.frame(old = N, new = (N+1), diff=diff, lambda=lambda)
    # junk1 <- rbind(junk1, junk)
    # print(junk1)
    # df <- rbind(df, junk1)
  }
  
  N = N-25
  diff_abs <- 1
  
  while(diff_abs > tol)
  {
    if(diff_new>0)
    {
      
      M <- bigmatrix_agestage_lrs(nn,s.params,r.params,g.params,d.params,N)
      mat <- M$G %*% M$S + M$D %*% M$R
      mat <- as.matrix(mat)
      
      lam <- Re(eigen(mat)$values[1])
      
      diff_abs <- abs(lam - 1) #abs(sum(nt.new)-sum(nt))
      diff_new <- (lam - 1) #abs(sum(nt.new)-sum(nt))
      junk_new <- data.frame(diff_new=diff_new, lam=lam, diff_abs=diff_abs, old=N, new=N+0.5)
      print(junk_new)
      N <- N+1
    }
    else
    {print("The diff is negative now. Stop")
      break
    }
    
  }
  
  N <- N
  return(N)
  
}

# For large matrices it is faster to calculate the dominant eigenvalue
# and eigenvectors associated with it via iteration rather than using eigen.
get.eigen.stuff <- function(mat){
  sz <- dim(mat)[1]
  t.now <- runif(sz)
  t.now <- t.now/sum(t.now)
  t.next <- mat%*%t.now
  t.next <- t.next/sum(t.next)
  i <- 0
  while (sum(abs(t.next-t.now))>0.0000001){
    i <- i+1
    print(i)
    t.now <- t.next
    t.next <- mat%*%t.now
    lambda <- sum(t.next)/sum(t.now)
    t.next <- t.next/sum(t.next)
  }
  r.now <- runif(sz)
  r.now <- r.now/sum(r.now)
  r.next <- r.now%*%mat
  r.next <- r.next/sum(r.next)
  while (sum(abs(r.next-r.now))>0.0000001){
    r.now <- r.next
    r.next <- r.now%*%mat
    r.next <- r.next/sum(r.next)
  }
  return(list(lambda,t.next,r.next))
}

# The function estimates generation time for species
gentime <- function(M){
  Pmat <- M$G%*%M$S
  Fmat <- M$D%*%M$R
  L <- A <- array(0,c(nn,nn,max.age))
  L[,,1] <- diag(nn)
  A[,,1] <- Fmat %*% L[,,1]
  for (i in 2:max.age){
    L[,,i] <- Pmat%*%L[,,i-1]
    A[,,i] <- Fmat%*%L[,,i]
  }
  Amat <- apply(A, c(1,2), sum)
  tt <- get.eigen.stuff(Amat)
  c <- tt[[2]]
  d <- tt[[3]]
  temp <- as.vector(d%*%c)
  d <- d/temp
  R0 <- tt[[1]]
  phi <- rep(0,max.age)
  for (i in 1:max.age){
    phi[i] <- (d%*%Fmat%*%L[,,i]%*%c)/R0
  }
  Tc <- sum(0:(max.age-1)*phi)
  return(Tc)
}

# Function to make block matrix for age and stage
block_matrix_func <- function(Fmat, Umat, age.at.death){
  
  n.stages <- nrow(Fmat)
  F_row <- matrix(c(rep(Fmat, age.at.death)), nrow=n.stages, ncol=n.stages*age.at.death)
  F_block_zeroes <- matrix(0, nrow = n.stages*age.at.death-n.stages, ncol=age.at.death*n.stages)
  F_block_mat <- as.matrix(rbind(F_row, F_block_zeroes))
  # dim(F_row)
  
  
  ### convert to single big matrices
  U_block <- matrix(0, nrow = age.at.death*n.stages-n.stages, ncol = age.at.death*n.stages)
  zero_mat <- matrix(0, ncol = n.stages, nrow=n.stages)
  # dim(U_block)
  # i=1
  for(i in seq_len(age.at.death-1)){
    # U_block[(i-1)*n.stages + (1:n.stages), (i-1)*n.stages + (1:n.stages)] <- Umat
    U_block[(i-1)*n.stages + (1:n.stages), (i-1)*n.stages + (1:n.stages)] <- Umat
    U_block[(age.at.death-2)*n.stages + (1:n.stages) , (age.at.death-2)*n.stages + (1:n.stages)] <- zero_mat
  }
  # dim(U_block)
  U_block_zeroes <- matrix(0, nrow = n.stages, ncol=age.at.death*n.stages)
  U_block_mat <- as.matrix(rbind(U_block_zeroes, U_block))
  block_mat <- as.matrix(rbind(F_row, U_block))
  
  # block_mat1 <- F_block_mat+U_block_mat
  # block_mat==block_mat1
  # dim(U_block_mat)
  # dim(block_mat)
  
  return(list(block_mat, U_block_mat, F_block_mat))
}

# For the block matrix of age and stage, the function estimates equilibrium carrying capacity and returns K
# for each parameter set
eqm_dist_block_matrix_func  <- function(toprms){
  
  s.params <- c(toprms[1], toprms[2],toprms[3],0,0)
  r.params <- c(toprms[4],toprms[5],toprms[6],0,0)
  g.params <- c(toprms[7],toprms[8],toprms[9],0,0,toprms[10],toprms[11],0,0,0)
  d.params <- c(toprms[12],toprms[13],toprms[14],0,0,toprms[15],toprms[16],0,0,0)
  
  # check the distribution of K for different starting population
  
  N <- 80
  diff <- 1
  diff_new <- 1
  junk1 <- data.frame()
  tol <- 0.0000001
  
  while(diff>0){
    # N <- sum(nt)
    M <- bigmatrix_agestage_lrs(nn,s.params,r.params,g.params,d.params,N)
    
    # matF <- matrix(NA, ncol=ncol(M$R), nrow=ncol(M$R))
    
    # mat <- M$G %*% M$S + M$D %*% matF
    mat <- M$G %*% M$S + M$D %*% M$R
    mat <- as.matrix(mat)
    
    # # using block matrix to compare R0 etc
    # matF <- M$D %*% M$R
    # matU <- M$G %*% M$S
    #
    
    # this is overall matrix
    matU <- M$G %*% M$S
    mat <- M$G %*% M$S +  M$D %*% M$R
    
    
    # using block matrix to compare R0 etc
    Fmat <- M$D %*% M$R
    Umat <- M$G %*% M$S
    
    block_mat <- block_matrix_func(Fmat, Umat, age.at.death = 16)
    # Re(eigen(block_mat)$values[1])
    (block_mat)
    # faster eigen calculation!
    lam <- lambda_block <- Re(RSpectra::eigs(block_mat[[1]], 1)$values[1])
    diff <- (lam - 1) #abs(sum(nt.new)-sum(nt))
    junk <- data.frame(diff=diff, lam=lam, old=N, new=N+10)
    print(junk)
    N <- N+10
  }
  
  N = N-25
  diff_abs <- 1
  
  while(diff_abs > tol)
  {
    if(diff_new>0)
    {M <- bigmatrix_agestage_lrs(nn,s.params,r.params,g.params,d.params,N)
    mat <- M$G %*% M$S + M$D %*% M$R
    mat <- as.matrix(mat)
    
    
    # lam <- Re(eigen(mat)$values[1])
    lam <- lambda_block <- Re(RSpectra::eigs(block_mat[[1]], 1)$values[1])
    diff_abs <- abs(lam - 1) #abs(sum(nt.new)-sum(nt))
    diff_new <- (lam - 1) #abs(sum(nt.new)-sum(nt))
    junk_new <- data.frame(diff_new=diff_new, lam=lam, diff_abs=diff_abs, old=N, new=N+0.5)
    print(junk_new)
    N <- N+1
    }
    else
    {print("The diff is negative now. Stop")
      break
    }
  }
  
  N <- N
  return(N)
  # return(junk1)
}


# Function to add K to original data frame
addK_func <- function(dat)
{
  # ready the data for eqm size calculation
  temp_eqm <- as.data.frame(t(dat))
  
  # get the eqm numbers for each sample dat
  newK <- apply(temp_eqm, 2, eqm_func)
  
  # temp <- unlist(DD.nt.res, recursive = F)
  # N_eq1 <- sapply(temp, sum)
  newK.df <- as.data.frame(newK)
  
  # add the eqm column to the data (input for the main function)
  dateqm  <- dat %>% as.data.frame() %>%
    mutate(newK = newK.df$newK) %>%
    t() %>%
    as.data.frame()
  return(t(dateqm))
}

# Function to add K to original data frame
addK_block_mat_func <- function(dat)
{
  # ready the data for eqm size calculation
  temp_eqm <- as.data.frame(t(dat))
  
  # get the eqm numbers for each sample dat
  newK <- apply(temp_eqm, 2, eqm_dist_block_matrix_func)
  
  # temp <- unlist(DD.nt.res, recursive = F)
  # N_eq1 <- sapply(temp, sum)
  newK.df <- as.data.frame(newK)
  
  # add the eqm column to the data (input for the main function)
  dateqm  <- dat %>% as.data.frame() %>%
    mutate(newK = newK.df$newK) %>%
    t() %>%
    as.data.frame()
  return(t(dateqm))
}

# Function returns covariates for the PCA analysis
# the input arguments are the dataset (sample.datK in our case) and
# eqm_ratio which gives the ratio of K at which we want the results
# We do the analysis for 0, 1/2, and 1 which corresponds to population
# N=0, N=K/2, and N=K.

get_vitals_func <- function(dat, eqm_ratio){
  
  tbl <- data.frame()
  tblnew <- data.frame()
  tblnew.sad <- data.frame()
  
  dateqm <- data.frame(nprms=1:17)
  
  dat <- t(dat)
  sim_dat = ncol(dat)
  N_vec <- c()
  # i=11
  # eqm_ratio = 0
  for(i in 1:sim_dat)
  {
    dateqm$values <- dat[,i]
    dateqm <- as.data.frame(dateqm)
    toprms <- dateqm
    
    s.params <- c(toprms$values[1], toprms$values[2],toprms$values[3],0,0)
    r.params <- c(toprms$values[4],toprms$values[5],toprms$values[6],0,0)
    g.params <- c(toprms$values[7],toprms$values[8],toprms$values[9],0,0,toprms$values[10],toprms$values[11],0,0,0)
    d.params <- c(toprms$values[12],toprms$values[13],toprms$values[14],0,0,toprms$values[15],toprms$values[16],0,0,0)
    # N <- toprms$values[17]
    
    N <- toprms$values[17]*(eqm_ratio)
    N_vec[i] <- N
    
    print(c(i,N))
    
    M <- bigmatrix_agestage_eqm(nn, s.params, r.params, g.params, d.params, N)
    
    #with development
    mat <- M$G %*% M$S + M$D %*% M$R
    mat <- as.matrix(mat)
    
    # if(colSums(M$G %*% M$S)>1){
    #   print("matrix sum exceeds 1")
    # }
    
    lambda <- eigen(mat)$values[1]
    sad1 <- Re(eigen(mat)$vectors[,1]) # Stable size distribution
    sad <- sad1/sum(sad1)
    
    
    resvec <-  M
    n_stage <- length(resvec$g.inc)
    # both proportional ratio and difference become zero at same index so keep one
    g_index <- which(resvec$g.inc<0)[1]-1
    
    # normalize partitions
    sad_small1 <- sad[1:(g_index)]
    sad_small <- sad_small1/sum(sad_small1)
    
    sad_big1 <- sad[(g_index+1):(n_stage)]
    sad_big <- sad_big1/sum(sad_big1)
    
    
    # growth increment for small
    ginc.small.sad <- sum(resvec$g.inc[1:g_index]*sad_small)
    
    # growth increment for big ssd
    ginc.big.sad <- sum(resvec$g.inc[(g_index+1):n_stage]*sad_big)
    
    # growth proportional ratios
    gprat.small.sad <- sum(resvec$g.prop.ratio.check3[1:g_index]*sad_small)
    
    gprat.big.sad <- sum(resvec$g.prop.ratio.check3[(g_index+1):n_stage]*sad_big)
    
    drat.sad <- sum(resvec$d.ratio*sad)
    
    # surv small
    s.small.sad <- sum(resvec$svec[1:(g_index)]*sad_small)
    
    # surv big
    s.big.sad <- sum(resvec$svec[(g_index+1):n_stage]*sad_big)
    
    # repro
    r.sad <- sum(resvec$rvec*sad)
    
    tbl1 <- cbind(
      s.small.sad, s.big.sad,
      r.sad,
      ginc.small.sad, ginc.big.sad,
      gprat.small.sad, gprat.big.sad,
      drat.sad
    )
    
    tbl <- rbind(tbl1, tbl)
    
  }
  
  tbl <- as.data.frame(tbl)
  tbl$N_vec <- N_vec
  
  colnames(tbl) <- c(
    "Ss", "Sb",
    "R",
    "Gsi", "Gbi",
    "Gspr", "Gbpr",
    "Drat",
    "PopN")
  return(tbl)
}

# Function binds the results for all equilibrium ratios in one table.
get_named_tbl <- function(dat, eqratio)
{
  all_tbl <- data.frame()
  for (i in 1:length(eqm_val)){
    temp <- get_vitals_func(dat, eqm_val[i])
    temp <- as.data.frame(temp)
    temp$eqratio <- rep(eqratio[i], nrow(temp))
    all_tbl <- rbind(all_tbl, temp)
  }
  return(all_tbl)
}

# Function to check if a matrix if singular
sing_check <- function(m) class(try(solve(m),silent=T))=="matrix"

# Set the age of death.For soay sheep we set it to 16 years.
age.at.death <- 16

# Function to calculate LRS distributions at different densities
get_lrs_func <- function(dat, eqm_ratio){
  
  tbl <- data.frame()
  tblnew <- data.frame()
  tblnew.sad <- data.frame()
  # newborn_stage <- 1
  matrix_list <- list()
  dateqm <- data.frame(nprms=1:17)
  
  lrs_gamma.distr <- data.frame()
  
  dat <- t(dat)
  
  nsim <- ncol(dat)
  # eqm_ratio <- 0
  # sad_temp <- data.frame(sno=1:50)
  # i=2
  for(i in 1:nsim)
  {
    
    dateqm$values <- dat[,i]
    dateqm <- as.data.frame(dateqm)
    toprms <- dateqm
    # list.stuff <- get_matrix_agestage_N(datcov, N[j])[[1]] # any [[1]] except 16 is ok because the matrcies are constant
    
    s.params <- c(toprms$values[1], toprms$values[2],toprms$values[3],0,0)
    r.params <- c(toprms$values[4],toprms$values[5],toprms$values[6],0,0)
    g.params <- c(toprms$values[7],toprms$values[8],toprms$values[9],0,0,toprms$values[10],toprms$values[11],0,0,0)
    d.params <- c(toprms$values[12],toprms$values[13],toprms$values[14],0,0,toprms$values[15],toprms$values[16],0,0,0)
    # N <- toprms$values[17]
    
    N <- toprms$values[17]*(eqm_ratio)
    
    M <- bigmatrix_agestage_lrs(nn, s.params, r.params, g.params, d.params, N)
    
    ## we have 16 age classes such that survival is same for each age as it is in stage
    ## and the last age survival is 0.
    
    age.at.death <- 16
    age.stage.list <- list()
    # nn <- 50 # no. of stages
    surv.age.at.death <- rep(0, length(s.params))
    
    age.stage.list[[1]] <- bigmatrix_agestage_lrs(nn, s.params, r.params, g.params, d.params,N)
    
    for (j in 2:(age.at.death))
    {
      age.stage.list[[j]] <-     age.stage.list[[1]]
    }
    
    
    # age.stage.list[[age.at.death]] <- bigmatrix_agestage(nn, surv.age.at.death, r.params, g.params, d.params,N)
    age.stage.list[[age.at.death]]$S <- matrix(rep(0, nrow(age.stage.list[[1]]$S)*ncol(age.stage.list[[1]]$S)),
                                               nrow= nrow(age.stage.list[[1]]$S), byrow=T)
    age.stage.list[[1]]$S
    
    # create fertility and survival matrices
    matF <- matrix(NA, ncol=ncol(M$R), nrow=ncol(M$R))
    #
    # # create the recruitment matrix (here only stage 1 being born into)
    # matF[1,] <- M$R[2,]
    matF[1,] <- diag(M$R)
    add.zero <- matrix(0, ncol = ncol(M$R), nrow=ncol(M$R)-1)
    matF[2:ncol(M$R),] <- add.zero
    
    # this is overall matrix
    matU <- M$G %*% M$S
    mat <- M$G %*% M$S +  M$D %*% M$R
    
    # using block matrix to compare R0 etc
    Fmat <- M$D %*% M$R
    Umat <- M$G %*% M$S
    
    block_mat <- block_matrix_func(Fmat, Umat, age.at.death = 16)
    # Re(eigen(block_mat)$values[1])
    
    # faster eigen calculation!
    lambda_block <- Re(RSpectra::eigs(block_mat[[1]], 1)$values[1])
    # Re(eigen(block_mat)$values[1])
    
    mat <- as.matrix(mat)
    lambda <- Re(eigen(mat)$values[1])
    
    diff_lambda = lambda-lambda_block
    diff_lambda
    sad1 <- Mod(eigen(mat)$vector[,1])
    sad <- sad1/sum(sad1)
    # sum(sad)
    
    newF <-  M$D %*% M$R
    matrix_list[[i]] <- newF
    
    # newF <-  M$D %*% diag(M$R) *(sad)
    newborns_dist1 <- colSums(t(newF))
    
    # normalize the newborn distribution
    newborns_dist <- newborns_dist1/sum(newborns_dist1)
    # plot(newF)
    
    matrices <- age.stage.list
    n.stages <- dim(matrices[[1]]$R)[2] # obtain number of stages
    kappamat <- r_to_kappamat1(matrices) # list of n.stages matrix
    G <- g_to_gmat(matrices) # list of n.stages matrix
    pvec <- s_to_pvec(matrices) #size conditional survival, n.stages stages, row is corresponded to stage
    # dim(pvec)
    # length(G)
    # length(kappamat)
    max.offspring <- 30
    
    ptm <- proc.time()
    
    gamma <- block_matrices_solve(Klist = kappamat, Glist = G, Pmatrix = pvec,
                                  max.offspring = max.offspring,
                                  initial.age = 1,
                                  end.age = 1, # not 2 anymore mother born each stage then sum
                                  stage_per_age = n.stages)
    
    new_gamma <- gamma %*% newborns_dist
    # proc.time() - ptm
    print(paste0(ncol(dat), ",", i))
    
    #plot the new LRS distribution weighted by new born distribution
    
    # gamma.distr1
    temp <- data.frame(density = new_gamma,
                       mids = seq_along(gamma[,1]) - 1)
    
    # mean and variance from the lrs distribution
    junk1 <- temp %>% summarise(R0=sum(mids*density),
                                var=sum(mids*density - sum(mids*density))^2,
                                dd=sum(density))
    # junk1$R0
    
    # #### get R0 from transition matrix
    I <- diag(1, nrow = nrow(mat))
    N <- solve(I-matU)
    R0.mat <- M$D %*% M$R %*% N
    R0_eig <- Re(eigen(R0.mat)$values[1]) # which is same as R0.mat[1,1] # Stage based R0
    
    # block mat R0
    I_block <- diag(1, nrow=nn*age.at.death)
    U_block <- block_mat[[2]]
    F_block <- block_mat[[3]]
    N_block <- solve(I_block-U_block)
    R0_block_mat <- F_block %*% N_block
    # dim(R0_block)
    R0_block <- Re(eigen(R0_block_mat)$values[1])
    
    # Get Generation Time using function above
    Tc <- gentime(M)
    
    temp$sim <- rep(i, nrow(temp))
    temp$lambda <- rep(lambda, nrow(temp))
    temp$mean.lrs.dist <- rep(junk1$R0, nrow(temp))
    temp$Tc <- rep(Tc, nrow(temp))
    temp$R0_eig <- rep(R0_eig, nrow(temp))
    temp$R0_block <- rep(R0_block, nrow(temp))
    lrs_gamma.distr <- rbind(lrs_gamma.distr, temp)
  }
  
  return(lrs_gamma.distr)
}

# RESCALE AND GET EXPECTATION

# Function to calculate probability of having no offspring
expectation_lrsnot0 <- function(dat){
  mean_vals <- data.frame()
  nsim <- nrow(dat)
  # i=1
  for (i in 1:nsim){
    junk <- dat %>% filter(sim==i) %>%
      mutate(prob_kid = mids*density) %>%
      tail(-1) %>%
      mutate(rescaled_den=density/sum(density),
             lrsnot0=mids*rescaled_den,
             v_lrsnot0=(mids*rescaled_den - sum(mids*rescaled_den))^2) %>%
      summarise(mean_lrsnot0= sum(lrsnot0),
                mean_lrs = sum(prob_kid),
                var_lrsnot0 = sum(v_lrsnot0))
    mean_vals <- rbind(mean_vals, junk)
  }
  return(mean_vals)
}

# This function merges attributes of expected probability of having no offspring,
# expected mean LRS (RO),
# probability of having more than one offspring
meanlrs_vital_func <- function(dat){
  lrs_plotdat <- data.frame()
  nsim <- nrow(dat)
  # i =1
  for(i in 1:length(eqm_val))
  {
    temp <- junk_temp <- get_lrs_func(dat, eqm_ratio = eqm_val[i])
    
    # arrange
    arr_temp <- temp %>% arrange(mids)
    lrs0 <- arr_temp[1:nsim,]
    lrs_not0 <- (1-lrs0$density)
    
    
    exp_dat <- expectation_lrsnot0(temp)
    exp_dat <- filter(exp_dat, mean_lrs!=0)
    
    lrs_plotdat1 <- data.frame(plrs0=lrs0$density,
                               exp_dat,
                               plrs_not0 = lrs_not0,
                               Tc = lrs0$Tc,
                               R0_block = lrs0$R0_block,
                               # life_exp_birth = lrs0$life_exp_birth,
                               R0_eig=lrs0$R0_eig,
                               R0_dist=lrs0$mean.lrs.dist,
                               lambda = lrs0$lambda
    )
    
    lrs_plotdat1$den <- eqm_val[i]
    lrs_plotdat <- rbind(lrs_plotdat, lrs_plotdat1)
    # cbind(lrs_plotdat, sad_df)
  }
  return(lrs_plotdat)
}

# Function for mother's LRS distribution
matF_func <- function(dat, eqm_ratio){
  
  tbl <- data.frame()
  tblnew <- data.frame()
  tblnew.sad <- data.frame()
  # newborn_stage <- 1
  matrix_list <- list()
  dateqm <- data.frame(nprms=1:17)
  
  lrs_gamma.distr <- data.frame()
  
  dat <- t(dat)
  
  nsim <- ncol(dat)
  for(i in 1:nsim)
  {
    
    dateqm$values <- dat[,i]
    dateqm <- as.data.frame(dateqm)
    toprms <- dateqm
    # list.stuff <- get_matrix_agestage_N(datcov, N[j])[[1]] # any [[1]] except 16 is ok because the matrcies are constant
    
    s.params <- c(toprms$values[1], toprms$values[2],toprms$values[3],0,0)
    r.params <- c(toprms$values[4],toprms$values[5],toprms$values[6],0,0)
    g.params <- c(toprms$values[7],toprms$values[8],toprms$values[9],0,0,toprms$values[10],toprms$values[11],0,0,0)
    d.params <- c(toprms$values[12],toprms$values[13],toprms$values[14],0,0,toprms$values[15],toprms$values[16],0,0,0)
    # N <- toprms$values[17]
    
    N <- toprms$values[17]*(eqm_ratio)
    
    M <- bigmatrix_agestage_lrs(nn, s.params, r.params, g.params, d.params, N)
    
    ## we have 16 age classes such that survival is same for each age as it is in stage
    ## and the last age survival is 0.
    
    age.at.death <- 16
    age.stage.list <- list()
    # nn <- 50 # no. of stages
    surv.age.at.death <- rep(0, length(s.params))
    
    age.stage.list[[1]] <- bigmatrix_agestage_lrs(nn, s.params, r.params, g.params, d.params,N)
    
    for (j in 2:(age.at.death))
    {
      age.stage.list[[j]] <-     age.stage.list[[1]]
    }
    
    
    # age.stage.list[[age.at.death]] <- bigmatrix_agestage(nn, surv.age.at.death, r.params, g.params, d.params,N)
    age.stage.list[[age.at.death]]$S <- matrix(rep(0, nrow(age.stage.list[[1]]$S)*ncol(age.stage.list[[1]]$S)),
                                               nrow= nrow(age.stage.list[[1]]$S), byrow=T)
    age.stage.list[[1]]$S
    
    # create fertility and survival matrices
    matF <- matrix(NA, ncol=ncol(M$R), nrow=ncol(M$R))
    
    # # create the recruitment matrix (here only stage 1 being born into)
    # matF[1,] <- M$R[2,]
    matF[1,] <- diag(M$R)
    add.zero <- matrix(0, ncol = ncol(M$R), nrow=ncol(M$R)-1)
    matF[2:ncol(M$R),] <- add.zero
    
    #
    # this is overall matrix
    matU <- M$G %*% M$S
    mat <- M$G %*% M$S +  M$D %*% M$R
    
    # using block matrix to compare R0 etc
    Fmat <- M$D %*% M$R
    Umat <- M$G %*% M$S
    
    block_mat <- block_matrix_func(Fmat, Umat, age.at.death = 16)
    # Re(eigen(block_mat)$values[1])
    
    # faster eigen calculation!
    lambda_block <- Re(RSpectra::eigs(block_mat[[1]], 1)$values[1])
    # Re(eigen(block_mat)$values[1])
    
    mat <- as.matrix(mat)
    lambda <- Re(eigen(mat)$values[1])
    
    diff_lambda = lambda-lambda_block
    diff_lambda
    sad1 <- Mod(eigen(mat)$vector[,1])
    sad <- sad1/sum(sad1)
    
    newF <-  M$D %*% (M$R*sad)
    matrix_list[[i]] <- newF
    
  }
  # Calculate the average matrix
  average_Fmat <- Reduce(`+`, matrix_list) / nsim
  
  return(average_Fmat)
}

# Function to combine results for different densities
matF_plot_func <- function(dat){
  matF_plot <- list()
  nsim <- nrow(dat)
  # i =1
  for(i in 1:length(eqm_val))
  {
    temp <- matF_func(dat, eqm_ratio = eqm_val[i])
    
    matF_plot[[i]] <- temp
  }
  return(matF_plot)
}

# Function for survival, recruitment, growth and combined changes
# The function calculates LRS for the perturbed data
eqm_val <- c(0, 0.5, 1)
lrs_calc_func <- function(dat)
{
  lrs_plotdat <- data.frame()
  
  for(i in 1:length(eqm_val))
  {
    temp <- get_lrs_func(dat, eqm_ratio = eqm_val[i])
    temp$den <- eqm_val[i]
    lrs_plotdat <- rbind(lrs_plotdat, temp)
  }
  return(lrs_plotdat)
}

# Function for evaluating teh SSD at different population densities
sad_func  <- function(toprms){
  
  s.params <- c(toprms[1], toprms[2],toprms[3],0,0)
  r.params <- c(toprms[4],toprms[5],toprms[6],0,0)
  g.params <- c(toprms[7],toprms[8],toprms[9],0,0,toprms[10],toprms[11],0,0,0)
  d.params <- c(toprms[12],toprms[13],toprms[14],0,0,toprms[15],toprms[16],0,0,0)
  
  eqm_ratio <- c(0, 1/2, 1)
  # eqm_ratio <- c(1)
  # i=1
  sad_tbl <- data.frame()
  N <- c()
  # sadtbl1 <- data.frame()
  
  for (i in 1:length(eqm_ratio))
  {
    N <- toprms[17]*eqm_ratio[i]
    
    M <- bigmatrix_agestage_lrs(nn, s.params,r.params,g.params,d.params,N)
    
    matF <- matrix(NA, ncol=ncol(M$R), nrow=ncol(M$R))
    
    # using block matrix
    matF <- M$D %*% M$R
    matU <- M$G %*% M$S
    
    
    # mat <- M$G %*% M$S + M$D %*% matF
    mat <- M$G %*% M$S + M$D %*% M$R
    mat <- as.matrix(mat)
    eig <- eigen(mat)
    
    # lambda <- eigen(mat)$values
    sad1 <- Re(eig$vectors[,1]) # Stable size distribution
    sad <- sad1/sum(sad1)
    plot(sad)
    
    # ssdtbl1  <- cbind(sad)
    sadtbl <- as.data.frame(sad)
    sad_tbl <- rbind(sad_tbl, sadtbl)
  }
  
  sad_tbl <- as.data.frame(sad_tbl)
  
  return(sad_tbl)
}
