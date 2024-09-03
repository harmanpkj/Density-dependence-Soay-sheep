rm(list=ls())
# Harman Jaggi
# Code for Ecology Letters manuscript on Sep 2024: Density dependence

# Install and load required packages
# install.packages("")
library(RColorBrewer)
library(ggpubr)
library(tidyverse)
library(reshape2)
library(devtools)
library(ggbiplot)
library(ggfortify)
library(corrplot)
library(ggpmisc)
library(patchwork)
library(vegan)
library(magick)
library(parallel)
library(factoextra)
library(cowplot)
library(hrbrthemes)
library(Hmisc)
library(RSpectra) # for faster eigenvalues for block matrices
library('plot.matrix')
library(patchwork)
library(stats)

# check your working directory
getwd()

# Set your working directory and save the functions load_LRS_functions
# in the directory along with the file rand.params.

# Source the R Script "load_LRS_functions" from the working directory
# using the code in the line below
source("load_LRS_functions.R")

# Read the params from covariance matrix samples.
# This file contains 10000 parameter sets after delifing.
# We have attached another code for delifing in order to carry out the
# analysis for a different/new species

dat_og1 <- read.csv("rand.params.csv")

# remove the serial number column
dat_og <- dat_og1[,-1]
head(dat_og)

# Set the number of parameter sets. For the results in paper we use 250.
# This is the sample of the parameter sets we would be working with.
# Each column corresponds to number of parameters. These are 16 for our case.
# Each row corresponds to number of samples. We set nsim = 250.
set.seed(120)
nsim = 10

# Sample nsim number of rows from the original data set
sample1.dat <- sample_n(dat_og, nsim)

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
    lam <- lambda_block <- Re(eigs(block_mat[[1]], 1)$values[1])
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
    lam <- lambda_block <- Re(eigs(block_mat[[1]], 1)$values[1])
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

# Sample dataset with K attached as a new column
sample1.datK <- as.data.frame(addK_func(sample1.dat))

# Plot for equilirbrium sizes distribution
# This plot generates Figure A1 in the Appendix.
ggplot(sample1.datK, aes(x=newK))+
  geom_histogram(fill="gray", color="black", bins=15)+
  theme_bw()+
  labs(x="Distribution of K for different life-histories")+
  scale_x_continuous(limits=c(400,500))+
  theme_bw()+
  theme(
    # axis.title.y=element_blank(),
    axis.title.x=element_text(size=22),
    axis.title.y=element_text(size=22),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=13),
    legend.text = element_text(size = 22),
    strip.text.x = element_text(size=14, face="bold"))

mean(sample1.datK$newK)

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

    if(colSums(M$G %*% M$S)>1){
      count <- count+1
    }

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

# Vector for equilibrium ratio or equilibrium values
# eqm_val = 0 corresponds to population size N = 0
# eqm_val = 1/2 corresponds to population size N = K/2
# eqm_val = 1 corresponds to population size N = K

eqm_val <- c(0, 1/2, 1)
# eqm_val <- c(0, 1/2, 1, 1.2, 1.4)

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

# The command to create a data frame for the PCA analysis
# for all equilibrium ratios, all covariates and all samples (nsim)

vitals_tbl1 <- get_named_tbl(sample1.datK, eqm_val)

# Add another column with proper population densities for filtering and plotting
vitals_tbl1$den_name <- factor(vitals_tbl1$eqratio, levels = c("0", "0.5", "1"),
                               labels = c("At N = 0", "At N = K/2", "At N = K"))

vitals_tbl1$den <- vitals_tbl1$eqratio

vitals_tbl1_named <- vitals_tbl1
names(vitals_tbl1_named)

#######################
# PCA
#######################
# We perform PCA first analyze at all densities and finally combine
# the plots together for the figure
# We get Figure 1 in the manuscript from this section of the code

# PCA Covariates names for the final figure
colnames(vitals_tbl1_named) <- c("Juv. Survival", "Ad. Survival", "Reproduction", "GIs", "GIb","Growth", "Gb", "D", "PopN", "eqratio", "den_name", "den")

# Create a data frame for population density N = 0
temp_n0 <- as.data.frame(vitals_tbl1_named %>% filter(eqratio==0))
temp_n0 <- temp_n0[,c(1:3,6)]
row.names(temp_n0) <- NULL

# Carry out PCA analysis at N = 0
pca_res_n0 <- prcomp(temp_n0, scale. = TRUE)
pca_res_n0 <- princomp(temp_n0, cor = T, scores = TRUE)

# Now plot the PCA results at N = 0
pca_plot_n0 <- autoplot(pca_res_n0, data = temp_n0, loadings.label = TRUE,
                        loadings = TRUE, loadings.colour = 'black',
                        loadings.label.size = 4.5,
                        alpha=0,
                        loadings.label.vjust = 0.8,
                        loadings.label.hjust = 0.7,
                        loadings.label.colour="cyan4")+
  theme_bw()+
  labs(title = "At N=0")+
  theme(
    plot.title = element_text(size=18))+
  geom_point(color = "orange2",  size = 1, alpha=0.55)


# Create a data frame for population density N = K/2
temp_k2 <- as.data.frame( vitals_tbl1_named %>% filter(eqratio==0.5))
temp_k2 <- temp_k2[,c(1:3,6)]
row.names(temp_k2) <- NULL

# Carry out PCA analysis at N = K/2
pca_res_k2 <- prcomp(temp_k2, scale. = TRUE)
pca_res_k2 <- princomp(temp_k2, cor = T, scores = TRUE)

# Now plot the PCA results at N = K/2
pca_plot_k2 <-autoplot(pca_res_k2, data = temp_k2, loadings.label = TRUE,
                       loadings = TRUE, loadings.colour = 'black',
                       loadings.label.size = 4.5,
                       alpha=0,
                       loadings.label.vjust = 1,
                       loadings.label.hjust = 0.7,
                       loadings.label.colour="cyan4")+
  theme_bw()+
  labs(title = "At N=K/2")+
  theme(
    plot.title = element_text(size=18))+
  geom_point(color = "orange2",  size = 1, alpha=0.55)

# Create a data frame for population density N = K
temp_k <- as.data.frame( vitals_tbl1_named %>% filter(eqratio==1))
temp_k <- temp_k[,c(1:3,6)]
row.names(temp_k) <- NULL

# Carry out PCA analysis at N = K
pca_res_k <- prcomp(temp_k, scale. = TRUE)
pca_res_k <- princomp(temp_k, cor = T, scores = TRUE)

# Now plot the PCA results at N = K
pca_plot_k <-autoplot(pca_res_k, data = temp_k,
                      loadings.colour = 'black',
                      # color=NULL,
                      loadings = TRUE,
                      alpha=0,
                      loadings.label = TRUE,
                      geom = "arrow",
                      loadings.label.size = 4.5,
                      loadings.label.vjust = 1,
                      loadings.label.hjust = 1,
                      loadings.label.colour="cyan4"
)+
  theme_bw()+
  labs(title = "At N=K")+
  theme(
    plot.title = element_text(size=18))+
  geom_point(color = "orange2",  size = 1, alpha=0.55)

# This command makes for Figure 1 in the manuscript
ggarrange(pca_plot_n0, pca_plot_k2, pca_plot_k, nrow=3)

#########################################################
# Repeat the analysis for densities higher than N = K.
# In the paper, we analyze for N = 1.2 K and 1.4 K.
# This commented section results in the right panels of Figure A7
#########################################################
# eqm_val <- c(0, 0.5, 1, 1.2, 1.4)
# vitals_tbl1 <- get_named_tbl(sample1.datK, eqm_val)
# vitals_tbl1$den_name <- factor(vitals_tbl1$eqratio, levels = c("0", "0.5", "1"),
#                               labels = c("At N = 0", "At N = K/2", "At N = K"))
# vitals_tbl1$den <- vitals_tbl1$eqratio
# vitals_tbl1_named <- vitals_tbl1
# names(vitals_tbl1_named)
# temp_2k <- as.data.frame( vitals_tbl1_named %>% filter(eqratio==1.2))
# temp_2k <- temp_2k[,c(1:3,6)]
# row.names(temp_2k) <- NULL

# tbl.pca_dep <- prcomp(temp, scale. = TRUE, center = T)

# pca_res_2k <- prcomp(temp_2k, scale. = TRUE, center = T)
# pca_res_2k <- princomp(temp_2k, cor = T, scores = TRUE)
# pca_plot_2k <-autoplot(pca_res_2k, data = temp_k,
#                        loadings.colour = 'black',
#                        loadings.size =3,
#                        # color=NULL,
#                        alpha=0,
#                        loadings.label = TRUE,
#                        geom = "arrow",
#                        loadings.label.size = 4.5,
#                        loadings.label.vjust = 1,
#                        loadings.label.hjust = 1,
#                        loadings.label.colour="cyan4"
# )+
#   theme_bw()+
#   labs(title = "At N=(1.2)K")+
#   theme(
#     plot.title = element_text(size=18))+
#   geom_point(color = "orange2",  size = 1, alpha=0.4)
#
#
# temp_2k1 <- as.data.frame( vitals_tbl1_named %>% filter(eqratio==1.4))
# temp_2k1 <- temp_2k1[,c(1:3,6)]
# row.names(temp_2k1) <- NULL
#
# # tbl.pca_dep <- prcomp(temp, scale. = TRUE, center = T)
#
# pca_res_2k1 <- prcomp(temp_2k1, scale. = TRUE, center = T)
# pca_res_2k1 <- princomp(temp_2k1, cor = T, scores = TRUE)
# pca_plot_2k1 <-autoplot(pca_res_2k1, data = temp_k,
#                         loadings.colour = 'black',
#                         loadings.size =3,
#                         # color=NULL,
#                         alpha=0,
#                         loadings.label = TRUE,
#                         geom = "arrow",
#                         loadings.label.size = 4.5,
#                         loadings.label.vjust = 1,
#                         loadings.label.hjust = 1,
#                         loadings.label.colour="cyan4"
# )+
#   theme_bw()+
#   labs(title = "At N=(1.4)K")+
#   theme(
#     plot.title = element_text(size=18))+
#   geom_point(color = "orange2",  size = 1, alpha=0.4)
# ggarrange(pca_plot_2k, pca_plot_2k1, nrow=2)

# Results from the PCA table that lead to Figure A3 in the manuscript
d1 <- as.data.frame(pca_res_n0$rotation)
d1$N <- rep("0", nrow(d1))
d2 <- as.data.frame(pca_res_k2$rotation)
d2$N <- rep("K/2", nrow(d2))
d3 <- as.data.frame(pca_res_k$rotation)
d3$N <- rep("K", nrow(d3))

pca_tbl_save <- rbind(d1, d2, d3)
pca_res_k$loadings

# Save the table
# write.csv(pca_tbl_save, "pca_tbl_save.csv")

# Code for broken-stick model to compare if PCA covariates are significant
# Figure A5 in the manuscript
# Rerun the lines by changing the dat for different densities
par(mar = c(0, 0, 0, 0))
dat <- vitals_tbl1 %>% filter(eqratio==0)
dat <- vitals_tbl1 %>% filter(eqratio==0.5)
dat <- vitals_tbl1 %>% filter(eqratio==1)

temp <- as.data.frame(dat)
temp <- temp[,1:5]
row.names(temp) <- NULL
# as.matrix(temp)
tbl.pca_dep <- prcomp(temp, scale. = TRUE, center = T)
pca_result <- tbl.pca_dep

# the observed eigenvalues from PCA
observed_eigenvalues <- pca_result$sdev^2

n_pca <- ncol(temp)

# Calculate broken stick values
broken_stick_values <- sapply(1:n_pca, function(k) sum(1/(k:n_pca))/n_pca)

# eigenvalues from PCA
eigenvalues <- pca_result$sdev^2

# Now compare observed eigenvalues to broken stick values
plot(1:n_pca, eigenvalues, type = 'b', pch = 19, col = 'blue', xlab = "Component", ylab = "Eigenvalue", main = "PCA Eigenvalues vs. Broken Stick")
lines(1:n_pca, broken_stick_values, type = 'b', pch = 17, col = 'red')
legend("topright", legend = c("Observed Eigenvalues", "Broken Stick"), col = c("blue", "red"), pch = c(19, 17))

# Commands for plotting the correlation matrices between covariates
# This line of code plots Figure A4 and left panel of Figure A7
corrdat1_N0 <- filter(vitals_tbl1, den==0)[,c(1:3,6,7)]
colnames(corrdat1_N0) <- c("Ss", "Sb", "R", "GIs", "D")

corrdat1_K_2 <- filter(vitals_tbl1, den==0.5)[,c(1:3,6,7)]
colnames(corrdat1_K_2) <- c("Ss", "Sb", "R", "GIs", "D")

corrdat1_K <- filter(vitals_tbl1, den==1)[,c(1:3,6,7)]
colnames(corrdat1_K) <- c("Ss", "Sb", "R", "GIs", "D")

cor_mat_sam1_N0 <- cor(corrdat1_N0)
cor_mat_sam1_K_2 <- cor(corrdat1_K_2)
cor_mat_sam1_K <- cor(corrdat1_K)

# Figure A4 and A7 in the manuscript
# pdf("corplot_sam1_ratio.pdf")
par(mfrow=c(3,1))
corrplot(cor_mat_sam1_N0, method="color", addCoef.col = "black", number.cex=0.8)
corrplot(cor_mat_sam1_K_2, method="color", addCoef.col = "black", number.cex=0.8)
corrplot(cor_mat_sam1_K, method="color", addCoef.col = "black", number.cex=0.8)
# dev.off()
# Run the above code for eqm values of 1.2K and 1.4K to get the correlatin matrices
# shown in the left panel of Figure A7

# pdf("corplot_sam1_beyondK.pdf")
# par(mfrow=c(1,2))
# corrplot(cor_mat_sam1_2K1, method="color", addCoef.col = "black", number.cex=0.8)
# # corrplot(cor_mat_sam1_2K2, method="color", addCoef.col = "black", number.cex=0.8)
# corrplot(cor_mat_sam1_2K3, method="color", addCoef.col = "black", number.cex=0.8)
# # corrplot(cor_mat_sam1_2K4, method="color", addCoef.col = "black", number.cex=0.8)
# dev.off()


# Plot for Figure A6 in the manuscript
# Plots the relationship between Average juvenile survival and Growth at different densities
# along with their slopes
ggplot(filter(vitals_tbl1), aes(y=Ss, x=Gspr, color=den_name, group=den_name)) +
  geom_point( alpha=1, size=1.5) +
  geom_smooth(method='lm', formula= y~x, se=F,
              size=0.5, linetype="dashed",
              aes(group = den_name, color=den_name))+
  scale_fill_discrete(name = "Dose", labels = c("A", "B", "C", "D", "E"))+
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
               parse = TRUE,
               rr.digits = 2, coef.digits = 2, size = 3, label.x = 0.97) +
theme_bw()+
  labs(y="Average Juvenile Survival at SSD",
       x="Average Growth for juveniles at SSD")+
  theme(axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        axis.text.y=element_text(size=20),
        axis.text.x=element_text(size=13),
        legend.title = element_text( size = 16),
        legend.text = element_text(size = 17),
        legend.position = "right")+
  guides(color = guide_legend(title = "Population N"))


###############################################
# LRS
###############################################

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
    lambda_block <- Re(eigs(block_mat[[1]], 1)$values[1])
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

# This section plots mean LRS at different populations sizes
# for all parameter sets.
# Results in Figure 2 in the manuscript
dat <- sample1.dat

# Set the population size of interest.
N_func <- N_num <- c(0, 100, 200, 300, 400, 500)
lrs_gamma.distr_N <- data.frame()
average_Fmat_list <- list()
for (k in 1:length(N_func))
{
  dat <- sample1.dat

  tbl <- data.frame()
  tblnew <- data.frame()
  tblnew.sad <- data.frame()
  # newborn_stage <- 1

  dateqm <- data.frame(nprms=1:16)

  lrs_gamma.distr <- data.frame()

  dat <- t(dat)

  nsim <- ncol(dat)
  matrix_list <- list()

  for(i in 1:nsim)
  {
    N <- N_func[k]
    print(k)
    dateqm$values <- dat[,i]
    dateqm <- as.data.frame(dateqm)
    toprms <- dateqm
    # test for different values of N: density- independence and density-dependent

    s.params <- c(toprms$values[1], toprms$values[2],toprms$values[3],0,0)
    r.params <- c(toprms$values[4],toprms$values[5],toprms$values[6],0,0)
    g.params <- c(toprms$values[7],toprms$values[8],toprms$values[9],0,0,toprms$values[10],toprms$values[11],0,0,0)
    d.params <- c(toprms$values[12],toprms$values[13],toprms$values[14],0,0,toprms$values[15],toprms$values[16],0,0,0)

    M <- bigmatrix_agestage_lrs(nn, s.params, r.params, g.params, d.params, N)

    ## we have 16 age classes such that survival is same for each age as it is in stage
    ## and the last age survival is 0.

    age.at.death <- 16
    age.stage.list <- list()
    # nn <- 50 # no. of stages
    # TEST FOR DENTIY
    surv.age.at.death <- rep(0, length(s.params))

    age.stage.list[[1]] <- bigmatrix_agestage_lrs(nn, s.params, r.params, g.params, d.params, N)

    for (j in 2:(age.at.death))
    {
      age.stage.list[[j]] <-     age.stage.list[[1]]
    }

    # age.stage.list[[age.at.death]] <- bigmatrix_agestage(nn, surv.age.at.death, r.params, g.params, d.params,N)
    age.stage.list[[age.at.death]]$S <- matrix(rep(0, nrow(age.stage.list[[1]]$S)*nrow(age.stage.list[[1]]$S)),
                                               nrow= nrow(age.stage.list[[1]]$S), byrow=T)

    # create fertility and survival matrices
    matF <- matrix(NA, ncol=ncol(M$R), nrow=ncol(M$R))

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
    lambda_block <- Re(eigs(block_mat[[1]], 1)$values[1])

    mat <- as.matrix(mat)
    lambda <- Re(eigen(mat)$values[1])

    diff_lambda = lambda-lambda_block
    diff_lambda
    sad1 <- Mod(eigen(mat)$vector[,1])
    sad <- sad1/sum(sad1)
    # sum(sad)
    # M$D %*% matF[1,]*sad

    newF <-  (M$D %*% M$R)
    matrix_list[[i]] <- newF
    # sum(sad)
    # plot(newF)

    # newF <-  M$D %*% diag(M$R) *(sad)
    newborns_dist1 <- colSums(t(newF))

    # normalize the newborn distribution
    newborns_dist <- newborns_dist1/sum(newborns_dist1)

    matrices <- age.stage.list
    n.stages <- dim(matrices[[1]]$R)[2] # obtain number of stages
    kappamat <- r_to_kappamat1(matrices) # list of n.stages matrix
    G <- g_to_gmat(matrices) # list of n.stages matrix
    pvec <- s_to_pvec(matrices) #size conditional survival, n.stages stages, row is corresponded to stage
    max.offspring <- 30

    ptm <- proc.time()

    gamma <- block_matrices_solve(Klist = kappamat, Glist = G, Pmatrix = pvec,
                                  max.offspring = max.offspring,
                                  initial.age = 1,
                                  end.age = 1, # not 2 anymore mother born each stage then sum
                                  stage_per_age = n.stages)

    new_gamma <- gamma %*% newborns_dist
    print(paste0(ncol(dat), ",", i))
    #plot the new LRS distribution weighted by new born distribution
    temp <- data.frame(density = new_gamma,
                       mids = seq_along(gamma[,1]) - 1)

    # mean and variance from the lrs distribution
    junk1 <- temp %>% summarise(R0=sum(mids*density),
                                var=sum(mids*density - sum(mids*density))^2,
                                dd=sum(density))

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
    # View(U_block)
    # dim(N_block)
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
    temp$N_val <- rep(N_func[k], nrow(temp))
    lrs_gamma.distr <- rbind(lrs_gamma.distr, temp)
  }

  # Calculate the average matrix
  average_Fmat <- Reduce(`+`, matrix_list) / nsim

  average_Fmat_list[[k]] <- average_Fmat

  lrs_gamma.distr_N <-    rbind(lrs_gamma.distr_N, lrs_gamma.distr)
}

# We remove the package for summarize to work properly
detach(package:ggbiplot)
detach(package:plyr)

mean_var_R0_at_N <- lrs_gamma.distr_N %>% group_by(N_val) %>%
  mutate(mean_mean_lrs=mean(mean.lrs.dist),var_mean_lrs=var(mean.lrs.dist))

unique(mean_var_R0_at_N$var_mean_lrs)
unique(mean_var_R0_at_N$var_mean_lrs/mean_var_R0_at_N$mean_mean_lrs)

# Command for Figure 2 in the manuscript
ggplot() +
  geom_point(lrs_gamma.distr_N, mapping=aes(x=N_val, y=mean.lrs.dist), size=1.7) +
  labs(y="R0", x="Population abundance: N")+
  geom_point(data=mean_var_R0_at_N, mapping=aes(x = N_val, y = mean_mean_lrs),
             col="dark red", size=4,
             shape=2)+
  theme_bw()+
  # labs(x="Juvenile Survival scaled by SSD")+
  theme(
    # axis.title.y=element_blank(),
    axis.title.x=element_text(size=22),
    axis.title.y=element_text(size=22),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=13),
    # strip.text.x = element_text(size = 15),
    # legend.title = element_text( size = 22),
    legend.text = element_text(size = 22),
    strip.text.x = element_text(size=14, face="bold"))+
  scale_y_continuous(breaks = seq(0, 100, by = 1))

# Save the file from above code
# write.csv(lrs_gamma.distr_N, "lrs_gamma.distr_N.csv")

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
    lambda_block <- Re(eigs(block_mat[[1]], 1)$values[1])
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

# Execute the function to get results for mother's distribution
matF_return <- matF_plot_func(sample1.datK)

# Combine matrices into one data frame
combined_df <- do.call(rbind, lapply(matF_return, function(mat) {
  melt(mat)
}))

# Overall range of values for colors
overall_range <- c(0, max(combined_df$value))

# Create a color palette
custom_palette <- c("gray90", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7",
                    "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061")


# Command for Figure 3 in the manuscript
# Plot each matrix with the same color scale breaks and custom palette
mother_plots <- lapply(matF_return, function(mat) {
  df <- melt(mat)
  # df <- df[longData1$value>0.01,]

  ggplot(df, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill = value)) +
    scale_fill_gradientn(colors = custom_palette, limits = overall_range) +  # Use the same color scale limits and custom palette
    labs(x="", y="") +
    scale_fill_distiller(palette = "Spectral",
                         limits = overall_range
                         # , values = c(0,0.05,0.1, 0.5,1)
    )+
    theme_classic() +
    theme(
      text = element_text(size = 20),
      panel.border = element_rect(colour = "black", fill = FALSE),
      aspect.ratio = 1,
      # legend.key.size = element_text(size=0.5),
      legend.position = "right"
    )
})

# Command for Figure 3 in the manuscript
combined_mother_plot <- ggarrange(mother_plots[[1]],
                                  mother_plots[[2]],
                                  mother_plots[[3]], nrow=1, common.legend = T, legend="right")
annotate_figure(combined_mother_plot, left = textGrob("Distribution of Offspring size", rot = 90, vjust = 1, gp = gpar(cex = 1.4)),
                bottom = textGrob("Distribution of Mothers size (scaled by SSD)", vjust=-6, gp = gpar(cex = 1.4)))

# combined_mother_plot <- mother_plots[[1]] +
#   mother_plots[[2]] +
#   mother_plots[[3]] & theme(legend.position = "right")
# combined_mother_plot + plot_layout(guides = "collect")

n_stage <- nn
stages <- 1:n_stage
xx <- colSums(matF_return[[1]])
xxx <- xx/sum(xx)
weight <- sum(xxx*stages)

sum((colSums(matF_return[[1]])/sum(colSums(matF_return[[1]])))*(stages*maxsize/n_stage))
sum((colSums(matF_return[[2]])/sum(colSums(matF_return[[2]])))*(stages*maxsize/n_stage))
sum((colSums(matF_return[[3]])/sum(colSums(matF_return[[3]])))*(stages*maxsize/n_stage))

sum((rowSums(matF_return[[1]])/sum(colSums(matF_return[[1]])))*(stages*maxsize/n_stage))
sum((rowSums(matF_return[[2]])/sum(colSums(matF_return[[2]])))*(stages*maxsize/n_stage))
sum((rowSums(matF_return[[3]])/sum(colSums(matF_return[[3]])))*(stages*maxsize/n_stage))

mean(colSums(matF_return[[2]]))
mean(colSums(matF_return[[3]]))

# Command for Figure 4 and Figure A9
# Distribution of offspring and mother body size with density

plot(y=colSums(matF_return[[1]]), x=seq(0.1,38, length=50),
     pch=20, xlab="Body size (in kg)", ylab="Proportion of mothers", xaxt='n')
axis(1, at = seq(0, 38, by=4))
lines(y=colSums(matF_return[[1]]), x=seq(0.1,38, length=50))
abline(v=sum((colSums(matF_return[[1]])/sum(colSums(matF_return[[1]])))*stages*maxsize/n_stage), lty=2)
points(y=colSums(matF_return[[2]]), x=seq(0.1,38, length=50), pch=20, col="cyan4", lty=2, xaxt='n')
lines(y=colSums(matF_return[[2]]), x=seq(0.1,38, length=50), col="cyan4")
abline(v=sum((colSums(matF_return[[2]])/sum(colSums(matF_return[[2]])))*stages*maxsize/n_stage), col="cyan4", lty=2, lwd=1.5)
points(y=colSums(matF_return[[3]]), x=seq(0.1,38, length=50), pch=20, col="purple4")
lines(y=colSums(matF_return[[3]]), x=seq(0.1,38, length=50), pch=20, col="purple4")
abline(v=sum((colSums(matF_return[[3]])/sum(colSums(matF_return[[3]])))*stages*maxsize/n_stage), col="purple4", lty=2, lwd=1.5)
legend(32, 0.03, legend=c("N = 0", "N = K/2", "N = K"),
       col=c("black","cyan4", "purple4"), lty=1, cex=0.9)


plot(y=rowSums(matF_return[[1]]), x=seq(0.1,38, length=50), pch=20, xlab="Body size (in kg)", ylab="Proportion of offspring", xaxt='n')
axis(1, at = seq(0, 38, by=4))
lines(y=rowSums(matF_return[[1]]), x=seq(0.1,38, length=50))
abline(v=sum((rowSums(matF_return[[1]])/sum(rowSums(matF_return[[1]])))*stages*maxsize/n_stage), lty=2)
points(y=rowSums(matF_return[[2]]), x=seq(0.1,38, length=50), pch=20, col="cyan4", lty=2, xaxt='n')
lines(y=rowSums(matF_return[[2]]), x=seq(0.1,38, length=50), col="cyan4")
abline(v=sum((rowSums(matF_return[[2]])/sum(rowSums(matF_return[[2]])))*stages*maxsize/n_stage), col="cyan4", lty=2, lwd=1.5)
points(y=rowSums(matF_return[[3]]), x=seq(0.1,38, length=50), pch=20, col="purple4")
lines(y=rowSums(matF_return[[3]]), x=seq(0.1,38, length=50), pch=20, col="purple4")
abline(v=sum((rowSums(matF_return[[3]])/sum(rowSums(matF_return[[3]])))*stages*maxsize/n_stage), col="purple4", lty=2, lwd=1.5)
legend(28, 0.05, legend=c("N = 0", "N = K/2", "N = K"),
       col=c("black","cyan4", "purple4"), lty=1, cex=0.9)

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

sad_res1 <- apply(sample1.datK, 1, sad_func)
sad_res <- as.data.frame(sad_res1)
dim(sad_res)

colnames(sad_res) <- c(1:nsim)

df <- data.frame()
df1<- data.frame(Sno=1:nrow(sad_res))
dim(sad_res)

for (i in 1:(length(sad_res)))
{
  df1$sad <-  sad_res[,i]
  df1 <- as.data.frame(df1)
  df1$sim <- rep(i, nrow(df1))
  df1$size <- rep(seq(1, 50, 1),3)
  df1$eqm <- c(rep(0, 50), rep(0.5, 50), rep(1, 50))
  df <- rbind(df, df1)
}

names(df)
detach(package:ggbiplot)
detach(package:plyr)
library(dplyr)

nrow(df)

# Command for Figure A8 in the appendix: plots SSD

df %>% group_by(eqm, size) %>% summarise(mean_sad=mean(sad)) %>%
  ggplot(., aes(x=size , y=mean_sad, group=as.factor(eqm), color=as.factor(eqm)))+
  geom_line(size=1.6)+
  geom_point(size=1.2)+
  # facet_wrap(.~eqm)+
  # xlim(c(2,50))+
  theme_bw()+
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"))+
  labs(y="Stable size distribution",
       x="Stages for body size",
       color="N")+
  theme(axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        axis.text.y=element_text(size=20),
        axis.text.x=element_text(size=13),
        # strip.text.x = element_text(size = 15),
        legend.title = element_text( size = 18),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size=15, face="bold"))+

  # scale_x_continuous(breaks = seq(0, 40, by = 4), limits=c(0,40))+
  geom_vline(xintercept=26.06,lwd=0.7,colour="#009E73", linetype="dashed")+
  geom_vline(xintercept=26.79,lwd=0.7,colour="#56B4E9", linetype="dashed")+
  geom_vline(xintercept=28.12,lwd=0.7,colour="#E69F00", linetype="dashed")+
  geom_vline(xintercept=17.78,lwd=0.7,colour="#009E73", linetype="twodash")+
  geom_vline(xintercept=18.24,lwd=0.7,colour="#56B4E9", linetype="twodash")+
  geom_vline(xintercept=19.21,lwd=0.7,colour="#E69F00", linetype="twodash")

df %>% group_by(eqm, size) %>% summarise(mean_sad=mean(sad)) %>%
  group_by(eqm) %>%
  summarise(mean_eqm=mean(mean_sad))

n <- nn
df <- df %>% mutate(bsize= rep(minsize+c(1:n)*(maxsize-minsize)/n, nsim*length(eqm_val)))
# Filter data based on eqm values
df_eqm_0 <- df %>% filter(eqm == 0)
df_eqm_05 <- df %>% filter(eqm == 0.5)
df_eqm_1 <- df %>% filter(eqm == 1)

# Calculate the average 'sad' value for each size for each eqm value
average_sad_eqm_0 <- df_eqm_0 %>% group_by(bsize) %>% summarise(sad = mean(sad))
average_sad_eqm_05 <- df_eqm_05 %>% group_by(bsize) %>% summarise(sad = mean(sad))
average_sad_eqm_1 <- df_eqm_1 %>% group_by(bsize) %>% summarise(sad = mean(sad))

# Define the average sizes for the vertical lines
v1 <- minsize+unique(vitals_tbl1$aver_juv_size)*(maxsize-minsize)/n
v2 <- minsize+unique(vitals_tbl1$aver_ad_size)*(maxsize-minsize)/n
new_average_values <- c(v1[1], v1[2], v1[3])
second_new_average_values <- c(v2[1], v2[2], v2[3])

minsize+unique(vitals_tbl1$aver_ad_size)*(maxsize-minsize)/n

colors <- c("N = 0" = "#E69F00", "N = K/2" = "#56B4E9", "N = K" = "#009E73")

# Create the plot
ggplot() +
  geom_point(data = df_eqm_0, aes(x = bsize, y = sad), alpha = 0.3, color = "#E69F00") +
  geom_line(data = average_sad_eqm_0, aes(x = bsize, y = sad), color = "#E69F00", size = 1.3, ) +
  geom_point(data = df_eqm_05, aes(x = bsize, y = sad), alpha = 0.3, color = "#56B4E9") +
  geom_line(data = average_sad_eqm_05, aes(x = bsize, y = sad), color = "#56B4E9", size = 1.3) +
  geom_point(data = df_eqm_1, aes(x = bsize, y = sad), alpha = 0.3, color = "#009E73") +
  geom_line(data = average_sad_eqm_1, aes(x = bsize, y = sad), color = "#009E73", size = 1.3) +
  geom_vline(xintercept = new_average_values[1], color = "#E69F00", linetype = "dashed") +
  geom_vline(xintercept = second_new_average_values[1], color = "#E69F00", linetype = "twodash") +
  geom_vline(xintercept = new_average_values[2], color = "#56B4E9", linetype = "dashed") +
  geom_vline(xintercept = second_new_average_values[2], color = "#56B4E9", linetype = "twodash") +
  geom_vline(xintercept = new_average_values[3], color = "#009E73", linetype = "dashed") +
  geom_vline(xintercept = second_new_average_values[3], color = "#009E73", linetype = "twodash") +
  scale_x_continuous(breaks = seq(0,38,length=6), labels = c(0.1, 8, 16, 24, 32, 38)) +
  labs(x = 'Body size', y = 'Stable stage distribution', color= "sth")+
  scale_color_manual(values = colors)+
  # scale_color_manual(name="average", values = c("N = 0" = "#E69F00", "N = K/2" = "#56B4E9", "N = K" = "#009E73")) +
  theme_classic()

# Code for overall LRS distribution plots
eqm_val <- c(0, 0.5, 1)
sam1_meanlrs <- meanlrs_vital_func(sample1.datK)
sam1_meanlrs$den_name <- factor(sam1_meanlrs$den, levels = c("0", "0.5", "1"),
                                labels = c("At N = 0", "At N = K/2", "At N = K"))

# nrow(sam1_meanlrs)

# Merge the vitals and meanlrs datasets used for PCA analysis and LRS analysis respectively
dat1_vitals_meanlrs <- cbind(vitals_tbl1, sam1_meanlrs)
dat1_vitals_meanlrs <- data.frame(dat1_vitals_meanlrs)

x <- filter(dat1_vitals_meanlrs, den==1)$PopN
y1 <- filter(dat1_vitals_meanlrs, den==0)$mean_lrs

# R0 result at density independence
# Figure A10 in the manuscript
ggplot(tibble(x,y1), aes(x,y1))+
  geom_point(alpha=0.9, size=2.5) +
  theme_bw()+
  labs(x="Carrying Capacity", y="R0 at density-independence (N=0)")+
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
               parse = TRUE, color = "darkblue",
               rr.digits = 2, coef.digits = 2, size = 5)+
  geom_smooth(method='lm', formula= y~x, se=F,
              size=0.3, linetype="dashed", color="red" )+
  theme(axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        axis.text.y=element_text(size=20),
        axis.text.x=element_text(size=13),
        # strip.text.x = element_text(size = 15),
        # legend.title = element_text( size = 22),
        legend.text = element_text(size = 22),
        strip.text.x = element_text(size=15, face="bold"),
        legend.position = "bottom")+
  guides(fill=guide_legend(title="New Legend Title"))

# Command for plotting gamma, that is, the probability of having no offspring,
# with average covariate functions

a <- ggplot(dat1_vitals_meanlrs, aes(Ss, plrs0))+
  geom_point(size=2.5, alpha=0.9, color="darkolivegreen3")+
  # geom_point(dat1_vitals_meanlrs, mapping=aes(Survival_small_ssd, mean_lrs), color="red", size=2.5)+
  # scale_y_continuous(sec.axis = sec_axis(~./10, name="R0")) +
  facet_wrap(~den_name)+
  theme_bw()+
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
               parse = TRUE, color = "blue",
               rr.digits = 2, coef.digits = 2, size = 3) +

  labs(y="Probability that LRS is 0",
       x="Average juvenile survival at ssd")+
  # stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
  #              parse = TRUE, color = "blue",
  #              rr.digits = 2, coef.digits = 2, size = 3) +
  geom_smooth(method='lm', formula= y~x, se=F, color="black", size=0.5, linetype="dashed")+
  theme(axis.title.y=element_text(size=20),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.text.x=element_text(size=13),
        strip.text.x = element_text(size = 20),
        strip.text = element_text(face = "bold"))+
  scale_y_continuous(breaks = seq(0, 1.1, by = 0.10))

b <- ggplot(data.frame(dat1_vitals_meanlrs), aes(x=Sb, y=mean_lrsnot0)) +

  geom_point( alpha=0.8, size=2.5, color="cyan4") +
  # stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
  #              parse = TRUE, color = "blue",
  #              rr.digits = 2, coef.digits = 2, size = 3) +
  labs(color="")+
  facet_wrap(~den_name)+
  geom_smooth(method='lm', formula= y~x, se=F,
              size=0.5, linetype="dashed", color="black")+
  # color="blue"  )+
  # scale_colour_manual(values = c("cyan4"))+
  theme_bw()+
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
               parse = TRUE, color = "blue",
               rr.digits = 2, coef.digits = 2, size = 3) +

  labs(x="Average Juvenile Survival at SSD", y = expression(gamma)) +
  # scale_color_manual(labels = c("R0", "R0 with non-zero offspring"), values = c("cyan4", "darkolivegreen3"))+
  # scale_y_continuous(limits = c(0, 11), breaks = seq(0, 11, by = 2))+
  theme(
    # axis.title.y=element_blank(),
    axis.title.y=element_text(size=20),
    axis.title.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=13),
    # strip.text.x = element_text(size = 15),
    # legend.title = element_text( size = 22),
    legend.text = element_text(size = 22),
    # strip.text.x = element_blank(),
    legend.position = "bottom",
    strip.text.x = element_text(size = 20),
    strip.text = element_text())

names(dat1_vitals_meanlrs)
d <- ggplot(data.frame(dat1_vitals_meanlrs), aes(x=R, y=mean_lrsnot0)) +

  geom_point( alpha=0.8, size=2.5, color="darkolivegreen3") +
  # stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
  #              parse = TRUE, color = "blue",
  #              rr.digits = 2, coef.digits = 2, size = 3) +
  labs(color="")+
  facet_wrap(~den_name)+
  geom_smooth(method='lm', formula= y~x, se=F,
              size=0.5, linetype="dashed", color="black")+
  # color="blue"  )+
  # scale_colour_manual(values = c("cyan4"))+
  theme_bw()+
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
               parse = TRUE, color = "blue",
               rr.digits = 2, coef.digits = 2, size = 3) +

  labs(x="Average Reproduction at SSD", y = expression(gamma)) +
  # scale_color_manual(labels = c("R0", "R0 with non-zero offspring"), values = c("cyan4", "darkolivegreen3"))+
  # scale_y_continuous(limits = c(0, 11), breaks = seq(0, 11, by = 2))+
  theme(
    # axis.title.y=element_blank(),
    axis.title.y=element_text(size=20),
    axis.title.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=13),
    # strip.text.x = element_text(size = 15),
    # legend.title = element_text( size = 22),
    legend.text = element_text(size = 22),
    # strip.text.x = element_blank(),
    legend.position = "bottom",
    strip.text.x = element_text(size = 20),
    strip.text = element_text())

e <- ggplot(data.frame(dat1_vitals_meanlrs), aes(x=Gspr, y=mean_lrsnot0)) +

  geom_point( alpha=0.8, size=2.5, color="darkslateblue") +
  # stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
  #              parse = TRUE, color = "blue",
  #              rr.digits = 2, coef.digits = 2, size = 3) +
  labs(color="")+
  facet_wrap(~den_name)+
  geom_smooth(method='lm', formula= y~x, se=F,
              size=0.5, linetype="dashed", color="black")+
  # color="blue"  )+
  # scale_colour_manual(values = c("cyan4"))+
  theme_bw()+
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
               parse = TRUE, color = "blue",
               rr.digits = 2, coef.digits = 2, size = 3) +

  labs(x="Average Growth at SSD", y = expression(gamma)) +
  # scale_color_manual(labels = c("R0", "R0 with non-zero offspring"), values = c("cyan4", "darkolivegreen3"))+
  # scale_y_continuous(limits = c(0, 11), breaks = seq(0, 11, by = 2))+
  theme(
    # axis.title.y=element_blank(),
    axis.title.y=element_text(size=20),
    axis.title.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=13),
    # strip.text.x = element_text(size = 15),
    # legend.title = element_text( size = 22),
    legend.text = element_text(size = 22),
    # strip.text.x = element_blank(),
    legend.position = "bottom",
    strip.text.x = element_text(size = 20),
    strip.text = element_text())

# Plot for Figure A11 in the manuscript
print(b/d/e)


# add_func_sample <- function(a, b)
# {
#   cc1 <- c()
#   for(i in 1:ncol(a))
#   {
#     cc <- as.numeric(b[,i]) + as.numeric(a[,i])
#     cc1 <- cbind(cc1, cc)
#   }
#   return(cc1)
# }


# Vital rate functions and mean lrs

# In the next line of code, the goal is to perturb only those parameter values that
# correspond to survival function, or reproduction or growth, one at a time.
# We also perturb all functions combined and compare LRS results

# First we find the mean of all parameter sets so that we can create deviations from that mean
# Mean from 10000 params (one way to do it!)
mean.prm <- c()
head(dat_og)
dim(dat_og)

for(i in 1:ncol(dat_og))
{mean.prm[i] <- mean(dat_og[,i])
}
mean.prm <- as.matrix(mean.prm)

# mean.prm <- as.data.frame(mean.prm)
nsim = 10

mean.prm_K <- eqm_func(as.matrix(mean.prm))
dat <- as.matrix(c(as.matrix(mean.prm), as.matrix(mean.prm_K)))

# seed 120
eqm_val <- c(0,0.5,1)

# Matrix of mean params where each col is sims and each row is mean of params
mean.prm.mat <- as.data.frame(matrix(as.matrix(mean.prm), nrow=nrow(mean.prm), ncol=nsim))
dim(mean.prm.mat)

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

# Make new params by adding to sample taken from data. Distrub only survival.
sample1.surv <- t(rbind(t(sample1.dat)[1:3,],
                        mean.prm.mat[4:nrow(t(sample1.dat)),]))

# Make new params by adding to sample taken from data. Distrub only recruitment.
sample1.rec <- t(rbind(mean.prm.mat[1:3,],
                       t(sample1.dat)[4:6,],
                       mean.prm.mat[7:nrow(t(sample1.dat)),]))

# Make new params by adding to sample taken from data. Distrub only growth.
sample1.gro <- t(rbind(mean.prm.mat[1:6,],
                       t(sample1.dat)[7:9,],
                       mean.prm.mat[10:nrow(t(sample1.dat)),]))

# Make new params by adding to sample taken from data. Distrub both survival and growth.
sample1.surv.gro <- t(rbind(t(sample1.dat)[1:3,],
                            mean.prm.mat[4:6,],
                            t(sample1.dat)[7:9,],
                            mean.prm.mat[10:nrow(t(sample1.dat)),]))

# Structure works!
sample1.test <- t(rbind(t(sample1.dat)[1:3,],
                        t(sample1.dat)[4:6,],
                        t(sample1.dat)[7:9,],
                        t(sample1.dat)[10:nrow(t(sample1.dat)),]))

dim(sample1.dat)
dim(sample1.surv)

# Attach N=K to sample surv etc before calculating the lrs!!!
# Attached eqm sizes to the sample data and work on them!!!

sample1.survK <- addK_func(sample1.surv)
sample1.survK <- as.data.frame(sample1.survK)

sample1.recK <- addK_func(as.data.frame(sample1.rec))
sample1.recK <- as.data.frame(sample1.recK)

sample1.groK <- addK_func(sample1.gro)
sample1.groK <- as.data.frame(sample1.groK)

sample1.surv.groK <- addK_func(sample1.surv.gro)
sample1.surv.groK <- as.data.frame(sample1.surv.groK)

# Calculate LRS for ALL perturbed params
sample.lrs.all.pert <- lrs_calc_func(sample1.datK)

# Calculate LRS for perturbed params
sample.surv.lrs <- lrs_calc_func(sample1.survK)
sample.rec.lrs <- lrs_calc_func(sample1.recK)
sample.gro.lrs <- lrs_calc_func(sample1.groK)
sample.surv.gro.lrs <- lrs_calc_func(sample1.surv.groK)

# add column name for plot
sample.lrs.all.pert$den_name <- factor(sample.lrs.all.pert$den, levels = c("0", "0.5", "1"),
                                       labels = c("At N = 0", "At N = K/2", "At N = K"))

sample.surv.lrs$den_name <- factor(sample.surv.lrs$den, levels = c("0", "0.5", "1"),
                                   labels = c("At N = 0", "At N = K/2", "At N = K"))

sample.rec.lrs$den_name <- factor(sample.rec.lrs$den, levels = c("0", "0.5", "1"),
                                  labels = c("At N = 0", "At N = K/2", "At N = K"))

sample.gro.lrs$den_name <- factor(sample.gro.lrs$den, levels = c("0", "0.5", "1"),
                                  labels = c("At N = 0", "At N = K/2", "At N = K"))

sample.surv.gro.lrs$den_name <- factor(sample.surv.gro.lrs$den, levels = c("0", "0.5", "1"),
                                       labels = c("At N = 0", "At N = K/2", "At N = K"))

# Plot for all parameters together perturbed! this is the plot for it!
ggplot(data=sample.lrs.all.pert, aes(x = mids, y = density,
                                     group=as.factor(sim),
                                     color=as.factor(sim)))+
  geom_point(size=1)+
  geom_line(data=sample.lrs.all.pert, aes(x = mids, y = density,
                                          group=sim,
                                          color=as.factor(sim)),
            size=1.2)+
  labs(x="Number of offspring", y="Probability",
       title="Change Recruitment")+
  theme_bw()+
  facet_grid(~den_name)+
  theme(legend.position = "none")+
  scale_color_brewer(palette ="Set2")+
  theme(
    axis.title.y=element_text(size=20),
    plot.title = element_text(size=25),
    axis.title.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=13),
    strip.text.x = element_text(size = 20),
    strip.text = element_text(face = "bold"))+
  scale_x_continuous(breaks = seq(0, 30, by = 1), limits=c(0,30))

# plots for recruitment, survival, growth change etc
sample_rec <- ggplot(data=sample.rec.lrs, aes(x = mids, y = density,
                                              group=as.factor(sim),
                                              color=as.factor(sim)))+
  geom_point(size=1)+
  geom_line(data=sample.rec.lrs, aes(x = mids, y = density,
                                     group=sim,
                                     color=as.factor(sim)),
            size=1.2)+
  labs(x="Number of offspring", y="Probability",
       title="Change Recruitment")+
  theme_bw()+
  facet_grid(~den)+
  theme(legend.position = "none")+
  scale_color_brewer(palette ="Set2")+
  theme(
    axis.title.y=element_text(size=20),
    plot.title = element_text(size=25),
    axis.title.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=13),
    strip.text.x = element_text(size = 20),
    strip.text = element_text(face = "bold"))

sample_surv <- ggplot(data=sample.surv.lrs, aes(x = mids, y = density,
                                                group=as.factor(sim),
                                                color=as.factor(sim)))+
  geom_point(size=1)+
  geom_line(data=sample.surv.lrs, aes(x = mids, y = density,
                                      group=sim,
                                      color=as.factor(sim)),
            size=1.2)+
  labs(x="Number of offspring", y="Probability", title="Change in Survival")+
  theme_bw()+
  theme(legend.position = "none")+
  facet_grid(~den)+
  scale_color_brewer(palette ="Set2")+
  theme(
    axis.title.y=element_text(size=20),
    plot.title = element_text(size=25),
    axis.title.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=13),
    strip.text.x = element_text(size = 20),
    strip.text = element_text(face = "bold"))

sample_gro <- ggplot(data=sample.gro.lrs, aes(x = mids, y = density,
                                              group=as.factor(sim),
                                              color=as.factor(sim)))+
  geom_point(size=1)+
  geom_line(data=sample.gro.lrs, aes(x = mids, y = density,
                                     group=sim,
                                     color=as.factor(sim)),
            size=1.2)+
  labs(x="Number of offspring", y="Probability",
       title="Change Growth")+
  theme_bw()+
  facet_grid(~den)+
  theme(legend.position = "none")+
  scale_color_brewer(palette ="Set2")+
  theme(
    axis.title.y=element_text(size=20),
    plot.title = element_text(size=25),
    axis.title.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=13),
    strip.text.x = element_text(size = 20),
    strip.text = element_text(face = "bold"))


sample_surv_gro <- ggplot(data=sample.surv.gro.lrs, aes(x = mids, y = density,
                                                        group=as.factor(sim),
                                                        color=as.factor(sim)))+
  geom_point(size=1)+
  geom_line(data=sample.surv.gro.lrs, aes(x = mids, y = density,
                                          group=sim,
                                          color=as.factor(sim)),
            size=1.2)+
  labs(x="Number of offspring", y="Probability",
       title="Change Survival and Growth")+
  theme_bw()+
  facet_grid(~den)+
  theme(legend.position = "none")+
  scale_color_brewer(palette ="Set2")+
  theme(
    axis.title.y=element_text(size=20),
    plot.title = element_text(size=25),
    axis.title.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=13),
    strip.text.x = element_text(size = 20),
    strip.text = element_text(face = "bold"))

# Compute average for lrs distributions
detach(package:plyr)
detach(package:ggbiplot)

# These are the commands for Figure 5 in the main manuscript,
# and Figure A12, Figure A13 in the main manuscript

test_all <- sample.lrs.all.pert %>% group_by(mids, den, den_name) %>%
  dplyr::summarise(mean_den= mean(density),
                   sd_den= sd(density),
                   n_den=n()) %>%
  mutate(se_den = sd_den / sqrt(n_den),
         lower.ci_den = mean_den - qt(1 - (0.05 / 2), n_den - 1) * se_den,
         upper.ci_den = mean_den + qt(1 - (0.05 / 2), n_den - 1) * se_den)

test_s <- sample.surv.lrs %>% group_by(mids, den, den_name) %>%
  dplyr::summarise(mean_den= mean(density),
                   sd_den= sd(density),
                   n_den=n()) %>%
  mutate(se_den = sd_den / sqrt(n_den),
         lower.ci_den = mean_den - qt(1 - (0.05 / 2), n_den - 1) * se_den,
         upper.ci_den = mean_den + qt(1 - (0.05 / 2), n_den - 1) * se_den)

test_r <- sample.rec.lrs %>% group_by(mids, den, den_name) %>%
  dplyr::summarise(mean_den= mean(density),
                   sd_den= sd(density),
                   n_den=n()) %>%
  mutate(se_den = sd_den / sqrt(n_den),
         lower.ci_den = mean_den - qt(1 - (0.05 / 2), n_den - 1) * se_den,
         upper.ci_den = mean_den + qt(1 - (0.05 / 2), n_den - 1) * se_den)

test_g <- sample.gro.lrs %>% group_by(mids, den, den_name) %>%
  dplyr::summarise(mean_den= mean(density),
                   sd_den= sd(density),
                   n_den=n()) %>%
  mutate(se_den = sd_den / sqrt(n_den),
         lower.ci_den = mean_den - qt(1 - (0.05 / 2), n_den - 1) * se_den,
         upper.ci_den = mean_den + qt(1 - (0.05 / 2), n_den - 1) * se_den)

test_sg <- sample.surv.gro.lrs %>% group_by(mids, den, den_name) %>%
  dplyr::summarise(mean_den= mean(density),
                   sd_den= sd(density),
                   n_den=n()) %>%
  dplyr::mutate(se_den = sd_den / sqrt(n_den),
                lower.ci_den = mean_den - qt(1 - (0.05 / 2), n_den - 1) * se_den,
                upper.ci_den = mean_den + qt(1 - (0.05 / 2), n_den - 1) * se_den)

plot_test_all <- ggplot(data=test_all, aes(x = mids, y = mean_den))+
  geom_point(size=1)+
  geom_line(data=test_all, aes(x = mids, y = mean_den),
            size=1)+
  geom_errorbar(aes(ymin=mean_den-sd_den, ymax=mean_den+sd_den), width=.5,
                position=position_dodge(0.05), color="deeppink3", size=0.5)+
  labs(x="Number of offspring", y="Probability",
       title="Perturb all parameters")+
  theme_bw()+
  facet_grid(~den_name)+
  theme(legend.position = "none")+
  scale_color_brewer(palette ="Set2")+
  theme(
    axis.title.y=element_text(size=20),
    plot.title = element_text(size=25),
    axis.title.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=13),
    strip.text.x = element_text(size = 20),
    strip.text = element_text(face = "bold"))+
  scale_x_continuous(breaks = seq(0, 20, by = 2), limits=c(0,20))

plot_test_all


plot_test_s <- ggplot(data=test_s, aes(x = mids, y = mean_den))+
  geom_point(size=1)+
  geom_line(data=test_s, aes(x = mids, y = mean_den),
            size=1)+
  geom_errorbar(aes(ymin=mean_den-sd_den, ymax=mean_den+sd_den), width=.5,
                position=position_dodge(0.05), color="deeppink3", size=0.5)+
  labs(x="Number of offspring", y="Probability",
       title="Change Survival")+
  theme_bw()+
  facet_grid(~den)+
  theme(legend.position = "none")+
  scale_color_brewer(palette ="Set2")+
  theme(
    axis.title.y=element_text(size=20),
    plot.title = element_text(size=25),
    axis.title.x=element_text(size=20),

    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=13),
    strip.text.x = element_text(size = 20),
    strip.text = element_text(face = "bold"))+
  scale_x_continuous(breaks = seq(0, 20, by = 2), limits=c(0,20))
plot_test_s


plot_test_g <- ggplot(data=test_g, aes(x = mids, y = mean_den))+
  geom_point(size=1)+
  geom_line(data=test_g, aes(x = mids, y = mean_den),
            size=1)+
  geom_errorbar(aes(ymin=mean_den-sd_den, ymax=mean_den+sd_den), width=.5,
                position=position_dodge(0.05), color="deeppink3", size=0.5)+
  labs(x="Number of offspring", y="Probability",
       title="Change Growth")+
  theme_bw()+
  facet_grid(~den)+
  theme(legend.position = "none")+
  scale_color_brewer(palette ="Set2")+
  theme(
    axis.title.y=element_text(size=20),
    plot.title = element_text(size=25),
    axis.title.x=element_text(size=20),
    # axis.title.x = element_blank(),
    # axis.title.y = element_blank(),

    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=13),
    strip.text.x = element_text(size = 20),
    strip.text = element_text(face = "bold"))+
  scale_x_continuous(breaks = seq(0, 16, by = 1), limits=c(0,16))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits=c(0,1))+
  scale_x_continuous(breaks = seq(0, 20, by = 2), limits=c(0,20))

plot_test_g

plot_test_r <- ggplot(data=test_r, aes(x = mids, y = mean_den))+
  geom_point(size=1)+
  geom_line(data=test_r, aes(x = mids, y = mean_den),
            size=1)+
  geom_errorbar(aes(ymin=mean_den-sd_den, ymax=mean_den+sd_den), width=.5,
                position=position_dodge(0.05), color="deeppink3", size=0.5)+
  labs(x="Number of offspring", y="Probability",
       title="Change Recruitment")+
  theme_bw()+
  # scale_colour_gradient(low = "darkgreen", high = "darkred")
  # scale_colour_viridis_d(option = "inferno")+
  facet_grid(~den)+
  theme(legend.position = "none")+
  scale_color_brewer(palette ="Set2")+
  theme(
    axis.title.y=element_text(size=20),
    plot.title = element_text(size=25),
    axis.title.x=element_text(size=20),
    # axis.title.x = element_blank(),
    # axis.title.y = element_blank(),

    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=13),
    strip.text.x = element_text(size = 20),
    strip.text = element_text(face = "bold"))+
  scale_x_continuous(breaks = seq(0, 20, by = 2), limits=c(0,20))

plot_test_r

plot_test_sg <- ggplot(data=test_sg, aes(x = mids, y = mean_den))+
  geom_point(size=1)+
  geom_line(data=test_sg, aes(x = mids, y = mean_den),
            size=1)+
  geom_errorbar(aes(ymin=mean_den-sd_den, ymax=mean_den+sd_den), width=.5,
                position=position_dodge(0.05), color="deeppink3", size=0.5)+
  labs(x="Number of offspring", y="Probability",
       title="Change Survival and Growth")+
  theme_bw()+
  # scale_colour_gradient(low = "darkgreen", high = "darkred")
  # scale_colour_viridis_d(option = "inferno")+
  facet_grid(~den)+
  theme(legend.position = "none")+
  scale_color_brewer(palette ="Set2")+
  theme(
    axis.title.y=element_text(size=20),
    plot.title = element_text(size=25),
    axis.title.x=element_text(size=20),
    # axis.title.x = element_blank(),
    # axis.title.y = element_blank(),

    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=13),
    strip.text.x = element_text(size = 20),
    strip.text = element_text(face = "bold"))+
  scale_x_continuous(breaks = seq(0, 20, by = 2), limits=c(0,20))

plot_test_sg

ggarrange(plot_test_all, plot_test_s, plot_test_sg, plot_test_r, nrow=2, ncol=2)

ggarrange(plot_test_s, plot_test_g, plot_test_sg, nrow=3, ncol=1)


# Define a range of N_t values
library(ggplot2)
library(dplyr)

# Code for Figure A2 in the manuscript
# Define the parameters
s0 <- -0.2392
s1 <- 0.1775
s2 <- -0.0043

r0 <- -1.7839
r1 <- 0.1090
r2 <- -0.0032
N_t_values <- c(0, 200, 400, 600)
z <- seq(0, 40, length.out = 400)

data <- expand.grid(z = z, N_t = N_t_values) %>%
  mutate(S = 1 / (1 + exp(-(s0 + s1 * z + s2 * N_t))))

# Create a data frame for plotting
data1 <- expand.grid(z = z, N_t = N_t_values) %>%
  mutate(R = 1 / (1 + exp(-(r0 + r1 * z + r2 * N_t))))

colors <- brewer.pal("Dark2", n = 4)

# Plot survival and recruitment function at different densities
p1 <- ggplot(data, aes(x = z, y = S, color = factor(N_t), group = N_t)) +
  geom_line(size=1.3) +
  scale_color_manual(values = colors) +
  labs(x = "Body size: z", y = "Survival: S(z)", color = "N") +
  theme_bw()+
  theme(axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        axis.text.y=element_text(size=20),
        axis.text.x=element_text(size=13),
        legend.position = "none")+
  ylim(c(0,1))

p2 <- ggplot(data1, aes(x = z, y = R, color = factor(N_t), group = N_t)) +
  geom_line(size=1.3) +
  scale_color_manual(values = colors) +
  labs(x = "Body size: z", y = "Reproduction: R(z)", color = "N") +
  theme_bw()+
  theme(axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        axis.text.y=element_text(size=20),
        axis.text.x=element_text(size=13),
        # strip.text.x = element_text(size = 15),
        legend.title = element_text( size = 22),
        legend.text = element_text(size = 14))+
  ylim(c(0,1))

ggarrange(p1, p2)

