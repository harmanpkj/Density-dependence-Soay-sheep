rm(list=ls())
library(tidyverse)
library(reshape2)
library(devtools)
library(magick)
library(ggbiplot)
library(parallel)
# install.packages("psych")
library(factoextra)
library(ggpubr)
library(cowplot)
library(corrplot)
library(ggpmisc)
library(patchwork)
library(hrbrthemes)
library(Hmisc)
library(RSpectra) # for faster eigenvalues for block matrices
library('plot.matrix')
library(patchwork)

getwd()
setwd("/Users/harmanjaggi/Documents/Research/LRS/FixedEnv/FixedEnv")
devtools::load_all(".")

# Read the params from covariance matrix samples. This file contains 10000 of them.
dat_og1 <- read.csv("rand.params.csv")

# remove the serial number column
dat_og <- dat_og1[,-1]
head(dat_og)
# hist.data.frame(dat_og)

nsim = 250
set.seed(120)
# set.seed(120)
sample1.dat <- sample_n(dat_og, nsim)
sample2.dat <- sample_n(dat_og, nsim)
sample3.dat <- sample_n(dat_og, nsim)

# Compare with optimal and suboptimal strategies
topdi <- as.data.frame(t(read.csv("top5di.csv")))
rownames(topdi) <- NULL
# topdi$params <- paste0(c(1:16),"DI")
topdd <- as.data.frame(t(read.csv("top5dd.csv")))
rownames(topdd) <- NULL
# topdd$params <- paste0(c(1:16),"DD")
head(topdd)


# sample1.dat == sample2.dat
# sample1.dat== sample3.dat
# which(c(dat_og[,1])==c(test[1]))
# dat_og[9839,]
# test
# sample2.dat <- apply(dat_og, 2, sample_func)
#
# sample3.dat <- apply(dat_og, 2, sample_func)
#
# sample4.dat <- apply(dat_og, 2, sample_func)

cv_sam1 <- sqrt(diag(var(sample1.dat)))/abs(colMeans(sample1.dat))
cv_sam2 <- sqrt(diag(var(sample2.dat)))/abs(colMeans(sample2.dat))
cv_sam3 <- sqrt(diag(var(sample3.dat)))/abs(colMeans(sample3.dat))

# make data frames
cvdat <- tibble(
  # cv_rose= cv_rose,
  cv_sam1=cv_sam1,
  cv_sam2=cv_sam2,
  cv_sam3=cv_sam3
)

# plots
ggplot(cvdat)+
  # geom_point(aes(1:16, cv_rose), color="red")+
  geom_line(aes(1:16, cv_sam2), color="red")+
  geom_line(aes(1:16, cv_sam3), color="blue")+
  geom_point(aes(1:16, cv_sam1), color="black")+
  geom_line(aes(1:16, cv_sam1), color="black")+
  theme_bw()

minsize <- 0.1  ### make this smaller than observed value
maxsize <- 38   ### make this larger than observed value
nn <- n <- 50
max.age <- 16

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
D.ratio.fun <- function(z, intercept.mu, z.slope.mu,
                  n.slope.mu,year.eff.mu,id.eff.mu,intercept.va,
                  z.slope.va,n.slope.va,year.eff.va,id.eff.va, N) {
  junk <- data.frame()
  mu.z <- c()
  ratio <- c()
  for (i in 1:length(z)) {
    mu.z[i] <- intercept.mu + z.slope.mu*z[i] + n.slope.mu*N + year.eff.mu + id.eff.mu # Eg(z)
    ratio[i] <- (mu.z[i])/(z[i])
  }
return(ratio)  
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

# n <-nn
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
  d.ratio.vec <- D.ratio.fun(y, d.params[1], d.params[2], d.params[3], d.params[4], d.params[5], d.params[6], d.params[7], d.params[8], d.params[9], d.params[10], N)
  

  
  # scale D and G so columns sum to 1
  G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  D <- D/matrix(as.vector(apply(D,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  return(list(S=S, R=R, R_kappa=R_kappa, G=G, D=D, svec=svec, rvec=rvec, 
              g.inc=g.inc.vec, d.ratio=d.ratio.vec,
              g.prop.ratio.check1=g.prop.ratio.check1, #Ez
              g.prop.ratio.check2=g.prop.ratio.check2, #z
              g.prop.ratio.check3=g.prop.ratio.check3, #ratio
              meshpts=y))
}

# n <- nn

# use truncated age-stage model
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



# toprms <- as.matrix(sample1.dat[1,])
# toprms <- as.matrix(topdd[1,])
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
    {M <- bigmatrix_agestage_lrs(nn,s.params,r.params,g.params,d.params,N)
    # matF <- matrix(NA, ncol=ncol(M$R), nrow=ncol(M$R))
    #
    # # create the recruitment matrix (here only stage 1 being born into)
    # matF[1,] <- M$R[2,]
    # add.zero <- matrix(0, ncol = ncol(M$R), nrow=ncol(M$R)-1)
    # matF[2:ncol(M$R),] <- add.zero
    
    # mat <- M$G %*% M$S + M$D %*% matF
    # mat <- M$G %*% M$S + matF
    
    mat <- M$G %*% M$S + M$D %*% M$R
    mat <- as.matrix(mat)
    
    lam <- Re(eigen(mat)$values[1])
    # matU <-  M$G %*% M$S
    # rec.mat <- M$D %*% M$R
    # print(lambda)
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
    # junk <- data.frame(old = N, new = (N+1), diff=diff, lambda=lambda)
    # junk1 <- rbind(junk1, junk)
    # print(junk1)
    # df <- rbind(df, junk1)
  }
  
  N <- N
  return(N)
  # return(junk1)
}

# For large matrices it is faster to calculate the dominant eigenvalue and eigenvectors associated with it via iteration rather than using eigen.
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

# Function to make block matrix
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


# Function to add newK to data frame
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

# Function to add newK to data frame
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

# Sample dataset with K attached
sample1.datK <- as.data.frame(addK_func(sample1.dat))

# read the file
# sample1.datK <- read.csv("sample1.datK.csv")

topdd.datK <- as.data.frame(addK_func(topdd))
topdi.datK <- as.data.frame(addK_func(topdi))

# hist(sample1.datK$newK)
# hist(topdd.datK$newK)
# hist(topdi.datK$newK)

# hist.data.frame(sample1.datK)
# hist.data.frame(sample2.datK)

# Plot for equilirbrium sizes distribution
ggplot(sample1.datK, aes(x=newK))+
  geom_histogram(fill="gray", color="black", bins=15)+
  theme_bw()+
  labs(x="Distribution of K for different life-histories")+
  scale_x_continuous(limits=c(400,500))+
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
    strip.text.x = element_text(size=14, face="bold"))

mean(sample1.datK$newK)

dat <-sample1.datK

# func returns both increment, pratio, ratio!!
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
  # names(tbl) <- NULL
  # colnames(tbl) <- c(
  #   "Survival_small_ssd", "Survival_big_ssd",
  #   "Reproduction_ssd",
  #   "Growth_Inc_small_ssd", "Growth_Inc_big_ssd", "PopN")
  
  colnames(tbl) <- c(
    "Ss", "Sb",
    "R",
    "Gsi", "Gbi", 
    "Gspr", "Gbpr", 
    "Drat",
    "PopN")
  return(tbl)
}

# Vector for eqm ratio eqm val
eqm_val <- c(0, 1/2, 1, 1.2, 1.4)

# dat <- sample1.datK
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


# Make the vitals data frame!
vitals_tbl1 <- get_named_tbl(sample1.datK, eqm_val)

vitals_tbl1$den_name <- factor(vitals_tbl1$eqratio, levels = c("0", "0.5", "1"),
                               labels = c("At N = 0", "At N = K/2", "At N = K"))
vitals_tbl1$den <- vitals_tbl1$eqratio


# regression plot to highlight growth survival tradeoff at high density

surv_growth_tradeoff <- ggplot(vitals_tbl1, aes(y=Ss, x=Gspr)) +
  
  geom_point( alpha=1, size=1.5, color="darksalmon") +
  labs(color="")+
  facet_wrap(~eqratio)+
  # geom_smooth(method='lm', formula= y~x, se=F,
  #             size=0.5, linetype="dashed", color="black")+
  # stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
  #              parse = TRUE, color = "darkblue",
  #              rr.digits = 2, coef.digits = 2, size = 4) +
  
  # # color="blue")+
  # scale_colour_manual(values = c( "darkolivegreen3", "orange2"))+
  theme_bw()+
  labs(y="Average Juvenile Survival at SSD",
       x="Average Growth Increment for juveniles at SSD")+
  # scale_y_continuous(limits = c(0, 11), breaks = seq(0, 11, by = 2))+
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

surv_growth_tradeoff

vitals_tbl1_named <- vitals_tbl1
names(vitals_tbl1_named)
temp_n0 <- as.data.frame(vitals_tbl1_named %>% filter(eqratio==0))
temp_n0 <- temp_n0[,c(1:3,6)]
row.names(temp_n0) <- NULL

# tbl.pca_dep <- prcomp(temp, scale. = TRUE, center = T)


pca_res_n0 <- prcomp(temp_n0, scale. = TRUE)
# pca_res_n0 <- prcomp(temp_n0)
pca_plot_n0 <- autoplot(pca_res_n0, data = temp_n0, loadings.label = TRUE,
                        loadings = TRUE, loadings.colour = 'black',
                        loadings.label.size = 7,
                        loadings.label.vjust = 1,
                        loadings.label.hjust = 1,
                        loadings.label.colour="red")+
  theme_minimal()+
  labs(title = "At N=0")+
  theme(
    plot.title = element_text(size=20, face="bold"))+
  geom_point(color = "orange2",  size = 1)

temp_k2 <- as.data.frame( vitals_tbl1_named %>% filter(eqratio==0.5))
temp_k2 <- temp_k2[,c(1:3,6)]
row.names(temp_k2) <- NULL

# tbl.pca_dep <- prcomp(temp, scale. = TRUE, center = T)

pca_res_k2 <- prcomp(temp_k2, scale. = TRUE)
pca_plot_k2 <-autoplot(pca_res_k2, data = temp_k2, loadings.label = TRUE,
                       loadings = TRUE, loadings.colour = 'black',
                       loadings.label.size = 7,
                       loadings.label.vjust = 1.4,
                       loadings.label.hjust = 0.1,
                       loadings.label.colour="red")+
  theme_minimal()+
  labs(title = "At N=K/2")+
  theme(
    plot.title = element_text(size=20, face="bold"))+
  geom_point(color = "orange2",  size = 1)

temp_k <- as.data.frame( vitals_tbl1_named %>% filter(eqratio==1))
temp_k <- temp_k[,c(1:3,6)]
row.names(temp_k) <- NULL

# tbl.pca_dep <- prcomp(temp, scale. = TRUE, center = T)

pca_res_k <- prcomp(temp_k, scale. = TRUE)
pca_plot_k <-autoplot(pca_res_k, data = temp_k,
                      loadings.colour = 'black',
                      # color=NULL,
                      loadings = TRUE,
                      loadings.label = TRUE,
                      geom = "arrow",
                      loadings.label.size = 7,
                      loadings.label.vjust = 1,
                      loadings.label.hjust = 1.4,
                      loadings.label.colour="red"
)+
  theme_minimal()+
  labs(title = "At N=K")+
  theme(
    plot.title = element_text(size=20, face="bold"))+
  geom_point(color = "orange2",  size = 1)

ggarrange(pca_plot_n0, pca_plot_k2, pca_plot_k, nrow=3)



x <- vitals_tbl1$Ss
logit_Ss <- log(x/(1-x))
names(vitals_tbl1)
ggplot(vitals_tbl1, aes(y=(logit_Ss), x=(log(Gspr)))) +
  
  geom_point( alpha=1, size=1.5, color="darksalmon") +
  labs(color="")+
  facet_wrap(~eqratio)+
  geom_smooth(method='lm', formula= y~x, se=F,
              size=0.5, linetype="dashed", color="black")+
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
               parse = TRUE, color = "darkblue",
               rr.digits = 2, coef.digits = 2, size = 4) +
  
  # # color="blue")+
  # scale_colour_manual(values = c( "darkolivegreen3", "orange2"))+
  theme_bw()+
  labs(x="Average Juvenile Survival at SSD",
       y="Average Growth Increment for juveniles at SSD")+
  # scale_y_continuous(limits = c(0, 11), breaks = seq(0, 11, by = 2))+
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

# vitals_tbldd <- get_named_tbl(topdd.datK, eqm_val)
# vitals_tbldi <- get_named_tbl(topdi.datK, eqm_val)


# check sigularity of (I-U)
sing_check <- function(m) class(try(solve(m),silent=T))=="matrix"

age.at.death <- 16
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
    
    # plot(cumprod(diag(M$S)))
    # par(mar=c(4,4,4,4))
    
    mat <- as.matrix(mat)
    lambda <- Re(eigen(mat)$values[1])
    
    diff_lambda = lambda-lambda_block
    diff_lambda
    sad1 <- Mod(eigen(mat)$vector[,1])
    sad <- sad1/sum(sad1)
    # sum(sad)
    # M$D %*% matF[1,]*sad
    
    newF <-  M$D %*% M$R
    matrix_list[[i]] <- newF
    
    # newF <-  M$D %*% diag(M$R) *(sad)
    newborns_dist1 <- colSums(t(newF))
    
    # normalize the newborn distribution
    newborns_dist <- newborns_dist1/sum(newborns_dist1)
    # plot(newF)
    
    # Use 'plot.matrix' as library and then plot the
    # plot(M$D %*% matF *sad)
    # plot(newF)
    # plot(newborns_dist)
    # unload('plot.matrix')
    # library('plot.matrix')
    # plot(newF)
    # plot(newborns_dist)
    # # distribution of mothers in population
    # plot(diag(M$R)*sad) # which is same as
    # #
    # plot(diag(M$R)*sad, pch=20)
    # lines(diag(M$R)*sad)
    
    # plot(M$D%*%matF[1,]*sad)
    # plot(M$D%*%diag(M$R)*sad)
    
    # plot(as.vector(M$D %*% matF[1,]*sad))
    # plot(M$D %*%M$R*sad)
    # matplot(diag(M$R)*sad, type = c("b"), pch=20)
    # matplot(newF, type = c("b"), pch=20,col = 1:50)
    
    
    
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
    # colSums(gamma)
    # colSums(new_gamma)
    # plot(new_gamma)
    # plot(gamma[,1])
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
    # View(U_block)
    # dim(N_block)
    R0_block_mat <- F_block %*% N_block
    # dim(R0_block)
    R0_block <- Re(eigen(R0_block_mat)$values[1])
    # R0_block
    
    #     # Get R0 and Tc by age and stage.
    #     age.at.death <- 100
    #     phi <- c()
    #     L <- I
    #     ns <- nrow(mat)
    #     B1.sum <- diag(0, ns, ns)
    #     Tc.sum <- 0
    #     n.age <- c(1:age.at.death)
    #
    #     for (k in 1:length(n.age)) {
    #       junk <- M$D %*% M$R %*% L
    #       junk <- junk[1,1]
    #       phi[k] <- junk
    #       L <- matU %*% L
    #       # B1.junk <- n.age[k] * matF %*% L
    #       # B1.sum <- B1.sum + B1.junk
    #       # Tc.junk <- n.age[k] * phi[k]
    #       # # print(Tc.junk)
    #       # Tc.sum <- Tc.sum + Tc.junk
    #     }
    #
    #     # check singularity here
    #     if(sing_check(I-matU)!= TRUE) next
    #     N1 <- solve(I - matU) %*% (I-matU^(length(n.age))) ## geometric finite sum formula
    #
    #     R01 <- M$D %*% M$R %*% N1
    #     R01 <- R01[1,1]
    #     # R01 <- matF %*% N1
    #     # R01[1,1]
    #
    #     phi <- phi
    #     Tc <- weighted.mean(n.age, phi)
    #
    # Get Generation Time using function above
    Tc <- gentime(M)
    
    temp$sim <- rep(i, nrow(temp))
    temp$lambda <- rep(lambda, nrow(temp))
    temp$mean.lrs.dist <- rep(junk1$R0, nrow(temp))
    temp$Tc <- rep(Tc, nrow(temp))
    temp$R0_eig <- rep(R0_eig, nrow(temp))
    temp$R0_block <- rep(R0_block, nrow(temp))
    lrs_gamma.distr <- rbind(lrs_gamma.distr, temp)
    
    # ggplot(temp, aes(mids, density))+
    #    geom_line()
    
  }
  
  return(lrs_gamma.distr)
}

dat <- sample1.dat
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
    # Re(eigen(block_mat)$values[1])
    
    # plot(cumprod(diag(M$S)))
    # par(mar=c(4,4,4,4))
    
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
    # colSums(gamma)
    # colSums(new_gamma)
    # plot(new_gamma)
    # plot(gamma[,1])
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

unique(lrs_gamma.distr_N$N_val)
lrs_gamma.distr_N$R0_eig

ggplot(lrs_gamma.distr_N, aes(x=N_val, y=R0_block)) +
  
  geom_point(alpha=0.65, size=1.8) +
  labs(y="R0", x="Population abundance: N")+
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
    strip.text.x = element_text(size=14, face="bold"))

# ggplot(lrs_gamma.distr_N, aes(x=N_val, y=R0_eig)) +
#   
#   geom_point(alpha=0.65, size=1.8) +
#   labs(y="R0", x="Population abundance: N")+
#   theme_bw()+
#   # labs(x="Juvenile Survival scaled by SSD")+
#   theme(
#     # axis.title.y=element_blank(),
#     axis.title.x=element_text(size=22),
#     axis.title.y=element_text(size=22),
#     axis.text.y=element_text(size=20),
#     axis.text.x=element_text(size=13),
#     # strip.text.x = element_text(size = 15),
#     # legend.title = element_text( size = 22),
#     legend.text = element_text(size = 22),
#     strip.text.x = element_text(size=14, face="bold"))

detach(package:ggbiplot)
detach(package:plyr)

# lrs_gamma.distr_N <- read.csv("lrs_gamma.distr_N.csv")
mean_var_R0_at_N <- lrs_gamma.distr_N %>% group_by(N_val) %>%
  mutate(mean_mean_lrs=mean(mean.lrs.dist),var_mean_lrs=var(mean.lrs.dist))

unique(mean_var_R0_at_N$var_mean_lrs)
unique(mean_var_R0_at_N$var_mean_lrs/mean_var_R0_at_N$mean_mean_lrs)

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
write.csv(lrs_gamma.distr_N, "lrs_gamma.distr_N.csv")


# RESCALE AND GET EXPECTATION
# dat <- temp
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

# dat <- as.data.frame(sample1.datK)

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

# Function for mother's distribution!
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
  # eqm_ratio <- 0
  # sad_temp <- data.frame(sno=1:50)
  # i=5
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
    
    # plot(cumprod(diag(M$S)))
    # par(mar=c(4,4,4,4))
    
    mat <- as.matrix(mat)
    lambda <- Re(eigen(mat)$values[1])
    
    diff_lambda = lambda-lambda_block
    diff_lambda
    sad1 <- Mod(eigen(mat)$vector[,1])
    sad <- sad1/sum(sad1)
    # sum(sad)
    # M$D %*% matF[1,]*sad
    
    newF <-  M$D %*% (M$R*sad)
    matrix_list[[i]] <- newF
    
  }
  # Calculate the average matrix
  average_Fmat <- Reduce(`+`, matrix_list) / nsim
  
  return(average_Fmat)
}

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


matF_return <- matF_plot_func(sample1.datK)

# Combine matrices into one data frame
combined_df <- do.call(rbind, lapply(matF_return, function(mat) {
  melt(mat)
}))

# Overall range of values for colors
# overall_range <- range(combined_df$value)
overall_range <- c(0, max(combined_df$value))

# Create a color palette
custom_palette <- c("gray90", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7",
                    "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061")


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

mother_plots

combined_mother_plot <- ggarrange(mother_plots[[1]],
                                  mother_plots[[2]],
                                  mother_plots[[3]], nrow=1, common.legend = T, legend="right")
annotate_figure(combined_mother_plot, left = textGrob("Distribution of Offspring size", rot = 90, vjust = 1, gp = gpar(cex = 1.4)),
                bottom = textGrob("Distribution of Mothers size (scaled by SSD)", vjust=-6, gp = gpar(cex = 1.4)))

combined_mother_plot <- mother_plots[[1]] +
  mother_plots[[2]] +
  mother_plots[[3]] & theme(legend.position = "right")
combined_mother_plot + plot_layout(guides = "collect")

ggplot(dat1_vitals_meanlrs, aes(x=PopN, y=mean_lrs)) +
  
  geom_point(alpha=0.65, size=1.8) +
  facet_wrap(~den_name)+
  labs(y="R0", x="Population abundance: N")+
  
  # geom_smooth(method='lm', formula= y~x, se=F,
  #             size=1,
  #             # linetype="dashed",
  #             aes(group = variable, color=variable))+
  # color="blue")+
  # scale_colour_manual(values = c(
  #   # "cyan4",
  #   #                              "darkolivegreen3",
  #   #
  #   "orange2"))+
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
  scale_y_continuous(breaks = seq(0, 100, by = 2))


eqm_val <- c(0, 0.5, 1)
sam1_meanlrs <- meanlrs_vital_func(sample1.datK)
sam1_meanlrs$den_name <- factor(sam1_meanlrs$den, levels = c("0", "0.5", "1"),
                                labels = c("At N = 0", "At N = K/2", "At N = K"))

nrow(sam1_meanlrs)
write.csv(sam1_meanlrs, "sam1_meanlrs.csv")
#Merge the vitals and meanlrs datasets
# dat1_vitals_meanlrs <- cbind(vitals_tbl1, filter(sam1_meanlrs, mean_lrsnot0!=0))
dat1_vitals_meanlrs <- cbind(vitals_tbl1, sam1_meanlrs)
dat1_vitals_meanlrs <- data.frame(dat1_vitals_meanlrs)

# sam1_meanlrs <- read.csv("sam1_meanlrs.csv")

# Check R0's from distribution and from stage structured matrix.
# Check R0 from block matrix and compare with stage only matrix.
names(dat1_vitals_meanlrs)
dat1_vitals_meanlrs$R0_block
R0_K_plot <- ggplot(dat1_vitals_meanlrs, aes(x=PopN, y=mean_lrs)) +
  
  geom_point(alpha=0.65, size=1.8) +
  facet_wrap(~den_name)+
  labs(y="R0", x="Population abundance: N")+
  
  # geom_smooth(method='lm', formula= y~x, se=F,
  #             size=1,
  #             # linetype="dashed",
  #             aes(group = variable, color=variable))+
  # color="blue")+
  # scale_colour_manual(values = c(
  #   # "cyan4",
  #   #                              "darkolivegreen3",
  #   #
  #   "orange2"))+
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
  scale_y_continuous(breaks = seq(0, 100, by = 2))

R0_K_plot

# par(mfrow=c(1,1))
R0_hist <- ggplot(sam1_meanlrs,
                  aes(x=mean_lrs))+
  geom_histogram(fill="cyan4", color="#e9ecef", alpha=0.8, bins = 15)+
  facet_wrap(~den_name)+
  theme_bw()+
  labs(x="Histogram for R0: mean lrs")+
  theme(axis.title.y=element_text(size=22),
        axis.title.x=element_text(size=22),
        axis.text.y=element_text(size=20),
        axis.text.x=element_text(size=20, angle = 90),
        strip.text.x = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 22))
# scale_x_continuous(breaks= seq(0,30,5), limits = c(0,30))

R0_hist

x <- filter(dat1_vitals_meanlrs, den==1)$PopN
y1 <- filter(dat1_vitals_meanlrs, den==0)$mean_lrs
y2 <- filter(dat1_vitals_meanlrs, den==0)$R0_block
y3 <- filter(dat1_vitals_meanlrs, den==0)$R0_eig

z <- filter(sam1_meanlrs, den==0)$Tc


# R0 result at density independence
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

ggplot(tibble(z,x), aes(z,x))+
  geom_point(alpha=0.9, size=2.5) +
  theme_bw()+
  # labs(x="Carrying Capacity", y="R0 at density-independence (N=0)")+
  # stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
  #              parse = TRUE, color = "darkblue",
  #              rr.digits = 2, coef.digits = 2, size = 5)+
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


ggplot(tibble(x,y2), aes(x,y2))+
  geom_point(alpha=0.9, size=2.5) +
  theme_bw()+
  labs(x="Carrying Capacity", y="R0 at density-independence (N=0)")+
  # stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
  #              parse = TRUE, color = "darkblue",
  #              rr.digits = 2, coef.digits = 2, size = 5)+
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

# Tc plot at density independence
ggplot(tibble(x,z), aes(x,z))+
  geom_point(alpha=0.9, size=2.5) +
  theme_bw()+
  labs(x="Carrying Capacity", y="Tc at density-independence (N=0)")+
  # stat_poly_eq(
  # formula = y~x,
  # formula = y~poly(x, 2, raw = TRUE),
  # aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
  # parse = TRUE, color = "darkblue",
  # rr.digits = 2, coef.digits = 2, size = 5)+
  
  # stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
  #              parse = TRUE, color = "darkblue",
  #              rr.digits = 2, coef.digits = 2, size = 5)+
  # geom_smooth(method='lm', formula= y~x, se=F,
#             size=0.3, linetype="dashed", color="red" )+
stat_smooth(method = "lm", formula = y ~ x, size = 0.2, color="red", se=F)+
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
# scale_y_continuous(breaks = seq(0, 15, by = 2.5))

# Life exp at density independence
# ggplot(tibble(x,z1), aes(x,z1))+
#   geom_point(alpha=0.9, size=2.5) +
#   theme_bw()+
#   labs(x="Carrying Capacity", y="Life Exp at density-independence (N=0)")+
#   stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
#                parse = TRUE, color = "darkblue",
#                rr.digits = 2, coef.digits = 2, size = 5)+
#   geom_smooth(method='lm', formula= y~x, se=F,
#               size=0.3, linetype="dashed", color="red" )+
#   theme(axis.title.y=element_text(size=18),
#         axis.title.x=element_text(size=18),
#         axis.text.y=element_text(size=20),
#         axis.text.x=element_text(size=13),
#         # strip.text.x = element_text(size = 15),
#         # legend.title = element_text( size = 22),
#         legend.text = element_text(size = 22),
#         strip.text.x = element_text(size=15, face="bold"),
#         legend.position = "bottom")+
#   guides(fill=guide_legend(title="New Legend Title"))
# scale_y_continuous(breaks = seq(0, 15, by = 2.5))

ggplot(tibble(y2,z,x), aes(x, log(y2)/z))+
  geom_point(alpha=0.9, size=2.5) +
  theme_bw()+
  labs(x="Carrying Capacity", y="log(R0)/Tc at density-independence (N=0)")+
  # stat_poly_eq(
  #             # formula = y~x,
  #             formula = y~poly(x, 2, raw = TRUE),
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
  #              parse = TRUE, color = "darkblue",
  #              rr.digits = 2, coef.digits = 2, size = 5)+
  
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
               parse = TRUE, color = "darkblue",
               rr.digits = 2, coef.digits = 2, size = 5)+
  # geom_smooth(method='lm', formula= y~x, se=F,
  #             size=0.3, linetype="dashed", color="red" )+
  stat_smooth(method = "lm", formula = y ~ x , size = 0.2, color="red")+
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


ggplot(tibble(y2,z,x), aes(y2, log(y2)/z))+
  geom_point(alpha=0.9, size=2.5) +
  theme_bw()+
  labs(x="R0 at density-independence (N=0)", y="log(R0)/Tc at density-independence (N=0)")+
  stat_poly_eq(
    formula = y~x,
    # formula = y~poly(x, 2, raw = TRUE),
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
    parse = TRUE, color = "darkblue",
    rr.digits = 2, coef.digits = 2, size = 5)+
  geom_smooth(method='lm', formula= y~x, se=F,
              size=0.3, linetype="dashed", color="red" )+
  # stat_smooth(method = "lm", formula = y ~ x , size = 0.2, color="red")+
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


# scale_y_continuous(breaks = seq(0, 15, by = 2.5))

# SSD at different densities
# toprms <- as.data.frame(as.matrix(sample1.datK[1,]))
# function for sad/ssd

# toprms <- sample1.datK[1,]
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
    
    # create the recruitment matrix (here only stage 1 being born into)
    # matF[1,] <- M$R[2,]
    # add.zero <- matrix(0, ncol = ncol(M$R), nrow=ncol(M$R)-1)
    # matF[2:ncol(M$R),] <- add.zero
    
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
    
    # rv <- as.vector(goat[[3]]/goat[[3]][1]) # Reproductive value
  }
  
  sad_tbl <- as.data.frame(sad_tbl)
  # sad_tbl$N <- c(rep(0,50), rep(0.5,50), rep(1,50))
  
  return(sad_tbl)
}

sad_res1 <- apply(sample1.datK, 1, sad_func)
sad_res <- as.data.frame(sad_res1)
dim(sad_res)

colnames(sad_res) <- c(1:100)

df <- data.frame()
df1<- data.frame(Sno=1:nrow(sad_res))
dim(sad_res)

for (i in 1:(length(sad_res)-2))
{
  df1$sad <-  sad_res[,i]
  df1 <- as.data.frame(df1)
  df1$sim <- rep(i, nrow(df1))
  df1$size <- rep(seq(1, 50, 1),3)
  df1$eqm <- c(rep(0, 50), rep(0.5, 50), rep(1, 50))
  df <- rbind(df, df1)
}

df
names(df)
detach(package:ggbiplot)
detach(package:plyr)
library(dplyr)

df %>% group_by(eqm, size) %>% summarise(mean_sad=mean(sad)) %>%
  ggplot(., aes(x=size , y=mean_sad, group=as.factor(eqm), color=as.factor(eqm)))+
  geom_line(size=1.6)+
  geom_point(size=1.2)+
  # facet_wrap(.~eqm)+
  # xlim(c(2,50))+
  theme_bw()+
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"))+
  labs(y="Stable Size Distribution",
       x="Number of Stages",
       color="N/K")+
  theme(axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        axis.text.y=element_text(size=20),
        axis.text.x=element_text(size=13),
        # strip.text.x = element_text(size = 15),
        legend.title = element_text( size = 18),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size=15, face="bold"))


write.csv(vitals_tbl1, "vitals_tbl1.csv")

# Make the correlation plots
# corrdat1_N0 <- filter(dat1_vitals_meanlrs, den==0)[,1:5]
# colnames(corrdat1_N0) <- c("Ss", "Sb", "R", "GIs", "GIb")
#
# corrdat1_K_2 <- filter(dat1_vitals_meanlrs, den==0.5)[,1:5]
# colnames(corrdat1_K_2) <- c("Ss", "Sb", "R", "GIs", "GIb")
#
# corrdat1_K <- filter(dat1_vitals_meanlrs, den==1)[,1:5]
# colnames(corrdat1_K) <- c("Ss", "Sb", "R", "GIs", "GIb")

corrdat1_N0 <- filter(vitals_tbl1, den==0)[,c(1:3,6,7)]
colnames(corrdat1_N0) <- c("Ss", "Sb", "R", "GIs", "D")

corrdat1_K_2 <- filter(vitals_tbl1, den==0.5)[,c(1:3,6,7)]
colnames(corrdat1_K_2) <- c("Ss", "Sb", "R", "GIs", "D")

corrdat1_K <- filter(vitals_tbl1, den==1)[,c(1:3,6,7)]
colnames(corrdat1_K) <- c("Ss", "Sb", "R", "GIs", "D")

cor_mat_sam1_N0 <- cor(corrdat1_N0)
cor_mat_sam1_K_2 <- cor(corrdat1_K_2)
cor_mat_sam1_K <- cor(corrdat1_K)


pdf("corplot_sam1_N0.png")
corrplot(cor_mat_sam1_N0, method="color", addCoef.col = "black")
dev.off()

pdf("corplot_sam1_K_2.png")
corrplot(cor_mat_sam1_K_2, method="color", addCoef.col = "black")
dev.off()

pdf("corplot_sam1_K.png")
corrplot(cor_mat_sam1_K, method="color", addCoef.col = "black")
dev.off()

pdf("corplot_sam1_new.pdf")
par(mfrow=c(3,1))
corrplot(cor_mat_sam1_N0, method="color", addCoef.col = "black", number.cex=0.8)
corrplot(cor_mat_sam1_K_2, method="color", addCoef.col = "black", number.cex=0.8)
corrplot(cor_mat_sam1_K, method="color", addCoef.col = "black", number.cex=0.8)
dev.off()

pdf("corplot_sam1_ratio.pdf")
par(mfrow=c(3,1))
corrplot(cor_mat_sam1_N0, method="color", addCoef.col = "black", number.cex=0.8)
corrplot(cor_mat_sam1_K_2, method="color", addCoef.col = "black", number.cex=0.8)
corrplot(cor_mat_sam1_K, method="color", addCoef.col = "black", number.cex=0.8)
dev.off()

pdf("corplot_sam1_beyondK.pdf")
par(mfrow=c(1,2))
corrplot(cor_mat_sam1_2K1, method="color", addCoef.col = "black", number.cex=0.8)
# corrplot(cor_mat_sam1_2K2, method="color", addCoef.col = "black", number.cex=0.8)
corrplot(cor_mat_sam1_2K3, method="color", addCoef.col = "black", number.cex=0.8)
# corrplot(cor_mat_sam1_2K4, method="color", addCoef.col = "black", number.cex=0.8)
dev.off()



# Check some basic plots
sam1_meanlrs %>% pivot_longer(cols = -c("plrs0",
                                        "plrs_not0",  "den", "den_name", "lambda",
                                        "R0_eig", "R0_dist", "var_lrsnot0"),
                              names_to = "variable", values_to = "value") %>%
  ggplot(aes(x=1:600, y=value)) +
  geom_point(aes(color=variable), alpha=0.8, size=3.2) +
  labs(color="")+
  facet_wrap(~den)+
  geom_smooth(method='lm', formula= y~x, se=F,
              size=0.5, linetype="dashed",
              # aes(group = variable, color=variable)
              color="blue")+
  scale_colour_manual(values = c("cyan4", "darkolivegreen3"))+
  theme_bw()+
  labs(x="Average Juvenile Survival at SSD")+
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
               parse = TRUE, color = "darkblue",
               rr.digits = 2, coef.digits = 2, size = 3) +
  theme(
    # axis.title.y=element_blank(),
    axis.title.y=element_text(size=22),
    axis.title.x=element_text(size=22),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=13),
    # strip.text.x = element_text(size = 15),
    legend.title = element_text( size = 22),
    legend.text = element_text(size = 22),
    strip.text.x = element_blank(),
    legend.position = "bottom")


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
a


b <- ggplot(data.frame(dat1_vitals_meanlrs), aes(x=Ss, y=mean_lrsnot0)) +
  
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

# Repeat and plot with Recruitment function

c <- ggplot(dat1_vitals_meanlrs, aes(R, plrs0))+
  geom_point(size=2.5, alpha=0.9)+
  # geom_point(dat1_vitals_meanlrs, mapping=aes(Survival_small_ssd, mean_lrs), color="red", size=2.5)+
  # scale_y_continuous(sec.axis = sec_axis(~./10, name="R0")) +
  facet_wrap(~den_name)+
  theme_bw()+
  labs(y="Probability that LRS is 0",
       x="Average reproduction at ssd")+
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
               parse = TRUE, color = "blue",
               rr.digits = 2, coef.digits = 2, size = 3) +
  geom_smooth(method='lm', formula= y~x, se=F, color="blue", size=0.5, linetype="dashed")+
  theme(axis.title.y=element_text(size=20),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.text.x=element_text(size=13),
        strip.text.x = element_text(size = 20),
        strip.text = element_text(face = "bold"))
# scale_y_continuous(breaks = seq(0, 1.1, by = 0.10))

c

# print(a/c)
# aaa <- dat1_vitals_meanlrs %>% filter(den==0)
# aaa$mean_lrsnot0

d <- dat1_vitals_meanlrs %>% pivot_longer(cols = -c("Ss", "Sb", "Reproduction_ssd",
                                                    "Growth_Inc_small_ssd", "Growth_Inc_big_ssd", "eqratio",  "plrs0",
                                                    "plrs_not0",  "den",  "den_name", "var_lrsnot0", "den_name",
                                                    "PopN", "R0_eig", "lambda", "R0_dist", "den.1", "den_name.1",
                                                    "Tc", "R0_block", "mean_lrs"),
                                          names_to = "variable", values_to = "value") %>%
  ggplot(aes(x=R, y=value)) +
  labs(color="")+
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
               parse = TRUE, color = "blue",
               rr.digits = 2, coef.digits = 2, size = 3) +
  
  facet_wrap(~den)+
  geom_point(color="darkolivegreen3", alpha=0.8, size=3.2) +
  
  geom_smooth(method='lm', formula= y~x, se=F,
              size=0.5, linetype="dashed",
              color="black"
  )+
  # scale_colour_manual(values = c("deeppink3", "orange2", "palegreen"))+
  theme_bw()+
  labs(x="Average Recruitment at SSD", y=expression(gamma))+
  theme(axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.text.x=element_text(size=13),
        # strip.text.x = element_text(size = 15),
        # legend.title = element_text( size = 22),
        legend.text = element_text(size = 22),
        strip.text.x = element_blank(),
        legend.position="bottom", legend.box = "horizontal")


e <- ggplot(data.frame(dat1_vitals_meanlrs), aes(x=Growth_Inc_small_ssd, y=mean_lrsnot0)) +
  
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

e
print(b/d)
print(c/d)
print(a/b)
print(e/f)
plot_R0_plrs0_comb <- (a/b)



#Vital rate functions and mean lrs
# Add function
add_func_sample <- function(a, b)
{
  cc1 <- c()
  for(i in 1:ncol(a))
  {
    cc <- as.numeric(b[,i]) + as.numeric(a[,i])
    cc1 <- cbind(cc1, cc)
  }
  return(cc1)
}

# Mean from 10000 params (one way to do it!)
mean.prm <- c()
head(dat_og)
dim(dat_og)

for(i in 1:ncol(dat_og))
{mean.prm[i] <- mean(dat_og[,i])
}
mean.prm <- as.matrix(mean.prm)

# mean.prm <- as.data.frame(mean.prm)
nsim =100

mean.prm_K <- eqm_func(as.matrix(mean.prm))
dat <- as.matrix(c(as.matrix(mean.prm), as.matrix(mean.prm_K)))
# seed 120
eqm_val <- c(0,0.5,1)

# Matrix of mean params where each col is sims and each row is mean of params
mean.prm.mat <- as.data.frame(matrix(as.matrix(mean.prm), nrow=nrow(mean.prm), ncol=nsim))
dim(mean.prm.mat)

# Function for survival, recruitment, growth and combined changes
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

# Structure works! not a problem
sample1.test <- t(rbind(t(sample1.dat)[1:3,],
                        t(sample1.dat)[4:6,],
                        t(sample1.dat)[7:9,],
                        t(sample1.dat)[10:nrow(t(sample1.dat)),]))

dim(sample1.dat)
dim(sample1.surv)
typeof(sample1.dat)
typeof(sample1.surv)
typeof(as.data.frame(as.data.frame(sample1.surv)))
sample1.surv
#Attach N=K to sample surv etc before calculating the lrs!!! Attached eqm sizes to the sample data and work on them!!!




sample1.survK <- addK_func(sample1.surv)
sample1.survK <- as.data.frame(sample1.survK)

sample1.recK <- addK_func(as.data.frame(sample1.rec))
sample1.recK <- as.data.frame(sample1.recK)

sample1.groK <- addK_func(sample1.gro)
sample1.groK <- as.data.frame(sample1.groK)

sample1.surv.groK <- addK_func(sample1.surv.gro)
sample1.surv.groK <- as.data.frame(sample1.surv.groK)

# sample1.testK <- addK_func(sample1.test)

# Write dates with their eqm values
write.csv(sample1.survK, "sample1.survK")
write.csv(sample1.recK, "sample1.recK")
write.csv(sample1.groK, "sample1.groK")
write.csv(sample1.surv.groK, "sample1.surv.groK")
write.csv(sample1.datK, "sample1.datK")

write.csv(sample1.survK, "sample1.survK.csv")
write.csv(sample1.recK, "sample1.recK.csv")
write.csv(sample1.groK, "sample1.groK.csv")
write.csv(sample1.surv.groK, "sample1.surv.groK.csv")
write.csv(sample1.datK, "sample1.datK.csv")


# sample1.survK <- read.csv("sample1.survK.csv")

# sample1.survK <- read.csv(file = "sample1.survK")[,2:18]

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

# Save the damn results and dont run again and again!!!
write.csv(sample.lrs.all.pert, "sample.lrs.all.pert.csv")
write.csv(sample.surv.lrs, "sample.surv.lrs.csv")
write.csv(sample.rec.lrs, "sample.rec.lrs.csv")
write.csv(sample.gro.lrs, "sample.gro.lrs.csv")
write.csv(sample.surv.gro.lrs, "sample.surv.gro.lrs.csv")

sample.surv.lrs <- read.csv("sample.surv.lrs.csv")
sample.rec.lrs <- read.csv("sample.rec.lrs.csv")
sample.gro.lrs <- read.csv("sample.gro.lrs.csv")
sample.surv.gro.lrs <- read.csv("sample.surv.gro.lrs.csv")


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
  # scale_colour_gradient(low = "darkgreen", high = "darkred")
  # scale_colour_viridis_d(option = "inferno")+
  facet_grid(~den_name)+
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
  
  # geom_line(data=filter(sample.lrs, model=="101"), aes(x = mids, y = density), color="black", size=0.5)+
  # labs(x = "Number of offspring (LRS)",
  #      y = expression(log[10](Probablity)))+
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
    strip.text = element_text(face = "bold"))
# scale_x_continuous(breaks = seq(0, 30, by = 1), limits=c(0,30))


sample_rec

# ggsave("sample_rec_zoom.png", sample_rec_zoom)

sample_surv <- ggplot(data=sample.surv.lrs, aes(x = mids, y = density,
                                                group=as.factor(sim),
                                                color=as.factor(sim)))+
  geom_point(size=1)+
  geom_line(data=sample.surv.lrs, aes(x = mids, y = density,
                                      group=sim,
                                      color=as.factor(sim)),
            size=1.2)+
  
  # geom_line(data=filter(sample.lrs, model=="101"), aes(x = mids, y = density), color="black", size=0.5)+
  # labs(x = "Number of offspring (LRS)",
  #      y = expression(log[10](Probablity)))+
  labs(x="Number of offspring", y="Probability", title="Change in Survival")+
  theme_bw()+
  theme(legend.position = "none")+
  facet_grid(~den)+
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
    strip.text = element_text(face = "bold"))
# scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits=c(0,1))

sample_surv


sample_gro <- ggplot(data=sample.gro.lrs, aes(x = mids, y = density,
                                              group=as.factor(sim),
                                              color=as.factor(sim)))+
  geom_point(size=1)+
  geom_line(data=sample.gro.lrs, aes(x = mids, y = density,
                                     group=sim,
                                     color=as.factor(sim)),
            size=1.2)+
  
  # geom_line(data=filter(sample.lrs, model=="101"), aes(x = mids, y = density), color="black", size=0.5)+
  # labs(x = "Number of offspring (LRS)",
  #      y = expression(log[10](Probablity)))+
  labs(x="Number of offspring", y="Probability",
       title="Change Growth")+
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
    strip.text = element_text(face = "bold"))
# scale_x_continuous(breaks = seq(0, 16, by = 1), limits=c(0,16))
# scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits=c(0,1))


sample_gro


# check for subset of 7
# unique(sample.gro.lrs.7$den)

sample_surv_gro <- ggplot(data=sample.surv.gro.lrs, aes(x = mids, y = density,
                                                        group=as.factor(sim),
                                                        color=as.factor(sim)))+
  geom_point(size=1)+
  geom_line(data=sample.surv.gro.lrs, aes(x = mids, y = density,
                                          group=sim,
                                          color=as.factor(sim)),
            size=1.2)+
  
  # geom_line(data=filter(sample.lrs, model=="101"), aes(x = mids, y = density), color="black", size=0.5)+
  # labs(x = "Number of offspring (LRS)",
  #      y = expression(log[10](Probablity)))+
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
    strip.text = element_text(face = "bold"))
# scale_x_continuous(breaks = seq(0, 16, by = 1), limits=c(0,16))
# scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits=c(0,1))

sample_surv_gro
# Plot error and standard deviations

# Compute average for lrs distributions
detach(package:plyr)
detach(package:ggbiplot)

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
  # scale_colour_gradient(low = "darkgreen", high = "darkred")
  # scale_colour_viridis_d(option = "inferno")+
  facet_grid(~den_name)+
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
# scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits=c(0,1))

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
# scale_x_continuous(breaks = seq(0, 16, by = 1), limits=c(0,16))+
# scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits=c(0,1))
plot_test_s

sample_surv

plot_test_g <- ggplot(data=test_g, aes(x = mids, y = mean_den))+
  geom_point(size=1)+
  geom_line(data=test_g, aes(x = mids, y = mean_den),
            size=1)+
  geom_errorbar(aes(ymin=mean_den-sd_den, ymax=mean_den+sd_den), width=.5,
                position=position_dodge(0.05), color="deeppink3", size=0.5)+
  labs(x="Number of offspring", y="Probability",
       title="Change Growth")+
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

sample_surv
xx1 <- test_s %>% filter(den==1)
xx2 <- test_g %>% filter(den==1)
xx3 <- test_sg %>% filter(den==1)

xxs_g <- xx1$mean_den+xx1$mean_den
xxsg <- xx3$mean_den

plot(xxsg, xxs_g)
lines(xxsg, xxs_g)
abline(a=0, b=1, col="blue")

dd=1
trade_off_mean_func <- function(df_s, df_g, df_sg, dd)
{
  s <- df_s %>% filter(den==dd)
  g <- df_g %>% filter(den==dd)
  sg <- df_sg %>% filter(den==dd)
  add_sg <- s$mean_den + g$mean_den
  sg <- sg$mean_den
  df <- tibble(add_sg, sg)
  
  diff <- c()
  
  for(i in 1:length(sg))
  {
    diff[i] <-  sg[i]- add_sg[i]
  }
  
  # plot(diff, pch=20)
  # lines(diff)
  # abline(0, 0, col="red")
  
  ggplot(as.data.frame(diff), aes(1:length(diff), diff))+
    geom_point()+
    geom_line()+
    theme_bw()+
    geom_hline(yintercept = 0, color="red")+
    ylim(c(-1,0))+
    labs(title=paste0("At den: ", dd))
  # plot!!
  # ggplot(df, aes(add_sg, sg))+
  #   geom_point()+
  #   geom_line()+
  #   geom_segment(x=0, y=0, color="red", xend=2, yend=2, linetype="dashed")+
  #   # geom_abline(slope=1, color="red")+
  #   theme_bw()+
  #   ylim(c(-0.4, 2))+
  #   xlim(c(-0.4, 2))+
  #   labs(title=paste0("At den: ", dd), x= "Add Survival and Growth", y= "Vary Survival and Growth together")
}

dd=1
df_s <- test_s
df_g <- test_g
df_sg <- test_sg

trade_off_var_func <- function(df_s, df_g, df_sg, dd)
{
  s <- df_s %>% filter(den==dd)
  g <- df_g %>% filter(den==dd)
  sg <- df_sg %>% filter(den==dd)
  add_var_sg <- (s$sd_den) + (g$sd_den)
  var_sg <- (sg$sd_den)
  df_var <- tibble(add_var_sg, var_sg)
  
  # add_var_sg <- (s$sd_den) + (g$sd_den)
  # var_sg <- (sg$sd_den)
  # df_var <- tibble(add_var_sg, var_sg)
  
  diff_var <- c()
  # print(length(s))
  
  for(i in 1:length(add_var_sg))
  {
    diff_var[i] <-  var_sg[i]- add_var_sg[i]
  }
  print(length(diff_var))
  
  # plot(diff, pch=20)
  # lines(diff)
  # abline(0, 0, col="red")
  
  ggplot(as.data.frame(diff_var), aes(0:(length(diff_var)-1), diff_var))+
    geom_point()+
    geom_line()+
    theme_bw()+
    geom_hline(yintercept = 0, color="red")+
    ylim(c(-0.015,0))+
    labs(title=paste0("At den: ", dd))+
    scale_x_continuous(breaks = seq(0, 15, by = 1), limits=c(0,3))
  # scale_y_continuous(breaks = seq(-0.0014, 0, by = 0.0005), limits=c(-0.0014,0))
  # plot!!
  # ggplot(df_var, aes(add_var_sg, var_sg))+
  #   geom_point()+
  #   geom_line()+
  #   geom_segment(x=-0.00001, y=-0.00001, color="red", xend=0.00001, yend=0.00001, linetype="dashed")+
  #   # geom_abline(slope=1, color="red")+
  #   theme_bw()+
  #   ylim(c(-0.0001, 0.00001))+
  #   xlim(c(-0.0001, 0.00001))+
  #   labs(title=paste0("At den: ", dd), x= "Add Survival and Growth", y= "Vary Survival and Growth together")
}


trade_off_cv_func <- function(df_s, df_g, df_sg, dd)
{
  s <- df_s %>% filter(den==dd)
  g <- df_g %>% filter(den==dd)
  sg <- df_sg %>% filter(den==dd)
  # add_var_sg <- (s$se_den) + (g$se_den)
  # var_sg <- (sg$se_den)
  # df_var <- tibble(add_var_sg, var_sg)
  
  add_var_sg <- (s$sd_den)/s$mean_den + (g$sd_den)/g$mean_den
  var_sg <- (sg$sd_den)/sg$mean_den
  df_var <- tibble(add_var_sg, var_sg)
  
  add_sg <- s$mean_den + g$mean_den
  sg <- sg$mean_den
  df <- tibble(add_sg, sg)
  
  diff_var <- c()
  diff <- c()
  # print(length(s))
  
  for(i in 1:length(add_var_sg))
  {
    diff_var[i] <-  var_sg[i]- add_var_sg[i]
    diff[i] <-  sg[i]- add_sg[i]
  }
  # print(length(diff_var))
  # dd
  # plot(diff, pch=20)
  # lines(diff)
  # abline(0, 0, col="red")
  # print(diff_var)
  ggplot(as.data.frame(diff_var), aes(0:(length(diff_var)-1), diff_var))+
    geom_point()+
    geom_line()+
    geom_line(df, mapping=aes(0:(length(diff_var)-1), diff), color="blue")+
    geom_point(df, mapping=aes(0:(length(diff)-1), diff), color="blue")+
    theme_bw()+
    geom_hline(yintercept = 0, color="red")+
    ylim(c(-0.5,0))+
    labs(title=paste0("At den: ", dd))+
    xlim(c(0,3))
  # +
  #   ylim(c(-0.001,0.001))
  # scale_x_continuous(breaks = seq(0, 15, by = 1), limits=c(0,15))+
  # scale_y_continuous(breaks = seq(-0.0001, 0.001, by = 0.0005), limits=c(-0.0001,0.001))
  # plot!!
  # ggplot(df_var, aes(add_var_sg, var_sg))+
  #   geom_point()+
  #   geom_line()+
  #   geom_segment(x=-0.00001, y=-0.00001, color="red", xend=0.00001, yend=0.00001, linetype="dashed")+
  #   # geom_abline(slope=1, color="red")+
  #   theme_bw()+
  #   ylim(c(-0.0001, 0.00001))+
  #   xlim(c(-0.0001, 0.00001))+
  #   labs(title=paste0("At den: ", dd), x= "Add Survival and Growth", y= "Vary Survival and Growth together")
}

trade_off_cv_func1 <- function(df_s, df_g, df_sg, dd)
{
  s <- df_s %>% filter(den==dd)
  g <- df_g %>% filter(den==dd)
  sg <- df_sg %>% filter(den==dd)
  # add_var_sg <- (s$se_den) + (g$se_den)
  # var_sg <- (sg$se_den)
  # df_var <- tibble(add_var_sg, var_sg)
  
  add_var_sg <- (s$sd_den)/s$mean_den + (g$sd_den)/g$mean_den
  var_sg <- (sg$sd_den)/sg$mean_den
  df_var <- tibble(add_var_sg, var_sg)
  
  add_sg <- s$mean_den + g$mean_den
  sg <- sg$mean_den
  df <- tibble(add_sg, sg)
  
  diff_var <- c()
  diff <- c()
  # print(length(s))
  
  for(i in 1:length(add_var_sg))
  {
    diff_var[i] <-  var_sg[i]- add_var_sg[i]
    diff[i] <-  sg[i]- add_sg[i]
  }
  
  return(diff_var)
  # print(length(diff_var))
  # dd
  # plot(diff, pch=20)
  # lines(diff)
  # abline(0, 0, col="red")
  # print(diff_var)
  # ggplot(as.data.frame(diff_var), aes(0:(length(diff_var)-1), diff_var))+
  #   geom_point()+
  #   geom_line()+
  #   geom_line(df, mapping=aes(0:(length(diff_var)-1), diff), color="blue")+
  #   geom_point(df, mapping=aes(0:(length(diff)-1), diff), color="blue")+
  #   theme_bw()+
  #   geom_hline(yintercept = 0, color="red")+
  #   ylim(c(-0.5,0))+
  #   labs(title=paste0("At den: ", dd))+
  #   xlim(c(0,15))
  # +
  #   ylim(c(-0.001,0.001))
  # scale_x_continuous(breaks = seq(0, 15, by = 1), limits=c(0,15))+
  # scale_y_continuous(breaks = seq(-0.0001, 0.001, by = 0.0005), limits=c(-0.0001,0.001))
  # plot!!
  # ggplot(df_var, aes(add_var_sg, var_sg))+
  #   geom_point()+
  #   geom_line()+
  #   geom_segment(x=-0.00001, y=-0.00001, color="red", xend=0.00001, yend=0.00001, linetype="dashed")+
  #   # geom_abline(slope=1, color="red")+
  #   theme_bw()+
  #   ylim(c(-0.0001, 0.00001))+
  #   xlim(c(-0.0001, 0.00001))+
  #   labs(title=paste0("At den: ", dd), x= "Add Survival and Growth", y= "Vary Survival and Growth together")
}


diff1 <- trade_off_cv_func1(test_s, test_g, test_sg, 0)
diff2 <- trade_off_cv_func1(test_s, test_g, test_sg, 0.5)
diff3 <- trade_off_cv_func1(test_s, test_g, test_sg, 1)

diff1 <- as.data.frame(diff1)
diff1$den <- rep(0, nrow(diff1))
colnames(diff1) <- c("diff", "den")

diff2 <- as.data.frame(diff2)
diff2$den <- rep(0.5, nrow(diff2))
colnames(diff2) <- c("diff", "den")


diff3 <- as.data.frame(diff3)
diff3$den <- rep(1, nrow(diff3))
colnames(diff3) <- c("diff", "den")


diff_df <- rbind(diff1, diff2, diff3)
ggplot(diff_df, aes(0:89, diff, color=as.factor(den), group=as.factor(den)))+
  geom_point()+
  theme_bw()

ggplot()+
  geom_point(diff1, mapping=aes(0:29, diff), color="red")+
  geom_line(diff1, mapping=aes(0:29, diff), color="red")+
  geom_point(diff2, mapping=aes(0:29, diff), color="blue")+
  geom_line(diff2, mapping=aes(0:29, diff), color="blue")+
  geom_point(diff3, mapping=aes(0:29, diff), color="black")+
  geom_line(diff3, mapping=aes(0:29, diff), color="black")+
  theme_bw()+
  xlim(c(0,10))+
  ylim(-0.8,0)

p11 <- trade_off_func(test_s, test_g, test_sg, 0)
p21 <- trade_off_func(test_s, test_g, test_sg, 0.5)
p31 <- trade_off_func(test_s, test_g, test_sg, 1)

pvar11 <- trade_off_var_func(test_s, test_g, test_sg, 0)
pvar21 <- trade_off_var_func(test_s, test_g, test_sg, 0.5)
pvar31 <- trade_off_var_func(test_s, test_g, test_sg, 1)

pvar111 <- trade_off_cv_func(test_s, test_g, test_sg, 0)
pvar211 <- trade_off_cv_func(test_s, test_g, test_sg, 0.5)
pvar311 <- trade_off_cv_func(test_s, test_g, test_sg, 1)

pvar1 <- trade_off_var_func(test_s, test_g, test_sg, 0)
pvar2 <- trade_off_var_func(test_s, test_g, test_sg, 0.5)
pvar3 <- trade_off_var_func(test_s, test_g, test_sg, 1)

pvar111
pvar211
pvar311

p1 <- trade_off_func(test_s, test_g, test_sg, 0)
p2 <- trade_off_func(test_s, test_g, test_sg, 0.5)
p3 <- trade_off_func(test_s, test_g, test_sg, 1)

sg_trade_plot <- ggarrange(p1, p2, p3, nrow=3, ncol=1)
sg_trade_plot

sg_trade_neg_plot <- ggarrange(p11 + rremove("ylab") + rremove("xlab"),
                               p21 + rremove("ylab") + rremove("xlab"),
                               p31+ rremove("ylab") + rremove("xlab"), nrow=3, ncol=1)
annotate_figure(sg_trade_neg_plot,
                left = textGrob("Varied together (SG) - Added (S+G)", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("No. of offspring", rot = 0, vjust = 0, gp = gpar(cex = 1.3)))
sg_trade_neg_plot

sg_trade_var_plot <- ggarrange(pvar11 + rremove("ylab") + rremove("xlab"),
                               pvar21 + rremove("ylab") + rremove("xlab"),
                               pvar31+ rremove("ylab") + rremove("xlab"), nrow=3, ncol=1)
annotate_figure(sg_trade_var_plot,
                left = textGrob("Varied together (SG) - Added (S+G)", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("No. of offspring", rot = 0, vjust = 0, gp = gpar(cex = 1.3)))

sg_trade_cv_plot <- ggarrange(pvar111 + rremove("ylab") + rremove("xlab"),
                              pvar211 + rremove("ylab") + rremove("xlab"),
                              pvar311+ rremove("ylab") + rremove("xlab"), nrow=3, ncol=1)
annotate_figure(sg_trade_cv_plot,
                left = textGrob("Varied together (SG) - Added (S+G)", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("No. of offspring", rot = 0, vjust = 0, gp = gpar(cex = 1.3)))

sg_trade_var_plot


sg_trade_var_plot1 <- ggarrange(pvar1 + rremove("ylab") + rremove("xlab"),
                                pvar2 + rremove("ylab") + rremove("xlab"),
                                pvar3+ rremove("ylab") + rremove("xlab"), nrow=3, ncol=1)
annotate_figure(sg_trade_var_plot1,
                left = textGrob("Varied together (SG) - Added (S+G)", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("No. of offspring", rot = 0, vjust = 0, gp = gpar(cex = 1.3)))
names(sam_lrs1)

sub_s <- filter(test_s, den==0.5)
sub_g <- filter(test_g, den==0.5)
sub_sg <- filter(test_sg, den==0.5)
tibble(y=sub_s$se_den-sub_sg$se_den)
plot(sub_s$se_den-sub_sg$se_den)

plot_test_r


sample_rec_zoom <- ggplot(data=sample.rec.lrs, aes(x = mids, y = density,
                                                   group=as.factor(sim),
                                                   color=as.factor(sim)))+
  geom_point(size=1)+
  geom_line(data=sample.rec.lrs, aes(x = mids, y = density,
                                     group=sim,
                                     color=as.factor(sim)),
            size=1.2)+
  
  # geom_line(data=filter(sample.lrs, model=="101"), aes(x = mids, y = density), color="black", size=0.5)+
  # labs(x = "Number of offspring (LRS)",
  #      y = expression(log[10](Probablity)))+
  labs(x="Number of offspring", y="Probability",
       title="zoomed")+
  theme_bw()+
  # scale_colour_gradient(low = "darkgreen", high = "darkred")
  # scale_colour_viridis_d(option = "inferno")+
  theme(legend.position = "none")+
  scale_color_brewer(palette ="Set2")+
  facet_grid(~den)+
  theme(
    # axis.title.y=element_text(size=20),
    plot.title = element_text(size=25),
    # axis.title.x=element_text(size=20),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=13),
    strip.text.x = element_text(size = 20),
    strip.text = element_text(face = "bold"))+
  scale_x_continuous(breaks = seq(0, 16, by = 1), limits=c(0,2))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits=c(0,1))

# scale_x_continuous(breaks = seq(0, 3, by = 1), limits=c(0,2))+
sample_surv_zoom <- ggplot(data=sample.surv.lrs, aes(x = mids, y = density,
                                                     group=as.factor(sim),
                                                     color=as.factor(sim)))+
  geom_point(size=1)+
  geom_line(data=sample.surv.lrs, aes(x = mids, y = density,
                                      group=sim,
                                      color=as.factor(sim)),
            size=1.2)+
  
  # geom_line(data=filter(sample.lrs, model=="101"), aes(x = mids, y = density), color="black", size=0.5)+
  # labs(x = "Number of offspring (LRS)",
  #      y = expression(log[10](Probablity)))+
  labs(x="Number of offspring", y="Probability", title="zoomed")+
  theme_bw()+
  # scale_colour_gradient(low = "darkgreen", high = "darkred")
  # scale_colour_viridis_d(option = "inferno")+
  theme(legend.position = "none")+
  facet_grid(~den)+
  scale_color_brewer(palette ="Set2")+
  theme(
    plot.title = element_text(size=25),
    # axis.title.x=element_text(size=20),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=13),
    strip.text.x = element_text(size = 20),
    strip.text = element_text(face = "bold"))+
  facet_grid(~den_name)+
  scale_x_continuous(breaks = seq(0, 16, by = 1), limits=c(0,2))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits=c(0,1))

sample_surv_zoom

sample_surv_gro_zoom <- ggplot(data=sample.surv.gro.lrs, aes(x = mids, y = density,
                                                             group=as.factor(sim),
                                                             color=as.factor(sim)))+
  geom_point(size=1)+
  geom_line(data=sample.surv.gro.lrs, aes(x = mids, y = density,
                                          group=sim,
                                          color=as.factor(sim)),
            size=1.2)+
  
  # geom_line(data=filter(sample.lrs, model=="101"), aes(x = mids, y = density), color="black", size=0.5)+
  # labs(x = "Number of offspring (LRS)",
  #      y = expression(log[10](Probablity)))+
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
  scale_x_continuous(breaks = seq(0, 16, by = 1), limits=c(0,2))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits=c(0,1))

sample_surv_zoom <- ggplot(data=sample.surv.lrs, aes(x = mids, y = density,
                                                     group=as.factor(sim),
                                                     color=as.factor(sim)))+
  geom_point(size=1)+
  geom_line(data=sample.surv.gro.lrs, aes(x = mids, y = density,
                                          group=sim,
                                          color=as.factor(sim)),
            size=1.2)+
  geom_line(data=sample.surv.lrs, aes(x = mids, y = density,
                                      group=sim,
                                      color=as.factor(sim)),
            size=1.2)+
  labs(x="Number of offspring", y="Probability",
       title="Change Survival")+
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
  scale_x_continuous(breaks = seq(0, 16, by = 1), limits=c(0,2))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits=c(0,1))

sample_surv
sample_surv_gro
sample_rec
sample_gro

sample_rec_zoom
sample_surv_gro_zoom
sample_surv_zoom

# sample_surv_gro
# sample_rec_zoom
# sample_surv_zoom


# Calculate KL divergence
# function for KL divergence
# this is divergence of x from y. x is numerator!!
kl_div <- function(x, y) {
  sum(x * log(x / y), na.rm = TRUE)
}


# sample.surv.gro.lrs.7
sub1 <- sample.surv.lrs %>% filter(sim==1)


dat <- sub1
# dat1 <- sample.surv.lrs %>%
# dat2 <- sample.surv.lrs %>% filter(sim==3)

kl_eqm <- data.frame()
# eqm_val <- c(0, 1)
i=1
get_kl_func <- function(dat)
{
  for(i in 1:length(eqm_val))
  {
    print(i)
    kl_dat <- data.frame()
    # k=2
    for (k in 1:(length(unique(dat$sim))-1))
    {
      count = k+1
      kl_col <- vector()
      print(count)
      
      x <- dat %>% filter(sim==k) %>% filter(den==eqm_val[i]) %>% dplyr::select(density)
      
      for (j in count:length(unique(dat$sim)))
      {
        y <- dat %>% filter(sim==j) %>% filter(den==eqm_val[i]) %>% dplyr::select(density)
        kl_col[j] <- kl_div(x,y)
      }
      kl_df <- as.data.frame(kl_col)
      kl_df$count <- rep(count, nrow(kl_df))
      kl_df
      kl_df$den <- eqm_val[i]
      kl_dat <- rbind(kl_dat, kl_df)
      kl_eqm <- rbind(kl_eqm, kl_dat)
      print(count)
      # colnames(kl_dat) <- "kl_val"
    }
  }
  return(kl_eqm)
}


sub1 <- sample.surv.lrs[1:500, ]
# a <- get_kl_func(sample.surv.lrs)
a <- get_kl_func(as.data.frame(sub1))
# b <- get_kl_func(sample.surv.gro.lrs)



# a$sim <- rep(6,nrow(a))
# b <- get_kl_func(sample.surv.gro.lrs)
# b$sim <- rep(3,nrow(b))
# nrow(a)
ggplot(a, aes(x=1:length(kl_col),y=kl_col))+
  geom_point()+
  theme_bw()+
  facet_wrap(~den)

View(a)
c <- rbind(a,b)
c$sno <- rep(1:3,2)

ggplot(c, mapping = aes(x=sno,y=kl_val, color=as.factor(sim), group=as.factor(sim)))+
  geom_point()+
  geom_line()+
  theme_bw()

surv_plot <- ggarrange(sample_surv, sample_surv_zoom, nrow=2)

annotate_figure(surv_plot, left = textGrob("Probability", rot = 90, vjust = 1, gp = gpar(cex = 1.7)),
                bottom = textGrob("Number of offspring", gp = gpar(cex = 1.7)))

rec_plot <- ggarrange(sample_rec, sample_rec_zoom, nrow=2)

annotate_figure(rec_plot, left = textGrob("Probability", rot = 90, vjust = 1, gp = gpar(cex = 1.7)),
                bottom = textGrob("Number of offspring", gp = gpar(cex = 1.7)))


sample_surv_zoom <- ggplot(data=sample.surv.lrs, aes(x = mids, y = density,
                                                     group=as.factor(model),
                                                     color=as.factor(model)))+
  geom_point(size=2.5)+
  geom_line(data=sample.surv.lrs, aes(x = mids, y = density,
                                      group=model,
                                      color=as.factor(model)),
            size=0.8)+
  
  # geom_line(data=filter(sample.lrs, model=="101"), aes(x = mids, y = density), color="black", size=0.5)+
  # labs(x = "Number of offspring (LRS)",
  #      y = expression(log[10](Probablity)))+
  xlim(c(0,2))+
  ylim(c(0,max(sample.surv.lrs$density)))+
  labs(x="LRS", y="Probability", title="Change in Surv only: zoomed")+
  theme_bw()+
  # scale_colour_gradient(low = "darkgreen", high = "darkred")
  # scale_colour_viridis_d(option = "inferno")+
  theme(legend.position = "none")+
  scale_color_brewer(palette ="Dark2")
sample_surv_zoom
ggsave("sample_surv_zoom.png", sample_surv_zoom)


# Variation in SSD and Rates!!!



lrs_at_K  <- get_lrs_eqm(sample1.dat, eqm_ratio = 1)
lrs_at_N0  <- get_lrs_eqm(sample1.dat, eqm_ratio = 0)
View(lrs_at_K)
View(lrs_at_N0)


lrs_at_K_expect <- expectation_lrsnot0(lrs_at_K)
lrs_at_N0_expect <- expectation_lrsnot0(lrs_at_N0)



lrsK.temp <- lrs_at_K %>% arrange(mids)
lrsN0.temp <- lrs_at_N0 %>% arrange(mids)
lrsK.lrs0 <- lrsK.temp[1:100,]
lrsN0.lrs0 <- lrsN0.temp[1:100,]
lrsK.lrs0$lrs_not0 <- (1-lrsK.lrs0$density)
lrsN0.lrs0$lrs_not0 <- (1-lrsN0.lrs0$density)

lrsK_plotdat <- data.frame(plrs0=lrsK.lrs0$density,
                           lrs_at_K_expect)
lrsK_plotdat$den <- rep("K", nrow(lrsK_plotdat))
lrsN0_plotdat <- data.frame(plrs0=lrsN0.lrs0$density, lrs_at_N0_expect)
lrsN0_plotdat$den <- rep("0", nrow(lrsN0_plotdat))
lrs_plotdat <- rbind(lrsK_plotdat, lrsN0_plotdat)
dim(lrsK_plotdat)
dim(lrsN0_plotdat)

ggplot(lrs_plotdat, aes(plrs0, mean_lrsnot0))+
  geom_point()+
  facet_wrap(~den)+
  theme_bw()

ggplot(lrsN0_plotdat, aes(plrs0, mean_lrsnot0))+
  geom_point()


# Make a function for this rescaling everything




View(lrsN0.lrs0)

mean(lrsK.lrs0$lrs_not0)
mean(lrsN0.lrs0$lrs_not0)
var(lrsK.lrs0$lrs_not0)
var(lrsN0.lrs0$lrs_not0)
hist(lrsK.lrs0$density)
hist(lrsK.lrs0$lrs_not0)

hist(lrsN0.lrs0$density)
hist(lrsN0.lrs0$lrs_not0)



# this is want to repeat for many different N like for R0
lrs_N <- data.frame()
newborn_stage <- 10
for(i in 1:length(N))
{
  temp  <- get_lrs_eqm(sample1.dat, eqm_ratio = 1)
  temp$N <- rep(N[i], nrow(temp))
  lrs_N <- rbind(lrs_N, temp)
  print(i)
}


View(lrs_N)

# plot PCA res from vitals
par(mar = c(0, 0, 0, 0))
dat <- vitals_tbl1 %>% filter(eqratio==0)
dat <- vitals_tbl1 %>% filter(eqratio==0.5)
dat <- vitals_tbl1 %>% filter(eqratio==1)
dat <- vitals_tbl1 %>% filter(eqratio==1.2)

dat <- vitals_tbl1_ratio %>% filter(eqratio==0)
dat <- vitals_tbl1_ratio %>% filter(eqratio==0.5)
dat <- vitals_tbl1_ratio %>% filter(eqratio==1)

dat <- vitals_tbl1
dat <- vitals_tbldd
dat <- vitals_tbldi

dat <- vitals_tbldi %>% filter(eqratio==0)
dat <- vitals_tbldi %>% filter(eqratio==0.5)
dat <- vitals_tbldi %>% filter(eqratio==1)


install.packages("vegan")
library(vegan)

plot_res_func <- function(dat)
{
  temp <- as.data.frame(dat)
  temp <- temp[,1:5]
  row.names(temp) <- NULL
  # as.matrix(temp)
  tbl.pca_dep <- prcomp(temp, scale. = TRUE, center = T)
  pca_plot <- ggbiplot::ggbiplot(tbl.pca_dep,
                                 scale=0.1)+
    theme_bw()+
    # labs(subtitle=paste0("density = ", "K/4"))+
    theme(axis.title.y=element_text(size=16),
          axis.title.x=element_text(size=16),
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14))
  print(pca_plot)
  return(pca_plot)
}

library(vegan)
pca_result <- tbl.pca_dep

# Get the observed eigenvalues from PCA
observed_eigenvalues <- pca_result$sdev^2

n <- ncol(temp)

# Calculate broken stick values
broken_stick_values <- sapply(1:n, function(k) sum(1/(k:n))/n)

# Get the eigenvalues from PCA
eigenvalues <- pca_result$sdev^2


# Compare observed eigenvalues to broken stick values
plot(1:n, eigenvalues, type = 'b', pch = 19, col = 'blue', xlab = "Component", ylab = "Eigenvalue", main = "PCA Eigenvalues vs. Broken Stick")
lines(1:n, broken_stick_values, type = 'b', pch = 17, col = 'red')
legend("topright", legend = c("Observed Eigenvalues", "Broken Stick"), col = c("blue", "red"), pch = c(19, 17))

library(ggbiplot)
library(stats)
library(ggrepel)
library(broom)  # devtools::install_github("tidymodels/broom")
library(cowplot)


pca_result <- prcomp(temp, scale. = TRUE, center =T)

biplot <- ggbiplot(pca_result,
                   obs.scale = 1,
                   var.scale = 1,
                   groups = NULL,
                   ellipse = TRUE,
                   circle = TRUE) +
  labs(title = "Biplot of PCA") +
  theme_minimal()
biplot
biplot +  geom_text_repel(aes(label = colnames(pca_result$rotation)),
                          box.padding = 0.5,
                          point.padding = 0.2)


# PCA using tidyverse
library(ggfortify)


vitals_tbl1_gro <- vitals_tbl1_ratio
colnames(vitals_tbl1_ratio) <- c("Ss", "Sb", "R", "GIs", "GIb","PopN", "eqratio", "den_name", "den")
colnames(vitals_tbl1_gro) <- c("Ss", "Sb", "R", "GIs", "GIb","PopN", "eqratio", "den_name", "den")
temp_n0 <- as.data.frame( vitals_tbl1_ratio %>% filter(eqratio==0))
temp_n0 <- temp_n0[,1:5]
row.names(temp_n0) <- NULL

# tbl.pca_dep <- prcomp(temp, scale. = TRUE, center = T)

pca_res_n0 <- prcomp(temp_n0, scale. = TRUE)
pca_plot_n0 <- autoplot(pca_res_n0, data = temp_n0, loadings.label = TRUE,
                        loadings = TRUE, loadings.colour = 'black',
                        loadings.label.size = 7,
                        loadings.label.vjust = 1,
                        loadings.label.hjust = 1,
                        loadings.label.colour="red")+
  theme_minimal()+
  labs(title = "At N=0")+
  theme(
    plot.title = element_text(size=20, face="bold"))+
  geom_point(color = "orange2",  size = 1)

temp_k2 <- as.data.frame( vitals_tbl1_ratio %>% filter(eqratio==0.5))
temp_k2 <- temp_k2[,1:5]
row.names(temp_k2) <- NULL

# tbl.pca_dep <- prcomp(temp, scale. = TRUE, center = T)

pca_res_k2 <- prcomp(temp_k2, scale. = TRUE)
pca_plot_k2 <-autoplot(pca_res_k2, data = temp_k2, loadings.label = TRUE,
                       loadings = TRUE, loadings.colour = 'black',
                       loadings.label.size = 7,
                       loadings.label.vjust = 1.4,
                       loadings.label.hjust = 0.1,
                       loadings.label.colour="red")+
  theme_minimal()+
  labs(title = "At N=K/2")+
  theme(
    plot.title = element_text(size=20, face="bold"))+
  geom_point(color = "orange2",  size = 1)

temp_k <- as.data.frame( vitals_tbl1_ratio %>% filter(eqratio==1))
temp_k <- temp_k[,1:5]
row.names(temp_k) <- NULL

# tbl.pca_dep <- prcomp(temp, scale. = TRUE, center = T)

pca_res_k <- prcomp(temp_k, scale. = TRUE)
pca_plot_k <-autoplot(pca_res_k, data = temp_k,
                      loadings.colour = 'black',
                      # color=NULL,
                      loadings = TRUE,
                      loadings.label = TRUE,
                      geom = "arrow",
                      loadings.label.size = 7,
                      loadings.label.vjust = 1,
                      loadings.label.hjust = 1.4,
                      loadings.label.colour="red"
)+
  theme_minimal()+
  labs(title = "At N=K")+
  theme(
    plot.title = element_text(size=20, face="bold"))+
  geom_point(color = "orange2",  size = 1)

gro_ratio_pca_plot <- ggarrange(pca_plot_n0,
                                pca_plot_k2,
                                pca_plot_k, nrow=3)

gro_ratio_pca_plot


components <- seq_len(length(pca_res_n0$sdev))
data_for_plot <- data.frame(
  Component = components,
  Eigenvalue = pca_res_n0$sdev^2,
  BrokenStick = broken_stick_values
)

# Using ggplot to plot observed eigenvalues and broken stick values
bs_n0 <- ggplot(data_for_plot, aes(x = Component)) +
  geom_line(aes(y = Eigenvalue, color = "Observed Eigenvalues")) +
  geom_point(aes(y = Eigenvalue, color = "Observed Eigenvalues")) +
  geom_line(aes(y = BrokenStick, color = "Broken Stick"), linetype = "dashed") +
  geom_point(aes(y = BrokenStick, color = "Broken Stick"), shape = 17) +
  scale_color_manual(values = c("Observed Eigenvalues" = "blue", "Broken Stick" = "red")) +
  labs(x = "Component", y = "Eigenvalue", title = "PCA Eigenvalues vs Broken Stick") +
  theme_minimal() +
  guides(color = guide_legend(title = "Legend"))

components <- seq_len(length(pca_res_k2$sdev))
data_for_plot <- data.frame(
  Component = components,
  Eigenvalue = pca_res_k2$sdev^2,
  BrokenStick = broken_stick_values
)

# Using ggplot to plot observed eigenvalues and broken stick values
bs_k2 <- ggplot(data_for_plot, aes(x = Component)) +
  geom_line(aes(y = Eigenvalue, color = "Observed Eigenvalues")) +
  geom_point(aes(y = Eigenvalue, color = "Observed Eigenvalues")) +
  geom_line(aes(y = BrokenStick, color = "Broken Stick"), linetype = "dashed") +
  geom_point(aes(y = BrokenStick, color = "Broken Stick"), shape = 17) +
  scale_color_manual(values = c("Observed Eigenvalues" = "blue", "Broken Stick" = "red")) +
  labs(x = "Component", y = "Eigenvalue",
       # title = "PCA Eigenvalues vs. Broken Stick"
  ) +
  theme_minimal() +
  guides(color = guide_legend(title = "Legend"))

components <- seq_len(length(pca_res_k$sdev))
data_for_plot <- data.frame(
  Component = components,
  Eigenvalue = pca_res_k$sdev^2,
  BrokenStick = broken_stick_values
)

# Using ggplot to plot observed eigenvalues and broken stick values
bs_k <- ggplot(data_for_plot, aes(x = Component)) +
  geom_line(aes(y = Eigenvalue, color = "Observed Eigenvalues")) +
  geom_point(aes(y = Eigenvalue, color = "Observed Eigenvalues")) +
  geom_line(aes(y = BrokenStick, color = "Broken Stick"), linetype = "dashed") +
  geom_point(aes(y = BrokenStick, color = "Broken Stick"), shape = 17) +
  scale_color_manual(values = c("Observed Eigenvalues" = "blue", "Broken Stick" = "red")) +
  labs(x = "Component", y = "Eigenvalue",
       # title = "PCA Eigenvalues vs. Broken Stick"
  ) +
  theme_minimal() +
  guides(color = guide_legend(title = "Legend"))

ggarrange(bs_n0, bs_k2, bs_k, nrow=3)

write.csv(summary(pca_res_n0)$rotation, 'pca_res_n0.csv')
write.csv(summary(pca_res_k2)$rotation, 'pca_res_k2.csv')
write.csv(summary(pca_res_k)$rotation, 'pca_res_k.csv')

summary(pca_res_n0)$rotation
summary(pca_res_k2)$rotation
summary(pca_res_k)$rotation

# pca_res_k$layers[[2]]$aes_params$size <- 0.5
# pca_res_k$layers[[2]]$geom_params$arrow$length <- unit(6, units = "points")
# pca_res_k

pca_plot <- ggarrange(pca_plot_n0,
                      pca_plot_k2,
                      pca_plot_k, nrow=3)

pca_plot
ggsave("pca_plot.png", pca_plot)
# For each subcatagory

sample.juv <- as.data.frame(sample.dat) %>%
  mutate(k_col = k.juv)

sample.sub.ad <- as.data.frame(sample.dat) %>%
  mutate(k_col = k.sub.ad)

sample.ad <- as.data.frame(sample.dat) %>%
  mutate(k_col = k.ad)

sample.old <- as.data.frame(sample.dat) %>%
  mutate(k_col = k.old)

sample.old1 <- as.data.frame(sample.dat) %>%
  mutate(k_col = k.old1)

View(sample.juv)

View(sample.dat)
View(sample11)

# exactltre just try with given matrix

testf1<- Fmat
testu1 <-Umat
testa1 <- mat

testf2<- Fmat
testu2 <-Umat
testa2 <- mat
library(exactLTRE)
library(popbio)
# Approximate LTRE:
# contributions to the difference in lambda
cont_diff<- approximateLTRE(list(testa1, testa2), method='fixed')
cont_diff; # matrix of contributions
sum(cont_diff); # sum of approximate contributions
lamDiff(list(A1,A2)); # true difference in lambda

# contributions to the variance of lambda
cont_var<- approximateLTRE(list(A1,A2,A3), method='random')
round(cont_var,digits=5); # matrix of contributions
sum(cont_var); # sum of contributions
lamVar(list(A1,A2,A3)); # true variance in lambda

# Exact LTRE:

# contributions to the difference in lambda
cont_diff<- exactLTRE(list(testu1,testu1), method='random')
cbind(as.character(cont_diff$varying.indices.list),round(cont_diff$epsilons,digits=5))
sum(cont_diff$epsilons); # sum of contributions
lamDiff(list(A1,A2)); # true difference in lambda

# contributions to the variance of lambda
cont_var<- exactLTRE(list(A1,A2,A3), method='random')
cbind(as.character(cont_var$varying.indices.list),round(cont_var$epsilons,digits=5));
sum(cont_var$epsilons); # sum of contributions
lamVar(list(A1,A2,A3)); # true variance in lambda

# install.packages("devtools")
devtools::install_github("chrissy3815/exactLTRE")
library(exactLTRE)


# carrying capacity for each sub-catagory
k.juv <- c()
k.sub.ad <- c()
k.ad <- c()
k.old <- c()
k.old1 <- c()
for (i  in 1:sims)
{k.juv[i] <- sum(DD.nt.res[[i]][1:10])
k.sub.ad[i] <- sum(DD.nt.res[[i]][11:20])
k.ad[i] <- sum(DD.nt.res[[i]][21:30])
k.old[i] <- sum(DD.nt.res[[i]][31:40])
k.old1[i] <- sum(DD.nt.res[[i]][41:50])
}


autoplot(pca_res)

dfpca<- data.frame()

df1 <- data.frame()
df <- data.frame()
temp1 <- data.frame()

test <- t(sample11)

test <- t(sample.juv)
test <- t(sample.sub.ad)
test <- t(sample.ad)
test <- t(sample.old)

test <- t(sample.old1)


# all plots at K/4
juv_K_4 <- plot_res_func(t(vitals_tbl1))
sub.ad_K_4 <-plot_res_func(t(df_pca_K_4(t(sample.sub.ad))))
ad_K_4 <-plot_res_func(t(df_pca_K_4(t(sample.ad))))
old_K_4 <-plot_res_func(t(df_pca_K_4(t(sample.old))))
old1_K_4 <-plot_res_func(t(df_pca_K_4(t(sample.old1))))

# # best params
# di1<- melt(topdi, na.rm = FALSE, value.name = "values",
#            id = "params")  %>%
#   filter(variable=="V1")
#
# di2<- melt(topdi, na.rm = FALSE, value.name = "values",
#            id = "params")  %>%
#   filter(variable=="V2")
#
# di3<- melt(topdi, na.rm = FALSE, value.name = "values",
#            id = "params")  %>%
#   filter(variable=="V3")
#
# di4<- melt(topdi, na.rm = FALSE, value.name = "values",
#            id = "params")  %>%
#   filter(variable=="V4")
#
# di5<- melt(topdi, na.rm = FALSE, value.name = "values",
#            id = "params")  %>%
#   filter(variable=="V5")
#
#
#
# dd1<- melt(topdd, na.rm = FALSE, value.name = "values",
#            id = "params")  %>%
#   filter(variable=="V1")
#
# # worst params
# dd5<- melt(topdd, na.rm = FALSE, value.name = "values",
#            id = "params")  %>%
#   filter(variable=="V5")
#
#
# mean.dd <- as.data.frame(rowMeans(topdd[,1:5]))
# colnames(mean.dd) <- "values"
# mean.dd
#

test <- sample1.dat[1,]
# Define the parameter values for S and R functions
s0 <- test[1]   # Change as needed
s1 <- test[2]    # Change as needed
s2 <- test[3]  # Change as needed

s0
s1
s2
r0 <- test[4]   # Change as needed
r1 <- test[5]  # Change as needed
r2 <- test[6]  # Change as needed

r0
r1
r2
# 1 + exp(s0 + s1 * z + s2 * 1)

# Define a range of N_t values
N_t_values <- seq(0, 500, by = 100)  # Adjust the range and step size as needed

# Define a range of z values
minsize <- 0.1
maxsize <- 38
z_values <- seq(minsize, maxsize, by = 0.1)

# Create a dataframe to store results
# Create a dataframe to store results
results <- data.frame()


# Calculate S and R for each combination of N_t and z
for (i in 1:length(N_t_values)) {
  N_t <- N_t_values[i]
  
  for (j in 1:length(z_values)) {
    z <- z_values[j]
    
    # Calculate S and R using the provided expressions
    # S <- 1 / (1 + exp(-(s0 + s1 * z + s2 * N_t)))
    # R <- 1 / (1 + exp(-(r0 + r1 * z + r2 * N_t)))
    
    S <- exp(s0 + s1 * z + s2 * N_t) / (1 + exp(s0 + s1 * z + s2 * N_t))
    R <- exp(r0 + r1 * z + r2 * N_t) / (1 + exp(r0 + r1 * z + r2 * N_t))
    
    # Add the results to the dataframe
    results <- rbind(results, data.frame(N_t = N_t, z = z, S = S, R = R))
  }
}

colnames(results) <- c("N_t", "z", "S", "R")
# Create separate plots for S and R
plot_S <- ggplot(data = results, aes(x = z, y = S, color=as.factor(N_t))) +
  geom_line(size=0.5) +
  labs(x = "Nt", y = "S(z, t)", color="Population") +
  labs(x="size z")+
  ggtitle("Plot of S(z, t) for values of N")+
  theme_bw()+
  theme(axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        axis.text.y=element_text(size=20),
        axis.text.x=element_text(size=13),
        legend.position = "none")
# strip.text.x = element_text(size = 15),
# legend.title = element_text( size = 22),
# legend.text = element_text(size = 14))


plot_R <- ggplot(data = results, aes(x = z, y = R, color=as.factor(N_t))) +
  geom_line(size=0.5) +
  labs(x = "Nt", y = "R(z, t)", color="Population") +
  labs(x="size z")+
  ggtitle("Plot of R(z, t) for values of N")+
  theme_bw()+
  theme(axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        axis.text.y=element_text(size=20),
        axis.text.x=element_text(size=13),
        # strip.text.x = element_text(size = 15),
        # legend.title = element_text( size = 22),
        legend.text = element_text(size = 14))


R_S_versus_N <- ggarrange(plot_S, plot_R, nrow=1)
R_S_versus_N


