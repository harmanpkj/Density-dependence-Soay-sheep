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
  p <- fftw::planFFT(max.offspring)
  pad.mat <- matrix(0, nrow = max.offspring - dim(Klist[[1]])[1], ncol = n.stages)
  Klist.pad <- Klist # set kappahat.pad as a list as Klist
  for(i in seq_len(n.ages)){
    Klist.pad[[i]] <- rbind(Klist[[i]], pad.mat) # pad the fertility matrix to desired size
  }
  kappahat.mat <- Klist.pad # set kappahat.mat as a list as Klist.pad
  for(i in seq_len(n.ages)){
    Klist.pad[[i]][is.na(Klist.pad[[i]])] <- 0
    for (j in seq_len(n.stages)) {
      kappahat.mat[[i]][,j] <- fftw::FFT(Klist.pad[[i]][,j], plan=p)
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

  gamma <- Re(fftw::IFFT(betahat.mat[, (initial.age - 1) * stage_per_age + 1], plan = p))
  if(end.age < initial.age){
    end.age <- initial.age
  }
  for (i in 2:(stage_per_age * end.age)) {
    gamma <- cbind(gamma, Re(fftw::IFFT(betahat.mat[,i], plan = p)))
  }
  #write.table(gamma, paste("age", (end.age- initial.age + 1), "X", stage_per_age, "stages.roedeer.block_fft", max.offspring,".txt", sep = ""), quote = F, row.names = F, col.names = F)
  return(gamma)
}
