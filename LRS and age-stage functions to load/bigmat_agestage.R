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
