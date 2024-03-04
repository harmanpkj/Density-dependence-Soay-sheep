# Calculate the lifetime reproductive success distribution for age + stage cases
# by treating each age stage combination as a stage, and then use the stage only method.


lrs_fft_agestage <- function(kappamat, U, stage_per_age,
                             initial.age = 1,
                             end.age = 1){
  ####===== FFT method to calculate the LRS distribution
  kappamat[is.na(kappamat)] <- 0
  max.offspring <- dim(kappamat)[1]
  n.stages <- dim(kappamat)[2]
  set.plan <- fftw::planFFT(max.offspring)
  alphahat.mat <- betahat.mat <-matrix(0, nrow = max.offspring, ncol = n.stages)
  for (i in 1:n.stages) {
    alphahat.mat[,i] <- fftw::FFT(kappamat[,i], plan = set.plan)
  }
  I <- diag(1, n.stages)
  e <- rep(1, n.stages)
  d <- (I - t(U)) %*% e
  for (i in 1:max.offspring) {
    W <- diag(alphahat.mat[i,])
    B1 <- solve(I - t(U) %*% W) # take about 39 seconds
    betahat.mat[i,] <- W %*% B1 %*% d
  }
  gamma <- Re(fftw::IFFT(betahat.mat[, (initial.age - 1) * stage_per_age + 1], plan = set.plan))
  if(end.age < initial.age){
    end.age <- initial.age
  }
  for (i in 2:(stage_per_age * end.age)) {
    gamma <- cbind(gamma, Re(fftw::IFFT(betahat.mat[,i], plan = set.plan)))
  }

  write.table(gamma, paste((end.age - initial.age + 1), "agesX", stage_per_age, "stages.roedeer.fft", max.offspring,".txt", sep = ""), quote = F, row.names = F, col.names = F)
  return(gamma) # gamma is a matrix including all initial stages
}
