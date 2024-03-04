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
