rm(list=ls())
# Harman Jaggi
# Code for Ecology Letters manuscript on Sep 2024: Density dependence
# This part executes the functions written in two separate files by loading them and examining the data.

# Install the required packages
install.packages("RSpectra") # for faster eigenvalues for block matrices
install.packages("ggpubr") # to arrange plots into grids and annotate figure
install.packages("tidyverse", dependencies = T) # for dataframe manipulation
# dependencies for tidyverse include dplyr, purrr, tidyr
install.packages("ggfortify") # for the pca plot, autoplot function uses ggfortify
install.packages("corrplot") # for the correlation plots
install.packages("RColorBrewer") # for nice colors on IPM function plot
install.packages("ggplot2") # for plotting
install.packages("reshape2") # to melt data frames

# Load the required packages
library(tidyverse)

#####################

# Check your working directory
getwd()

# Set your working directory and save the functions "load_LRS_functions_file.R" and "load_PCA_SSD_functions_file.R"
# in the directory along with the file rand.params.csv
setwd("/Users/harmanjaggi/Documents/Research/LRS/FixedEnv/FixedEnv")

# Source the R Script "load_LRS_functions_file.E" from the working directory
# using the code in the line below
source("load_LRS_functions_file.R")
source("load_PCA_SSD_functions_file.R")

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
# The columns in dat_og correspond to number of parameters. These are 16 for our case.
# The rows in dat_og correspond to number of samples. We set nsim = 250.
# It is set to nsim = 5 to check the code and reduce runtime.
# Please edit it to 250 to replicate the figures in the manuscript.
set.seed(120)
nsim = 5

# Sample nsim number of rows from the original data set
sample1.dat <- dplyr::sample_n(dat_og, nsim)

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

# Vector for equilibrium ratio or equilibrium values
# eqm_val = 0 corresponds to population size N = 0
# eqm_val = 1/2 corresponds to population size N = K/2
# eqm_val = 1 corresponds to population size N = K

eqm_val <- c(0, 1/2, 1)
# eqm_val <- c(0, 1/2, 1, 1.2, 1.4)

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

library(ggfortify)
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
ggpubr::ggarrange(pca_plot_n0, pca_plot_k2, pca_plot_k, nrow=3)

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
corrplot::corrplot(cor_mat_sam1_N0, method="color", addCoef.col = "black", number.cex=0.8)
corrplot::corrplot(cor_mat_sam1_K_2, method="color", addCoef.col = "black", number.cex=0.8)
corrplot::corrplot(cor_mat_sam1_K, method="color", addCoef.col = "black", number.cex=0.8)
# dev.off()
# Run the above code for eqm values of 1.2K and 1.4K to get the correlatin matrices
# shown in the left panel of Figure A7

# pdf("corplot_sam1_beyondK.pdf")
# par(mfrow=c(1,2))
# corrplot::corrplot(cor_mat_sam1_2K1, method="color", addCoef.col = "black", number.cex=0.8)
# corrplot::corrplot(cor_mat_sam1_2K2, method="color", addCoef.col = "black", number.cex=0.8)
# corrplot::corrplot(cor_mat_sam1_2K3, method="color", addCoef.col = "black", number.cex=0.8)
# corrplot::corrplot(cor_mat_sam1_2K4, method="color", addCoef.col = "black", number.cex=0.8)
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
  # stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
  #              parse = TRUE,
  #              rr.digits = 2, coef.digits = 2, size = 3, label.x = 0.97) +
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

# detach(package:ggbiplot)
# detach(package:plyr)

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
    lambda_block <- Re(RSpectra::eigs(block_mat[[1]], 1)$values[1])

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
# detach(package:ggbiplot)
# detach(package:plyr)
# library(tidyverse)

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

# Execute the function to get results for mother's distribution
matF_return <- matF_plot_func(sample1.datK)

# Combine matrices into one data frame
combined_df <- do.call(rbind, lapply(matF_return, function(mat) {
  reshape2::melt(mat)
}))

# Overall range of values for colors
overall_range <- c(0, max(combined_df$value))

# Create a color palette
custom_palette <- c("gray90", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7",
                    "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061")


# Command for Figure 3 in the manuscript
# Plot each matrix with the same color scale breaks and custom palette
mother_plots <- lapply(matF_return, function(mat) {
  df <- reshape2::melt(mat)
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
combined_mother_plot <- ggpubr::ggarrange(mother_plots[[1]],
                                  mother_plots[[2]],
                                  mother_plots[[3]], nrow=1, common.legend = T, legend="right")
ggpubr::annotate_figure(combined_mother_plot, left = grid::textGrob("Distribution of Offspring size", rot = 90, vjust = 1, gp = grid::gpar(cex = 1.4)),
                 bottom = grid::textGrob("Distribution of Mothers size (scaled by SSD)", vjust=-6, gp = grid::gpar(cex = 1.4)))
combined_mother_plot

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

# Command for plotting gamma, that is, the probability of having no offspring,
# with average covariate functions

a <- ggplot(dat1_vitals_meanlrs, aes(Ss, plrs0))+
  geom_point(size=2.5, alpha=0.9, color="darkolivegreen3")+
  # geom_point(dat1_vitals_meanlrs, mapping=aes(Survival_small_ssd, mean_lrs), color="red", size=2.5)+
  # scale_y_continuous(sec.axis = sec_axis(~./10, name="R0")) +
  facet_wrap(~den_name)+
  theme_bw()+
  # stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
  #              parse = TRUE, color = "blue",
  #              rr.digits = 2, coef.digits = 2, size = 3) +

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
  # stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
  #              parse = TRUE, color = "blue",
  #              rr.digits = 2, coef.digits = 2, size = 3) +

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
  # stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
  #              parse = TRUE, color = "blue",
  #              rr.digits = 2, coef.digits = 2, size = 3) +

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
  labs(color="")+
  facet_wrap(~den_name)+
  geom_smooth(method='lm', formula= y~x, se=F,
              size=0.5, linetype="dashed", color="black")+
  # color="blue"  )+
  # scale_colour_manual(values = c("cyan4"))+
  theme_bw()+
  # stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
  #              parse = TRUE, color = "blue",
  #              rr.digits = 2, coef.digits = 2, size = 3) +

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
library(patchwork)
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
print(nsim)
# Set nsim to the same as before. We use nsim = 250 in the paper.
nsim = 5

mean.prm_K <- eqm_func(as.matrix(mean.prm))
dat <- as.matrix(c(as.matrix(mean.prm), as.matrix(mean.prm_K)))

# seed 120
eqm_val <- c(0,0.5,1)

# Matrix of mean params where each col is sims and each row is mean of params
mean.prm.mat <- as.data.frame(matrix(as.matrix(mean.prm), nrow=nrow(mean.prm), ncol=nsim))
dim(mean.prm.mat)

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

# detach(package:plyr)
# detach(package:ggbiplot)
## Note: This is incase group_by doesnt work and is masked by a package. You can try detaching the plyr and ggbiplot.

# Compute average for lrs distributions
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

ggpubr::ggarrange(plot_test_all, plot_test_s, plot_test_sg, plot_test_r, nrow=2, ncol=2)

ggpubr::ggarrange(plot_test_s, plot_test_g, plot_test_sg, nrow=3, ncol=1)


# Define a range of N_t values
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

library(RColorBrewer)
colors <- brewer.pal("Dark2", n = 4)

# Plot survival and recruitment function at different densities
p1 <- ggplot2::ggplot(data, aes(x = z, y = S, color = factor(N_t), group = N_t)) +
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

p2 <- ggplot2::ggplot(data1, aes(x = z, y = R, color = factor(N_t), group = N_t)) +
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

ggpubr::ggarrange(p1, p2)

