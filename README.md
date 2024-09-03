# Density-dependence-Soay-sheep
The repository contains code to generate plots in the manuscript. Here we provide a detailed description to replicate the results in the manuscript.

Data file and RScripts:
1. "rand.params.csv"
2. "load_LRS_functions.R"
3. "Density_dep_soaysheep_code_Sep2024.R"

Download the above files and save in your working directory. Brief summary of files below:

"rand.params.csv" contains 10000 parameter sets (sampled from the covariance matrix) as rows and has 16 columns each of whom correspond to parameters of interest.

"load_LRS_functions.R" contains necessary functions to execute the LRS analysis.

"Density_dep_soaysheep_code_Sep2024.R" is the main file to run the analysis in the manuscript. For ease of replication, we have commented each segment of code and the what Figure it generates. This file sources "load_LRS_functions.R". Please ensure all three files are in the same folder.

Now launch the main code for the manuscript titled "Density_dep_soaysheep_code_Sep2024.R" which provides the step by step commands for PCA and LRS analysis. We work with a small sample of the original dataset in the code and thus sample 250 rows from the original dataset of 10000 rows. The code begins with defining IPM functions and generates equilibrium carrying capacities for each of the parameter sets. Covariates for PCA are defined and we examine results from PCA and LRS distributions to examine tradeoffs at different densities.

Additional Data and RScript for delifing. This is to replicate the analysis for another species of interest and create a rand.params for your species. After this, you can follow the proceedure in the main code "Density_dep_soaysheep_code_Sep2024.R" to carry out further analysis.
1. "sheep data 1986 to 1996.csv"
2. "delifing_code.R"

