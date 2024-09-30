# install libraries
library(tidyverse)
library(cmdstanr)
library(posterior)

# utility functions related to running model
source("scripts/model_functions.R")
# cell types used in analysis
clusters = c("EC", "LC", "PTM", "SPC", "SPG", "St")

output_path = paste0("model_output/")

# processed data for running model (see prep file)
load(paste0("data/suff_stats.RData"))

# ############################################################################
# ############################ MODEL PREP AND EDA ############################ 
# ############################################################################
# run model prep
prep_list = model_prep(mrna_suff, protein_suff, n_protein = 5, n_mrna = 5, clusters = clusters)
save(prep_list, file = paste0(output_path, "prepped_data.RData"))

# compile model
file <- file.path(cmdstan_path(), "hs.stan")
mod <- cmdstan_model(file)
mod$print()

# prepare stan input
stan_data = output_stan_data(prep_list)

# run model and save output to output_path
fit <- mod$sample(
  data = stan_data,
  seed = 123,
  chains = 10,
  parallel_chains = 10,
  init = 0,
  refresh = 30,
  iter_warmup = 600,
  iter_sampling = 200
)
fit$save_object(file = paste0(output_path, "fit.rds"))
