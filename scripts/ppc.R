# load libraries
library(tidyverse)
library(cmdstanr)
library(posterior)
source("hs/model_functions.R")
source("hs/ppc_functions.R")
clusters = c("EC", "LC", "PTM", "SPC", "SPG", "St")

# specify output path and folder to put ppc checks
output_path = paste0("model_output/")
ppc_path = paste0("model_output/")

# upload fit
fit = readRDS(file = paste0(output_path, "fit.rds"))
fit$code()
load(paste0("data/suff_stats.RData"))

# model prepped values
prep_list = model_prep(mrna_suff, protein_suff, n_protein = 5, n_mrna = 5, clusters = clusters)
mrna = prep_list$mrna %>% mutate(scaled_counts = mrna_sum/counts)
stan_data = output_stan_data(prep_list)

# compile mrna posterior predictive stan file
file_mrna_gq <- file.path(cmdstan_path(), "mrna_ppc_prior.stan")
mod_mrna_gq <- cmdstan_model(file_mrna_gq)

# generate mrna posterior predictive draws
fit_mrna_gq <- mod_mrna_gq$generate_quantities(fit, data = stan_data, seed = 123)
fit_mrna_gq$save_object(file = paste0(ppc_path, "mrna_gq.rds"))

# extract mrna pp draws, format
mrna_ppc = as_draws_df(fit_mrna_gq$draws()) %>%
  pivot_longer(cols = setdiff(colnames(.), c(".iteration", ".chain", ".draw")), values_to = "mrna_sum_rep") %>%
  dplyr::transmute(mrna_sum_rep, .chain, .iteration,
                   mrna_obs_id = unlist(str_extract_all(name, "\\[[^()]+\\]"))) %>%
  dplyr::mutate(mrna_obs_id = substring(mrna_obs_id, 2, nchar(mrna_obs_id)-1)) %>%
  merge(transmute(prep_list$mrna, UNIPROT, ct, SYMBOL,
                  n_cells, mrna_sum, scaled_counts = mrna_sum/counts,
                  counts, pop_mrna)) %>%
  dplyr::mutate(mrna_bar_rep = mrna_sum_rep/counts,
                mrna_bar_rep_counts = mrna_sum_rep/counts)
save(mrna_ppc, file = paste0(ppc_path, "mrna_ppc.RData"))
# load(file = paste0(ppc_path, "mrna_ppc.RData"))

# compute z scores of mrna pp draws and output summary
mrna_ppc_z = prepare_zscore_ppc(mrna_ppc, type = "mrna")
save(mrna_ppc_z, file = paste0(ppc_path, "mrna_ppc_z.RData"))

# compute posterior predictive across cluster variance intervals
mrna_ppc_variance = compute_ppc_variance(mrna_ppc, type = "mrna")
save(mrna_ppc_variance, file = paste0(ppc_path, "mrna_ppc_variance.RData"))

# compile protein posterior predictive stan file
file_protein_gq <- file.path(cmdstan_path(), "protein_ppc.stan")
mod_protein_gq <- cmdstan_model(file_protein_gq)

# generate protein posterior predictive draws
fit_protein_gq <- mod_protein_gq$generate_quantities(fit, data = stan_data, seed = 123)
fit_protein_gq$save_object(file = paste0(ppc_path, "protein_gq.rds"))
rm(fit)

# extract protein pp draws, format
protein_ppc = as_draws_df(fit_protein_gq$draws()) %>%
  pivot_longer(cols = setdiff(colnames(.), c(".iteration", ".chain", ".draw")), values_to = "protein_bar_rep") %>%
  dplyr::transmute(protein_bar_rep, .chain, .iteration,
                   protein_obs_id = unlist(str_extract_all(name, "\\[[^()]+\\]"))) %>%
  dplyr::mutate(protein_obs_id = substring(protein_obs_id, 2, nchar(protein_obs_id)-1)) %>%
  merge(mutate(prep_list$protein))
save(protein_ppc, file = paste0(ppc_path, "protein_ppc.RData"))
rm(fit_protein_gq)

# compute z scores of protein pp draws and output summary
protein_ppc_z = prepare_zscore_ppc(protein_ppc, type = "protein")
save(protein_ppc_z, file = paste0(ppc_path, "protein_ppc_z.RData"))

# compute posterior predictive across cluster variance intervals
protein_ppc_variance = compute_ppc_variance(protein_ppc, type = "protein")
save(protein_ppc_variance, file = paste0(ppc_path, "protein_ppc_variance.RData"))

load(file = paste0(ppc_path, "mrna_ppc.RData"))
# posterior predictive across cluster mrna, protein correlations
ppc_cor_obj = compute_ppc_correlations(mrna_ppc, protein_ppc)
save(ppc_cor_obj, file = paste0(ppc_path, "ppc_cor.RData"))

# posterior predictive ratio of standard deviations
scale_ppc = extract_scale_ppc(mrna_ppc, protein_ppc)
save(scale_ppc, file = paste0(ppc_path, "scale_ppc.RData"))