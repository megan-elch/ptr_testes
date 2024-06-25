library(tidyverse)
library(cmdstanr)
library(posterior)
source("hs/model_functions.R")
source("hs/go_functions.R")
clusters = c("EC", "LC", "PTM", "SPC", "SPG", "St")
clusters_df = data.frame(ct = clusters, celltype_num = 1:length(clusters))

# specify output files, load in prepped data
output_path = paste0("/model_output/")
load(paste0("/processed_data/processed_data_hs.RData"))
prep_list = model_prep(mrna_suff, protein_suff, n_protein = 5, n_mrna = 5, 
                       clusters = clusters)

# load in stan fit
fit = readRDS(file = paste0(output_path, "fit.rds"))
fit$code()
fit$diagnostic_summary()

###############################################################################
########################### EXTRACT POSTERIOR DRAWS ###########################
###############################################################################
# extract posterior draws: mrna
mu = as_draws_df(fit$draws("mu")) %>%
  pivot_longer(cols = setdiff(colnames(.), c(".iteration", ".chain", ".draw")), values_to = "mu")  %>%
  dplyr::transmute(mu = mu, .chain = .chain, .iteration = .iteration,
                   param_id = unlist(str_extract_all(name, "\\[[^()]+\\]"))) %>%
  dplyr::mutate(param_id = substring(param_id, 2, nchar(param_id)-1)) %>%
  separate(param_id, c("uni_map","celltype_num"), sep = ',') %>%
  merge(prep_list$gene_map) %>%
  merge(clusters_df)

# extract posterior draws: rptr
r = as_draws_df(fit$draws("r")) %>%
  pivot_longer(cols = setdiff(colnames(.), c(".iteration", ".chain", ".draw")), values_to = "r")  %>%
  dplyr::transmute(r = r, .chain = .chain, .iteration = .iteration,
                   param_id = unlist(str_extract_all(name, "\\[[^()]+\\]"))) %>%
  dplyr::mutate(param_id = substring(param_id, 2, nchar(param_id)-1)) %>%
  separate(param_id, c("uni_map","celltype_num"), sep = ',')

# merge, define protein as mu + r and save
posterior_draws = merge(mu, r) %>%
  dplyr::mutate(prot = mu + r)
save(posterior_draws, file = paste0(output_path, "posterior_draws.RData"))

############################################################################
########################### GENE-LEVEL ANALYSIS ############################
############################################################################
# run gene-level test
gene_res = test_genes_centered(posterior_draws)
save(gene_res, file = paste0(output_path, "gene_res.RData"))

# compute across cluster mrna, protein correlation
cor_comparison = compute_correlation_comparison(gene_res = gene_res, posterior_draws = posterior_draws, prep_list = prep_list)
save(cor_comparison, file = paste0(output_path, "cor_comparison.RData"))

# compute statistics about variance of mrna and protein across clusters
posterior_signal = posterior_draws %>%
  dplyr::group_by(UNIPROT, .iteration, .chain) %>%
  dplyr::summarise(mu_var = var(mu, na.rm = T),
                   prot_var = var(prot, na.rm = T),
                   signal_cor = cor(mu, prot, use = "pairwise.complete.obs"),
                   r_var = var(r, na.rm = T)) %>%
  ungroup() %>%
  dplyr::group_by(UNIPROT) %>%
  dplyr::summarise(signal_var_mean = mean(mu_var/prot_var, na.rm = T),
                   signal_var_var = var(mu_var/prot_var, na.rm = T),
                   cor_mean = mean(signal_cor, na.rm = T),
                   cor_var = var(signal_cor, na.rm = T),
                   cor_sq_mean = mean(signal_cor^2, na.rm = T),
                   cor_sq_var = var(signal_cor^2, na.rm = T),
                   mu_var_mean = mean(mu_var, na.rm = T),
                   mu_var_var = var(mu_var, na.rm = T),
                   prot_var_mean = mean(prot_var, na.rm = T),
                   prot_var_var = var(prot_var, na.rm = T),
                   r_var_mean = mean(r_var, na.rm = T),
                   r_var_var = var(r_var, na.rm = T))
save(posterior_signal, file = paste0(output_path, "posterior_signal.RData"))

############################################################################
########################### GO TESTING AND PLOTS ###########################
############################################################################
# load predefined list of UNIPROT IDs and GO terms
load(paste0("/processed_data/GO_prep.RData"))

# run go testing
go_results = run_go_test_hs(posterior_draws = posterior_draws, GO_info, GO_info2, GO_tab)
save(go_results, file = paste0(output_path, "go_results_centered.RData"))

# additional summary information about go testing that helps to draw plots
posterior_draws_go = go_results$posterior_draws_go
save(posterior_draws_go, file = paste0(output_path, "go_results_draws_centered.RData"))
test_res = go_results$test_res
save(test_res, file = paste0(output_path, "go_test_results_centered.RData"))

for(i in 1:length(clusters)){
  GO_dat = compute_ticks_draws(posterior_draws_go = posterior_draws_go, test_res = test_res, celltype = clusters[i])
  save(GO_dat, file = paste0(output_path, "GO_dat_", clusters[i], ".RData"))
}

####################################################################################
########################### COMPUTE COMPLEX CORRELATIONS ###########################
####################################################################################
# load complex reference
complex_info = readr::read_delim("processed_data/coreComplexes.txt", delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)

# format complex df
complex_info = complex_info %>% 
  dplyr::group_by(ComplexName) %>% 
  dplyr::reframe(UNIPROT = unlist(strsplit(`subunits(UniProt IDs)`, ";", fixed = T)))

# load relevant output
load(paste0(output_path, "posterior_draws.RData"))
load(paste0(output_path, "gene_res.RData"))

# only use complexes associated with gene products that we observe
uni_unique = posterior_draws %>% pull(UNIPROT) %>% unique()
complex_info = complex_info %>% filter(UNIPROT %in% uni_unique)

# correlate mRNA-protein complex-level averages across cluters
complex_correlations = posterior_draws %>%
  merge(complex_info) %>%
  dplyr::group_by(ComplexName, ct, .iteration, .chain) %>%
  dplyr::summarise(mrna_go = mean(mu, na.rm = T),
                   prot_go = mean(prot, na.rm = T)) %>%
  ungroup() %>%
  dplyr::group_by(ComplexName, .iteration, .chain) %>%
  dplyr::summarise(cor_go = cor(mrna_go, prot_go, use = "pairwise.complete.obs"),
                   mrna_var = var(mrna_go),
                   prot_var = var(prot_go)) %>%
  ungroup() 

save(complex_correlations, file = paste0(output_path, "complex_correlations.RData"))

# compute within complex correlations
within_complex_correlations_a = posterior_draws %>%
  merge(complex_info) %>%
  dplyr::mutate(uniprot_a = UNIPROT,
                mu_a = mu,
                prot_a = prot) %>%
  dplyr::group_by(ComplexName) %>%
  dplyr::mutate(nuni = length(unique(uniprot_a))) %>%
  filter(nuni > 1) %>%
  ungroup() %>%
  dplyr::select(-c(UNIPROT, mu, prot, r, celltype_num, uni_map))

# for each complex, merge complex output with itself and compute pairwise correlations
complex_list = unique(within_complex_correlations_a$ComplexName)
within_complex_correlations = list()
for(i in 1:length(complex_list)){
  complex = within_complex_correlations_a %>%
    filter(ComplexName == complex_list[i])
  
  # merge same df
  within_complex_correlations[[i]] =  complex %>%
    mutate(uniprot_b = uniprot_a,
           mu_b = mu_a,
           prot_b = prot_a) %>%
    dplyr::select(-c(uniprot_a, mu_a, prot_a,)) %>%
    merge(complex) %>%
    # remove pairs where gene product is being correlated with itself
    filter(uniprot_a != uniprot_b) %>%
    dplyr::group_by(ComplexName, uniprot_a, uniprot_b, .iteration, .chain) %>%
    # compute correlation across clusters
    dplyr::summarise(mrna_cor = cor(mu_a, mu_b, use = "pairwise.complete.obs"),
                     prot_cor = cor(prot_a, prot_b, use = "pairwise.complete.obs")) %>%
    ungroup() %>%
    # compute median across subunits
    dplyr::group_by(ComplexName, .iteration, .chain) %>%
    dplyr::summarise(mrna_cor = median(mrna_cor),
                     prot_cor = median(prot_cor))
  print(i)
  
  }

# save object
within_complex_correlations = rlist::list.rbind(within_complex_correlations)  
save(within_complex_correlations, file = paste0(output_path, "within_complex_correlations.RData"))