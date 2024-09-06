library(tidyverse)
# library(cmdstanr)
# library(posterior)
source("hs/model_functions.R")
source("hs/go_functions.R")
clusters = c("EC", "LC", "PTM", "SPC", "SPG", "St")
clusters_df = data.frame(ct = clusters, celltype_num = 1:length(clusters))

# specify output files, load in prepped data
output_path = paste0("model_output/")
load(paste0("data/suff_stats.RData"))
prep_list = model_prep(mrna_suff, protein_suff, n_protein = 5, n_mrna = 5, 
                       clusters = clusters)

# since mu, and r, are formed as rectangular matrices, make sure that posterior samples
# for imputated gene products are not used for downstream analysis
present_combos_prot = prep_list$protein %>%
  dplyr::group_by(UNIPROT) %>%
  dplyr::reframe(ct = unique(ct))

present_combos = prep_list$mrna %>%
  dplyr::group_by(UNIPROT) %>%
  dplyr::reframe(ct = unique(ct)) %>%
  merge(present_combos_prot)

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
  merge(present_combos) %>%
  dplyr::mutate(prot = mu + r)
save(posterior_draws, file = paste0(output_path, "posterior_draws.RData"))

############################################################################
########################### GENE-LEVEL ANALYSIS ############################
############################################################################
# run gene-level test
gene_res = test_genes_centered(posterior_draws)
save(gene_res, file = paste0(output_path, "gene_res.RData"))

# posterior intervals for mrna and protein
posterior_intervals = posterior_draws %>%
  dplyr::group_by(UNIPROT, .iteration, .chain) %>%
  dplyr::mutate(mu_mean = mean(mu, na.rm = T),
                prot_mean = mean(prot, na.rm = T)) %>%
  ungroup() %>%
  dplyr::mutate(mu_ctr = mu - mu_mean, # center mrna by across ct average 
                prot_ctr = prot - prot_mean) %>% # center protein by across ct average
  dplyr::group_by(UNIPROT, ct) %>%
  dplyr::summarise(# summary stats for mrna
                   mu_av = mean(mu),
                   mu_med = median(mu),
                   mu_lwr = quantile(mu, 0.025),
                   mu_upr = quantile(mu, 0.975),
                   mu_lwr2 = quantile(mu, 0.25),
                   mu_upr2 = quantile(mu, 0.75),
                   # summary stats for centered mrna
                   mu_av_ctr = mean(mu_ctr),
                   mu_med_ctr = median(mu_ctr),
                   mu_lwr_ctr = quantile(mu_ctr, 0.025),
                   mu_upr_ctr = quantile(mu_ctr, 0.975),
                   mu_lwr_ctr2 = quantile(mu_ctr, 0.25),
                   mu_upr_ctr2 = quantile(mu_ctr, 0.75),
                   # summary stats for protein
                   prot_av = mean(prot),
                   prot_med = median(prot),
                   prot_lwr = quantile(prot, 0.025),
                   prot_upr = quantile(prot, 0.975),
                   prot_lwr2 = quantile(prot, 0.25),
                   prot_upr2 = quantile(prot, 0.75),
                   # summary stats for centered protein
                   prot_av_ctr = mean(prot_ctr),
                   prot_med_ctr = median(prot_ctr),
                   prot_lwr_ctr = quantile(prot_ctr, 0.025),
                   prot_upr_ctr = quantile(prot_ctr, 0.975),
                   prot_lwr_ctr2 = quantile(prot_ctr, 0.25),
                   prot_upr_ctr2 = quantile(prot_ctr, 0.75))
save(posterior_intervals, file = paste0(output_path, "mu_prot_intervals.RData"))

# test differences between pairs of genes
test_gene_pairs = posterior_draws %>%
  transmute(UNIPROT, ct_a = ct, r_a = r, .iteration, .chain) %>%
  merge(transmute(posterior_draws, UNIPROT, ct_b = ct, r_b = r, .iteration, .chain)) %>%
  filter(ct_a != ct_b) %>%
  mutate(ct_pair = purrr::map2_chr(ct_a, ct_b, ~toString(sort(c(.x, .y))))) %>%
  arrange(ct_a) %>%
  distinct(UNIPROT, ct_pair, .iteration, .chain, .keep_all = TRUE) %>%
  dplyr::mutate(r_diff = r_a - r_b)
save(test_gene_pairs, file = paste0(output_path, "matched_posterior_draws.RData"))
head(test_gene_pairs)

# summary stats for difference in ratios
test_gene_pairs = test_gene_pairs %>%
  dplyr::group_by(UNIPROT, ct_pair) %>%
  dplyr::summarise(r_diff_lwr = quantile(r_diff, 0.025, na.rm = T),
                   r_diff_upr = quantile(r_diff, 0.975, na.rm = T),
                   r_diff_med = median(r_diff, na.rm = T),
                   r_diff_mean = mean(r_diff, na.rm = T),
                   significant_pair = !(r_diff_lwr < 0 & r_diff_upr >= 0),
                   pep = ifelse(r_diff_mean > 0, mean(r_diff < 0), mean(r_diff > 0))) %>%
  ungroup() %>%
  dplyr::group_by(UNIPROT) %>%
  dplyr::mutate(significant_gene = sum(significant_pair) > 0) %>%
  ungroup() %>%
  arrange(pep) %>%
  dplyr::mutate(fdr = cummean(pep[order(pep)]), # expected proportion
                significant_fdr = fdr <= 0.05)
save(test_gene_pairs, file = paste0(output_path, "test_gene_pairs.RData"))

# compute across cluster mrna, protein correlation
cor_comparison = compute_correlation_comparison(gene_res = gene_res, posterior_draws = posterior_draws, prep_list = prep_list)
save(cor_comparison, file = paste0(output_path, "cor_comparison.RData"))

# compute statistics about variance of mrna and protein across clusters
posterior_signal = posterior_draws %>%
  dplyr::group_by(UNIPROT, .iteration, .chain) %>%
  dplyr::summarise(mu_var = var(mu, na.rm = T),
                   prot_var = var(prot, na.rm = T),
                   # mrna protein correlattion
                   signal_cor = cor(mu, prot, use = "pairwise.complete.obs"),
                   r_var = var(r, na.rm = T)) %>%
  ungroup() %>%
  dplyr::group_by(.iteration, .chain) %>%
  dplyr::mutate(across_uni_med = median(signal_cor, na.rm = T),
                across_uni_mean = mean(signal_cor, na.rm = T)) %>%
  ungroup() %>%
  dplyr::group_by(UNIPROT) %>%
  dplyr::summarise(signal_var_mean = mean(mu_var/prot_var, na.rm = T),
                   signal_var_var = var(mu_var/prot_var, na.rm = T),
                   cor_mean = mean(signal_cor, na.rm = T),
                   pep_cor_med = mean(signal_cor < across_uni_med),
                   pep_cor_mean = mean(signal_cor < across_uni_mean),
                   cor_lwr = quantile(signal_cor, 0.025, na.rm = T),
                   cor_upr = quantile(signal_cor, 0.975, na.rm = T),
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

##########################################################################
########################### GO TESTING AND PLOTS ###########################
############################################################################
# # load relevant output
load(paste0(output_path, "posterior_draws.RData"))
# posterior_draws = posterior_draws %>% merge(present_combos)
load(paste0(output_path, "gene_res.RData"))

# load predefined list of UNIPROT IDs and GO terms
load(paste0("/home/m.elcheikhali/hs/data/GO_prep.RData"))

# run go testing
go_results = run_go_test_hs(posterior_draws = posterior_draws, GO_info, GO_info2, GO_tab)
save(go_results, file = paste0(output_path, "go_results_centered.RData"))

# additional summary information about go testing that helps to draw plots
posterior_draws_go = go_results$posterior_draws_go
save(posterior_draws_go, file = paste0(output_path, "go_results_draws_centered.RData"))
test_res = go_results$test_res
save(test_res, file = paste0(output_path, "go_test_av.RData"))

# go level summary information
for(i in 1:length(clusters)){
  GO_dat = compute_ticks_draws(posterior_draws_go = posterior_draws_go, test_res = test_res, celltype = clusters[i])
  save(GO_dat, file = paste0(output_path, "GO_dat_", clusters[i], ".RData"))
}

# test mrna protein correlations at go level
go_correlations = posterior_draws_go %>%
  dplyr::group_by(UNIPROT, GO, .iteration, .chain) %>%
  dplyr::summarise(go_cor = cor(mu, prot, use = "pairwise.complete.obs")) %>%
  ungroup() %>%
  dplyr::group_by(GO, .iteration, .chain) %>%
  dplyr::summarise(go_cor = median(go_cor, na.rm = T)) %>%
  ungroup()

save(go_correlations, file = paste0(output_path, "go_test_corr.RData"))

# test across gene product variance
var_test_go = posterior_draws %>%
  filter(.chain <= 5) %>%
  merge(GO_info) %>%
  dplyr::group_by(GO, ct, .iteration, .chain) %>%
  dplyr::summarise(group_var = var(r, na.rm = T)) %>%
  ungroup() %>%
  dplyr::group_by(GO, ct) %>%
  dplyr::summarise(var_av = mean(group_var, na.rm = T),
                   var_med = median(group_var, na.rm = T),
                   var_lwr = quantile(group_var, 0.025, na.rm = T),
                   var_upr = quantile(group_var, 0.975, na.rm = T))
save(var_test_go, file = paste0(output_path, "var_test_go.RData"))

####################################################################################
########################### COMPUTE COMPLEX CORRELATIONS ###########################
####################################################################################
# load complex reference
complex_info = readr::read_delim("hs/data/coreComplexes.txt", delim = "\t", escape_double = FALSE,
                                 trim_ws = TRUE)

# format complex df
complex_info = complex_info %>%
  dplyr::group_by(ComplexName) %>%
  dplyr::reframe(UNIPROT = unlist(strsplit(`subunits(UniProt IDs)`, ";", fixed = T)))

# only use complexes associated with gene products that we observe
uni_unique = posterior_draws %>% pull(UNIPROT) %>% unique()
complex_info = complex_info %>% filter(UNIPROT %in% uni_unique)

# correlate mRNA-protein complex-level averages across cluters
complex_correlations = posterior_draws %>%
  dplyr::group_by(UNIPROT, .iteration, .chain) %>%
  dplyr::summarise(complex_cor = cor(mu, prot, use = "pairwise.complete.obs")) %>%
  ungroup() %>%
  merge(complex_info) %>%
  dplyr::group_by(ComplexName, .iteration, .chain) %>%
  dplyr::summarise(cor_complex = median(complex_cor, na.rm = T)) %>%
  ungroup()

save(complex_correlations, file = paste0(output_path, "complex_test_corr.RData"))

# run complex level test on average rptr
test_res_complex = posterior_draws %>%
  merge(complex_info) %>%
  dplyr::group_by(ComplexName, ct, .iteration, .chain) %>%
  dplyr::summarise(r_mean = mean(r, na.rm = T),
                   mu_mean = mean(mu, na.rm = T),
                   prot_mean = mean(prot, na.rm = T)) %>%
  ungroup() %>%
  dplyr::group_by(ComplexName, .iteration, .chain) %>%
  dplyr::mutate(r_complex_av = mean(r_mean, na.rm = T),
                r_centered_mean = r_mean - r_complex_av, # center across subunits
                mu_complex_av = mean(mu_mean, na.rm = T),
                mu_centered_mean = mu_mean - mu_complex_av,
                prot_complex_av = mean(prot_mean, na.rm = T),
                prot_centered_mean = prot_mean - prot_complex_av) %>%
  ungroup() %>%
  dplyr::group_by(ComplexName, ct) %>%
  dplyr::summarise(r_av = mean(r_mean, na.rm = T), # average r in group
                   r_centered_av = mean(r_centered_mean, na.rm = T),
                   r_lwr = quantile(r_mean, 0.025, na.rm = T), # 95 pct interval r
                   r_med = median(r_mean, na.rm = T),
                   r_upr = quantile(r_mean, 0.975, na.rm = T),
                   r_centered_lwr = quantile(r_centered_mean, 0.025, na.rm = T), # 95 pct interval r
                   r_centered_med = median(r_centered_mean, na.rm = T),
                   r_centered_upr = quantile(r_centered_mean, 0.975, na.rm = T),
                   p_r = ifelse(mean(r_centered_mean < 0) <= 1 - mean(r_centered_mean < 0), mean(r_centered_mean < 0), 1 - mean(r_centered_mean < 0)),
                   significant = !(r_centered_lwr < 0 & r_centered_upr > 0),
                   mu_av = mean(mu_mean, na.rm = T), # average mrna in group
                   mu_centered_av = mean(mu_centered_mean, na.rm = T),
                   mu_med = median(mu_mean, na.rm = T),
                   mu_centered_med = median(mu_centered_mean, na.rm = T),
                   p_mu = ifelse(mean(mu_centered_mean < 0) <= 1 - mean(mu_centered_mean < 0), mean(mu_centered_mean < 0), 1 - mean(mu_centered_mean < 0)),
                   prot_av = mean(prot_mean, na.rm = T), # average prot in group
                   prot_centered_av = mean(prot_centered_mean, na.rm = T),
                   prot_med = median(prot_mean, na.rm = T),
                   prot_centered_med = median(prot_centered_mean, na.rm = T),
                   p_prot = ifelse(mean(prot_centered_mean < 0) <= 1 - mean(prot_centered_mean < 0), mean(prot_centered_mean < 0), 1 - mean(prot_centered_mean < 0))) %>%
  ungroup()
save(test_res_complex, file = paste0(output_path, "complex_test_av.RData"))

# test complex level variance
var_test_complex = posterior_draws %>%
  merge(complex_info) %>%
  dplyr::group_by(ComplexName, ct, .iteration, .chain) %>%
  dplyr::summarise(group_var = var(r, na.rm = T)) %>%
  ungroup() %>%
  dplyr::group_by(ComplexName, ct) %>%
  dplyr::summarise(var_av = mean(group_var, na.rm = T),
                   var_med = median(group_var, na.rm = T),
                   var_lwr = quantile(group_var, 0.025, na.rm = T),
                   var_upr = quantile(group_var, 0.975, na.rm = T))
save(var_test_complex, file = paste0(output_path, "var_test_complex.RData"))

# collect pairs of genes for go groups
go_pairs = GO_info %>%
  filter(GO %in% test_res$GO & UNIPROT %in% posterior_draws$UNIPROT) %>%
  transmute(GO, UNIPROT_a = UNIPROT) %>%
  merge(transmute(GO_info, GO, UNIPROT_b = UNIPROT)) %>%
  filter(UNIPROT_a != UNIPROT_b) %>%
  mutate(gene_pair = purrr::map2_chr(UNIPROT_a, UNIPROT_b, ~toString(sort(c(.x, .y))))) %>%
  # keep only distinct pairs of genes
  distinct(gene_pair, GO, .keep_all = TRUE)

# function to compute correlation within proteins in each go group
compute_within_cor_go = function(go_set, posterior_draws){
  within_cor = go_set %>%
    merge(transmute(posterior_draws, UNIPROT_a = UNIPROT, prot_a = prot, .iteration, .chain, ct)) %>%
    merge(transmute(posterior_draws, UNIPROT_b = UNIPROT, prot_b = prot, .iteration, .chain, ct)) %>%
    dplyr::group_by(UNIPROT_a, UNIPROT_b, GO, .iteration, .chain) %>%
    dplyr::summarise(prot_cor = cor(prot_a, prot_b)) %>%
    ungroup() %>%
    dplyr::group_by(GO, .iteration, .chain) %>%
    dplyr::summarise(prot_cor = median(prot_cor, na.rm = T)) %>%
    ungroup() %>%
    dplyr::group_by(GO) %>%
    dplyr::summarise(prot_cor = median(prot_cor, na.rm = T))

  return(within_cor)
}

# compute correlation within proteins in go groups 
within_go_correlations = go_pairs %>%
  group_split(GO) %>%
  lapply(., function(x) compute_within_cor_go(x, posterior_draws)) %>%
  rlist::list.rbind()
save(within_go_correlations, file = paste0(output_path, "within_go_correlations.RData"))

# pairs of complex subunits to compute correlation 
go_pairs = complex_info %>%
  filter(ComplexName %in% test_res_complex$ComplexName & UNIPROT %in% posterior_draws$UNIPROT) %>%
  transmute(ComplexName, UNIPROT_a = UNIPROT) %>%
  merge(transmute(complex_info, ComplexName, UNIPROT_b = UNIPROT)) %>%
  filter(UNIPROT_a != UNIPROT_b) %>%
  mutate(gene_pair = purrr::map2_chr(UNIPROT_a, UNIPROT_b, ~toString(sort(c(.x, .y))))) %>%
  distinct(gene_pair, ComplexName, .keep_all = TRUE)

# function to compute correlation within protein complex subunits 
compute_within_cor_go = function(go_set, posterior_draws){
  within_cor = go_set %>%
    merge(transmute(posterior_draws, UNIPROT_a = UNIPROT, prot_a = prot, .iteration, .chain, ct)) %>%
    merge(transmute(posterior_draws, UNIPROT_b = UNIPROT, prot_b = prot, .iteration, .chain, ct)) %>%
    dplyr::group_by(UNIPROT_a, UNIPROT_b, ComplexName, .iteration, .chain) %>%
    dplyr::summarise(prot_cor = cor(prot_a, prot_b)) %>%
    ungroup() %>%
    dplyr::group_by(ComplexName, .iteration, .chain) %>%
    dplyr::summarise(prot_cor = median(prot_cor, na.rm = T)) %>%
    ungroup() %>%
    dplyr::group_by(ComplexName) %>%
    dplyr::summarise(prot_cor = median(prot_cor, na.rm = T))

  return(within_cor)
}

# compute correlation within proteins in protein complexes 
within_complex_correlations = go_pairs %>%
  group_split(ComplexName) %>%
  lapply(., function(x) compute_within_cor_go(x, posterior_draws)) %>%
  rlist::list.rbind()
save(within_complex_correlations, file = paste0(output_path, "within_complex_correlations.RData"))