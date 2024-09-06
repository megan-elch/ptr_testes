#########################################################################################################################################
#                                              FUNCTIONS RELATED TO POSTERIOR PREDICTIVE CHECKING 
#########################################################################################################################################
# utility functions for posterior predictive checks
prepare_zscore_ppc = function(ppc_df, type = "mrna"){
  if(tolower(type) == "mrna"){
    means = ppc_df %>%
      dplyr::mutate(mrna_bar_rep = log2(mrna_bar_rep + 0.0001)) %>%
      dplyr::group_by(UNIPROT, pop_mrna, .iteration, .chain) %>%
      dplyr::mutate(mrna_zscore = scale(mrna_bar_rep)[,1]) %>%
      ungroup() %>%
      dplyr::group_by(UNIPROT, ct, pop_mrna) %>%
      dplyr::summarise(mrna_ppc_zscore = mean(mrna_zscore, na.rm = T),
                       mrna_ppc_lwr = quantile(mrna_zscore, 0.025, na.rm = T),
                       mrna_ppc_upr = quantile(mrna_zscore, 0.975, na.rm = T),
                       mrna_ppc_med = median(mrna_zscore, na.rm = T),
                       mrna_ppc_unscaled = mean(mrna_bar_rep, na.rm = T),
                       mrna_ppc_lwr_unscaled = quantile(mrna_bar_rep, 0.025, na.rm = T),
                       mrna_ppc_upr_unscaled = quantile(mrna_bar_rep, 0.975, na.rm = T),
                       mrna_ppc_med_unscaled = median(mrna_bar_rep, na.rm = T)) %>%
      ungroup()
  }
  if(tolower(type) == "protein"){
    means = ppc_df %>%
      dplyr::group_by(UNIPROT, pop_protein, .iteration, .chain) %>%
      dplyr::mutate(protein_zscore = scale(protein_bar_rep)[,1]) %>%
      ungroup() %>%
      dplyr::group_by(UNIPROT, ct, pop_protein) %>%
      dplyr::summarise(protein_ppc_zscore = mean(protein_zscore, na.rm = T),
                       protein_ppc_lwr = quantile(protein_zscore, 0.025, na.rm = T),
                       protein_ppc_upr = quantile(protein_zscore, 0.975, na.rm = T),
                       protein_ppc_med = median(protein_zscore, na.rm = T),
                       protein_ppc_unscaled = mean(protein_bar_rep, na.rm = T),
                       protein_ppc_lwr_unscaled = quantile(protein_bar_rep, 0.025, na.rm = T),
                       protein_ppc_upr_unscaled = quantile(protein_bar_rep, 0.975, na.rm = T),
                       protein_ppc_med_unscaled = median(protein_bar_rep, na.rm = T)) %>%
      ungroup()
  }
  return(means)    
}

# compute across cluster variance of posterior predictive draws, return posterior intervals
compute_ppc_variance = function(ppc_df, type = "mrna"){
  if(tolower(type) == "mrna"){
    vars = ppc_df %>%
      dplyr::mutate(mrna_av_ppc = log2(mrna_bar_rep_counts + 0.0001),
                    mrna_av_obs = log2(scaled_counts + 0.0001)) %>%
      dplyr::group_by(UNIPROT, pop_mrna, .iteration, .chain) %>%
      dplyr::summarise(ppc_var = var(mrna_av_ppc, na.rm = T),
                       obs_var = var(mrna_av_obs, na.rm = T)) %>%
      ungroup() %>%
      dplyr::group_by(UNIPROT, pop_mrna) %>%
      dplyr::summarise(ppc_var_lwr = quantile(ppc_var, 0.025, na.rm = T),
                       ppc_var_med = median(ppc_var, na.rm = T),
                       ppc_var_upr = quantile(ppc_var, 0.975, na.rm = T),
                       obs_var = median(obs_var, na.rm = T)) %>%
      ungroup()
  }
  
  if(tolower(type) == "protein"){
    vars = ppc_df %>%
      #  dplyr::group_by(UNIPROT, ct, pop_protein, .iteration, .chain) %>%
      #  dplyr::summarise(protein_bar_rep = mean(protein_bar_rep, na.rm = T),
      #                   protein_av = mean(pep_av, na.rm = T)) %>%
      ungroup() %>%
      dplyr::group_by(UNIPROT, pop_protein, .iteration, .chain) %>%
      dplyr::summarise(ppc_var = var(protein_bar_rep, na.rm = T),
                       obs_var = var(protein_av, na.rm = T)) %>%
      ungroup() %>%
      dplyr::group_by(UNIPROT, pop_protein) %>%
      dplyr::summarise(ppc_var_lwr = quantile(ppc_var, 0.025, na.rm = T),
                       ppc_var_med = median(ppc_var, na.rm = T),
                       ppc_var_upr = quantile(ppc_var, 0.975, na.rm = T),
                       obs_var = median(obs_var, na.rm = T)) %>%
      ungroup()            
  }
  return(vars)    
}

# extract tls slope and ratio of across clusters sd
extract_scale_ppc = function(mrna_ppc, protein_ppc){
  tls_fun = function(x, y){
    r = prcomp(cbind(x, y))
    v = r$rotation
    beta <- v[2,1]/v[1,1]
    intercept = r$center[2] - beta*r$center[1]
    list(beta = beta, intercept = intercept)
  }
  
  mrna_ppc = mrna_ppc %>%
    dplyr::group_by(pop_mrna, UNIPROT, .iteration, .chain) %>%
    dplyr::transmute(ct = ct, mrna_value = log2(mrna_sum_rep/counts + 0.0001), 
                     mrna_av = mean(mrna_value, na.rm = T), mrna_relative = mrna_value - mrna_av) %>%
    ungroup() %>%
    dplyr::mutate(pop_mrna = paste("mRNA", pop_mrna)) 
  
  scale_ppc = protein_ppc %>%
    dplyr::group_by(UNIPROT, ct, pop_protein, .iteration, .chain) %>%
    dplyr::summarise(protein_value = mean(protein_bar_rep, na.rm = T)) %>%
    ungroup() %>%
    dplyr::group_by(UNIPROT, pop_protein, .iteration, .chain) %>%
    dplyr::mutate(protein_av = mean(protein_value, na.rm = T), 
                  protein_relative = protein_value - protein_av) %>%
    merge(mrna_ppc) %>%
    filter(is.finite(mrna_relative) & is.finite(protein_relative)) %>%
    ungroup() %>%
    dplyr::group_by(UNIPROT, pop_protein, pop_mrna, .iteration, .chain) %>%
    dplyr::summarise(mrna_sd = sd(mrna_relative, na.rm = T),
                     protein_sd = sd(protein_relative, na.rm = T),
                     sd_ppc = protein_sd/mrna_sd,
                     tls_ppc = tls_fun(mrna_relative, protein_relative)$beta) %>%
    ungroup() %>%
    filter(is.finite(mrna_sd) & is.finite(protein_sd) & 
             is.finite(sd_ppc) & is.finite(tls_ppc)) %>%
    dplyr::group_by(UNIPROT, pop_protein, pop_mrna) %>%
    dplyr::summarise(mrna_sd_lwr = quantile(mrna_sd, 0.025, na.rm = T),
                     mrna_sd_upr = quantile(mrna_sd, 0.975, na.rm = T),
                     mrna_sd_med = median(mrna_sd, na.rm = T),
                     mrna_ppc_av = mean(mrna_sd, na.rm = T),  
                     
                     prot_sd_lwr = quantile(protein_sd, 0.025, na.rm = T),
                     prot_sd_upr = quantile(protein_sd, 0.975, na.rm = T),
                     prot_sd_med = median(protein_sd, na.rm = T),
                     prot_ppc_av = mean(protein_sd, na.rm = T), 
                     
                     sd_ppc_lwr = quantile(sd_ppc, 0.025, na.rm = T),
                     sd_ppc_upr = quantile(sd_ppc, 0.975, na.rm = T),
                     sd_ppc_med = median(sd_ppc, na.rm = T),
                     sd_ppc_av = mean(sd_ppc, na.rm = T),  
                     
                     tls_ppc_lwr = quantile(tls_ppc, 0.025, na.rm = T),
                     tls_ppc_upr = quantile(tls_ppc, 0.975, na.rm = T),
                     tls_ppc_med = median(tls_ppc, na.rm = T),
                     tls_ppc_av = mean(tls_ppc, na.rm = T))
  
  return(scale_ppc)                             
}

# posterior predictive across cluster mrna protein correlation
compute_ppc_correlations = function(mrna_ppc, protein_ppc){
  mrna_ppc = mrna_ppc %>%
    dplyr::mutate(mrna_bar_rep = log2(mrna_bar_rep_counts + 0.0001),
                  mrna_bar = log2(scaled_counts + 0.0001)) %>%
    dplyr::group_by(UNIPROT, pop_mrna, .iteration, .chain) %>%
    dplyr::summarise(mrna_ppc_zscore = scale(mrna_bar_rep)[,1],
                     mrna_obs_zscore = scale(mrna_bar)[,1],
                     ct = ct) %>%
    ungroup()
  
  protein_ppc = protein_ppc %>%
    # dplyr::group_by(UNIPROT, ct, pop_protein, .iteration, .chain) %>%
    # dplyr::summarise(protein_bar_rep = mean(protein_bar_rep, na.rm = T),
    #                  prot_av = mean(pep_av, na.rm = T)) %>%
    ungroup() %>%
    dplyr::group_by(UNIPROT, pop_protein, .iteration, .chain) %>%
    dplyr::summarise(protein_ppc_zscore = scale(protein_bar_rep)[,1],
                     protein_obs_zscore = scale(protein_av)[,1],
                     ct = ct) %>%
    ungroup() 
  
  ppc_cor = merge(mrna_ppc, protein_ppc) %>%
    dplyr::group_by(UNIPROT, pop_protein, pop_mrna, .iteration, .chain) %>%
    dplyr::summarise(ppc_cor = ifelse(n() > 3,
                                      cor(mrna_ppc_zscore, protein_ppc_zscore, use = "pairwise.complete.obs"),
                                      NA)) %>%
    na.omit() %>%
    ungroup() 
  
  mrna_ppc = mrna_ppc %>% 
    dplyr::group_by(UNIPROT, pop_mrna, ct) %>%
    dplyr::summarise(mrna_obs_zscore = mean(mrna_obs_zscore, na.rm = T)) %>%
    ungroup()
  
  protein_ppc = protein_ppc %>%
    dplyr::group_by(UNIPROT, pop_protein, ct) %>%
    dplyr::summarise(protein_obs_zscore = mean(protein_obs_zscore, na.rm = T)) %>%
    ungroup()
  
  obs_cor = merge(mrna_ppc, protein_ppc) %>%
    dplyr::group_by(UNIPROT, pop_mrna, pop_protein) %>%
    dplyr::summarise(pop_cor = ifelse(n() > 3, 
                                      cor(mrna_obs_zscore, protein_obs_zscore, use = "pairwise.complete.obs"),
                                      NA)) %>%
    ungroup() %>%
    na.omit() 
  list(obs_cor = obs_cor, ppc_cor = ppc_cor) 
}
