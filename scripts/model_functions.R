#########################################################################################################################################
#                                                    FUNCTIONS RELATED TO RUNNING MODEL 
#########################################################################################################################################

# function to compute necessary quantities to run model
compute_suff_stats = function(df, df_info, type = "mRNA", clusters = NULL, pop_label = NULL, 
                              already_logged = T, protein_label = "prot", peptide_label = "pep", 
                              omit_zeros = FALSE){
 df_info = na.omit(df_info)
  if(tolower(type) == "mrna"){
    df = df %>%
      as.matrix() %>%
      data.frame()
    if(omit_zeros == TRUE){
        df = replace(df, df==0, NA) # 0s treated as NAs
    }
    df$SYMBOL = rownames(df)
    # compute sum of transcript counts
    means = df %>%
      pivot_longer(cols = intersect(df_info$id, colnames(df)), names_to = "id") %>%
      dplyr::mutate(value = as.numeric(value)) %>%
      base::merge(df_info) %>% # merge cell type labels
      filter(ct %in% clusters) %>% # filter observations to be in correct clusters
      dplyr::group_by(id) %>%
      dplyr::mutate(ts_counts = sum(value[is.finite(value)])) %>% # total number of transcripts in each cell
      ungroup() %>%
      group_by(SYMBOL, ct) %>% 
      dplyr::summarise(n_cells = sum(value > 0),
                       mrna_sum = sum(value, na.rm = T), # models sum of transcripts counts for each gene, cell type
                       mrna_av = mean(value, na.rm = T),
                       mrna_mean_log2 = mean(log2(value + 1), na.rm = T),
                       counts = sum(ts_counts)) %>% # total number of transcript counts across cells in each cluster
      ungroup()
      
      if(!is.null(pop_label)){
          means = means %>% mutate(pop_mrna = pop_label)
      }
  }

  if(tolower(type) == "protein"){
    means = df %>% 
      pivot_longer(cols = intersect(df_info$id, colnames(df)), names_to = "id")
    if(already_logged == FALSE){ # take log2 if needed
        means = means %>% mutate(value = log2(value))
    }  
    means[["UNIPROT"]] = means[[protein_label]]
    means[["pep"]] = means[[peptide_label]]
    # compute peptide level averages
    means = means %>%
            base::merge(df_info) %>%
            dplyr::group_by(UNIPROT, ct) %>% 
            dplyr::mutate(n_cells_prot = length(unique(id[is.finite(value)]))) %>% # number of cells observed for each protein cell type
            ungroup() %>%
            dplyr::group_by(pep, ct) %>%
            dplyr::summarise(pep_av = mean(value, na.rm = T), # models average log2 peptide intensity (original data already log2 transformed)
                             n_cells = sum(is.finite(value)), # number of cells for each peptide, cell type
                             n_cells_prot = unique(n_cells_prot), # number of cells for each protein, cell type
                             UNIPROT = unique(as.character(UNIPROT)),
                             pep_sum = sum(value, na.rm = T),
                             pep_sum_sq = sum(value^2, na.rm = T)) %>% # sum and sum sq if we want to increment target likelihood directly
            filter(n_cells > 5, ct %in% clusters) %>% # filter observations to be in correct clusters
            dplyr::group_by(UNIPROT) %>%
            dplyr::mutate(npep_ct = length(unique(pep))) %>% # number of peptides observed across cell types
            ungroup()
      
    if(!is.null(pop_label)){
        means = means %>% mutate(pop_protein = pop_label)
     }
  }
  return(means)
}

# organize sufficient statistics
model_prep = function(mrna_suff, protein_suff, 
                      n_mrna = 2, n_protein = 2, min_ct = 3,
                      clusters, specific_genes = NULL){
    options(mc.cores = parallel::detectCores())
    
  # organize peptide level averages
    protein_suff = protein_suff %>% 
        filter(ct %in% clusters) %>%
        mutate(pop_factor_protein = as.numeric(factor(pop_protein, levels = paste0("Pop", 1:n_protein)))) %>%
        mutate(measure_factor = as.numeric(factor(measure, levels = c("scope", "plex"))))
  # organize transcript counts   
    mrna_suff = mrna_suff %>% 
        filter(ct %in% clusters) %>%
        mutate(pop_factor_mrna = as.numeric(factor(pop_mrna, levels = paste0("Pop", 1:n_mrna))))
    
    # filter mrna to make sure its observed in enough cell types, nonzero variance
    mrna_suff = na.omit(mrna_suff) %>%
                dplyr::group_by(SYMBOL, pop_mrna)  %>%
                dplyr::mutate(n_ct = length(unique(ct)), zero_all = sum(mrna_sum) == 0) %>%
                ungroup() %>%
                filter(n_ct >= min_ct & zero_all == FALSE) %>%
                dplyr::select(-c(n_ct, zero_all))
                
    print(head(mrna_suff))

    # identify unique cell types in mrna data set for each gen
    mrna_cts = mrna_suff %>%
               dplyr::group_by(UNIPROT) %>%
               dplyr::reframe(ct_m = unique(ct)) %>% # keeping track of number of cell types
               ungroup()

    # filter protein to make sure its observed in enough cell types, nonzero variance
    protein_suff = na.omit(protein_suff) %>%
                   dplyr::group_by(UNIPROT, pop_protein) %>%
                   dplyr::mutate(n_ct = length(unique(ct)),
                                 var_av = var(pep_av, na.rm = T)) %>% # compute number of cell types observed for each gene, data set
                   ungroup() %>%
                   dplyr::filter(n_ct >= min_ct & is.finite(log(var_av))) %>% # must have enough cell types observed
                   dplyr::select(-n_ct)
                   
    print(head(protein_suff))

    # identify unique cell types in protein data set for each gen
    protein_cts = protein_suff %>%
                  dplyr::group_by(UNIPROT) %>%
                  dplyr::reframe(ct_p = unique(ct)) %>% # keeping track of number of cell types
                  ungroup()

    # keep proteins observed with enough cell types in both modalities
    cts = merge(mrna_cts, protein_cts) %>%
          dplyr::group_by(UNIPROT) %>%
          dplyr::summarise(n_ct_int = length(intersect(unique(ct_p), unique(ct_m)))) %>% # compute number of cell types shared across modalities
          ungroup() %>%
          dplyr::filter(n_ct_int >= min_ct) %>% # keep only genes observed in enough cell types 
          pull(UNIPROT)
          
    celltypes_union = intersect(mrna_suff$ct, protein_suff$ct) # keep shared cell types
          
    # identify proteins observed in both modalities
    proteins_union = protein_suff %>%
                     filter(UNIPROT %in% mrna_suff$UNIPROT, UNIPROT %in% cts) %>% # keep shared proteins in enough cell types
                     pull(UNIPROT)
                     
    print(length(proteins_union))
              
    # or select specific features to model       
    if(!is.null(specific_genes)){
        proteins_union = intersect(proteins_union, specific_genes)
    }
  
    # map observations to UNIPROT labels
    UNIPROT_map = mrna_suff %>% # make reference df to keep track of gene identities in model
                  dplyr::filter(UNIPROT %in% proteins_union, ct %in% celltypes_union, UNIPROT %in% cts) %>%
                  dplyr::group_by(UNIPROT) %>%
                  dplyr::summarise() %>%
                  ungroup() %>%
                  dplyr::mutate(uni_map = 1:length(unique(UNIPROT))) # id number to go with each gene      
                  
    print(paste("number of genes in model:", nrow(UNIPROT_map)))              
  
    # filter to shared cell types
    mrna_suff = mrna_suff %>% filter(ct %in% celltypes_union) # cell type filtering 
    protein_suff = protein_suff %>% filter(ct %in% celltypes_union) # cell type filtering 

    # celltype mapping, population mapping
    protein_suff_sample = merge(protein_suff, UNIPROT_map) %>% # merge to remove all filtered out obs
                          arrange(pop_protein) %>%
                          mutate(ct_factor = as.numeric(factor(ct, levels = clusters)),
                                 protein_obs_id = as.character(row_number())) # make factor for cell type ids
    print(paste("Number of Peptide Observations:", nrow(protein_suff_sample)))
         
    mrna_suff_sample = merge(mrna_suff, UNIPROT_map) %>%
                       arrange(pop_mrna) %>%
                       mutate(ct_factor = as.numeric(factor(ct, levels = clusters)),
                              mrna_obs_id = as.character(row_number()))
                       
    print(paste("Number of mRNA Observations:", nrow(mrna_suff_sample)))
    # output list with prepared data sets
    list(mrna = mrna_suff_sample, protein = protein_suff_sample, gene_map = UNIPROT_map, clusters = clusters)
}

# prepare data to be input to stan model
output_stan_data = function(prep_list){
    # use relevant quantities from prep function
    mrna = prep_list$mrna
    protein = prep_list$protein
    gene_map = prep_list$gene_map
    clusters = prep_list$clusters
    total_obs_mrna = nrow(mrna)
    total_obs_protein = nrow(protein)
    
    protein$n_cells = protein$npep_ct
    
    # organized data to input in stan model
    stan_input = list(n_g = nrow(gene_map), # number of genes
                  total_obs_mrna = total_obs_mrna, # total number of observations (mrna)
                  total_obs_protein = total_obs_protein, # total number of observations (protein)
                  mrna_sum = mrna$mrna_sum, # observed mrna transcript sum                       
                  protein_av = protein$pep_av, # observed log2 peptide averages
                  mrna_gene_rec = as.numeric(mrna$uni_map), # keep track of gene associated w each observation      
                  protein_gene_rec = as.numeric(protein$uni_map),  # keep track of protein associated w each observation
                  mrna_celltype_rec = as.numeric(mrna$ct_factor), # keep track of cell types mrna
                  protein_celltype_rec = as.numeric(protein$ct_factor), # keep track of cell types protein
                  mrna_pop_rec = as.numeric(mrna$pop_factor_mrna), # keep track of mrna data set
                  protein_pop_rec = as.numeric(protein$pop_factor_protein), # keep track of protein data sets
                  protein_measure_rec = as.numeric(protein$measure_factor),
                  n_mrna = mrna$counts, # total transcript counts
                  n_protein = protein$n_cells) # number of cells protein
    return(stan_input)
}
                      
# test for extreme r at gene level, using 95 pct interval or expected proportion of false discoveries
# put in posterior draws of mu, r, mu + r, or stan fit directly
# outputs data frame containing posterior summaries and gene, celltype test results
test_genes_centered = function(posterior_draws, threshold = 0.01){
    gene_res = posterior_draws %>%
               dplyr::group_by(UNIPROT, .iteration, .chain) %>%
               dplyr::mutate(r_gene = mean(r)) %>%
               ungroup() %>%
               dplyr::mutate(r_unscaled = r,
                             r = r - r_gene) %>%
               dplyr::group_by(UNIPROT, ct) %>%
               dplyr::summarise(p_r = ifelse(mean(r < 0) <= 1 - mean(r < 0), mean(r < 0), 1 - mean(r < 0)),
                                p_ru = ifelse(mean(r_unscaled < 0) <= 1 - mean(r_unscaled < 0), mean(r_unscaled < 0), 1 - mean(r_unscaled < 0)),                                mu_av = mean(mu, na.rm = T),
                                mu_lwr = quantile(mu, 0.025, na.rm = T),
                                mu_upr = quantile(mu, 0.975, na.rm = T),
                                mu_med = median(mu, na.rm = T),
                                r_av = mean(r, na.rm = T), # posterior summaries r
                                r_lwr = quantile(r, 0.025, na.rm = T),
                                r_upr = quantile(r, 0.975, na.rm = T),
                                r_med = median(r, na.rm = T),
                                r_sd = sd(r, na.rm = T),
                                ru_av = mean(r_unscaled, na.rm = T), # posterior summaries r
                                ru_lwr = quantile(r_unscaled, 0.025, na.rm = T),
                                ru_upr = quantile(r_unscaled, 0.975, na.rm = T),
                                ru_med = median(r_unscaled, na.rm = T),       
                                ru_sd = sd(r_unscaled, na.rm = T),
                                prot_av = mean(prot, na.rm = T), # posterior samples prot
                                prot_lwr = quantile(prot, 0.025, na.rm = T),
                                prot_upr = quantile(prot, 0.975, na.rm = T),
                                prot_med = median(prot, na.rm = T),
                                significant = !(r_lwr <= 0 & r_upr >= 0),
                                significant_u = !(ru_lwr <= 0 & ru_upr >= 0)) %>% # interval size
                   ungroup() %>%
                   dplyr::group_by(UNIPROT) %>%
                   dplyr::mutate(n_ct = length(intersect(unique(ct[is.finite(mu_av)]), unique(ct[is.finite(prot_av)])))) %>%
                   ungroup() %>%
                   dplyr::group_by(ct) %>%
                   arrange(p_r, .by_group = T) %>% # compute expected proportion of false discoveries for r, mu, mu + r
                   dplyr::mutate(fdr = cummean(p_r[order(p_r)]), # expected proportion
                                 significant_fdr = fdr <= threshold) %>%
                   ungroup()
    return(gene_res)
}

# function that computes the across cell types correlation of mrna, protein for both observed and model fit data 
compute_correlation_comparison = function(gene_res, posterior_draws, prep_list, age = F, mrna_log2 = F, type = 1){
  
  # identify posterior intervals for mu and mu + r          
  fit_intervals = posterior_draws %>%
    dplyr::group_by(UNIPROT, .iteration, .chain) %>%
    dplyr::summarise(fit_cor = cor(mu, prot, use = "pairwise.complete.obs")) %>%
    ungroup()
  
  if(type == 1){
    fit_intervals = fit_intervals %>%
      dplyr::group_by(UNIPROT) %>%
      dplyr::summarise(param_cor = mean(fit_cor, na.rm = T),
                       param_cor_med = median(fit_cor, na.rm = T),
                       param_cor_lwr = quantile(fit_cor, 0.025, na.rm = T),
                       param_cor_upr = quantile(fit_cor, 0.975, na.rm = T),
                       sd_cor = sd(fit_cor, na.rm = T))
  }
  
  if(type == 2){
    fit_intervals = fit_intervals %>%
      dplyr::group_by(.chain, .iteration) %>%
      dplyr::summarise(param_cor = mean(fit_cor, na.rm = T),
                       param_cor_med = median(fit_cor, na.rm = T),
                       param_cor_lwr = quantile(fit_cor, 0.025, na.rm = T),
                       param_cor_upr = quantile(fit_cor, 0.975, na.rm = T))
  }
  # use protein and mrna averages for observed correlation                                 
  protein = prep_list$protein %>%
    dplyr::group_by(UNIPROT, ct, pop_protein) %>%
    dplyr::summarise(prot_av = mean(pep_av, na.rm = T))
  mrna = prep_list$mrna %>% 
    dplyr::mutate(mrna_av = log2(mrna_sum/counts + 0.0001)) %>%
    dplyr::select(UNIPROT, ct, pop_mrna, mrna_av)
  cor_comparison = protein %>% # merge modalities to compute correlation
    merge(mrna) %>% 
    distinct(.keep_all = TRUE) %>%
    dplyr::group_by(UNIPROT, pop_protein, pop_mrna) %>%
    dplyr::summarise(n_ct = length(unique(ct)), # record number of cell types
                     pop_cor = ifelse(n_ct > 2, # observed data correlation
                                      cor(mrna_av, prot_av, use = "pairwise.complete.obs"),
                                      NA)) %>%
    ungroup() %>%
    na.omit()

  if(type == 1){
    cor_comparison = cor_comparison %>%
      merge(fit_intervals) # combine observed and model fit correlations into one df
  }
  if(type == 2){
    cor_comparison = cor_comparison %>% 
      dplyr::group_by(pop_protein, pop_mrna) %>%
      dplyr::summarise(pop_cor = median(pop_cor, na.rm = T)) %>%
      merge(fit_intervals)
  }
  return(cor_comparison)
}

# function to compute total least squares slope
tls_fun = function(x, y){
 r = prcomp(cbind(x, y))
 v = r$rotation
 beta <- v[2,1]/v[1,1]
 intercept = r$center[2] - beta*r$center[1]
 list(beta = beta, intercept = intercept)
}
