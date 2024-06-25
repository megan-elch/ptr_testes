# source p functions
source("scripts/model_functions.R")

# point plot of mrna, protein log fold changes for each cell type
obs_mrna_prot = function(prep_list){
  mrna = prep_list$mrna
  protein = prep_list$protein
  
  # compute protein-level averages based on peptide-level log2 intensities
  protein = protein %>%
    dplyr::group_by(UNIPROT, ct) %>%
    dplyr::summarise(prot_av = mean(pep_av, na.rm = T))
  
  # compute mrna lfcs
  mrna = mrna %>%
    dplyr::group_by(UNIPROT, ct) %>%
    dplyr::summarise(mrna_av = mean(log2(mrna_av), na.rm = T))
  
  df = merge(mrna, protein)
  
  # draw plots
  g = ggplot() +
    geom_point(data = df, mapping = aes(x = mrna_av, y = prot_av)) +
    facet_wrap(~ct, ncol = 3) +
    xlab("Mean Log2 Transcript Count") +
    ylab("Mean Log2 Protein Intensity") +
    ggtitle("mRNA and Protein Observed Averages") +
    theme(text = element_text(size = 35),
          axis.text = element_text(size = 22),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# histogram of mrna averages
obs_hist_mrna = function(prep_list){
  mrna = prep_list$mrna
  # compute averages
  g_vline = mrna %>% ungroup() %>% dplyr::group_by(pop_mrna) %>% dplyr::summarise(obs_med = median(mrna_av, na.rm = T))
  
  # draw histogram
  g = ggplot(mrna, aes(x = mrna_av)) +
    geom_histogram(bins = 100) +
    geom_vline(data = g_vline, aes(xintercept = obs_med)) +
    ggtitle("Observed Transcript Level Averages") +
    facet_wrap(pop_mrna ~ ct) +
    xlab("Transcript Average Across Cells") +
    theme(text = element_text(size = 35),
          axis.text = element_text(size = 22),
          legend.key.size = unit(2, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# histgram of protein level averages
obs_hist_protein = function(prep_list){
  protein = prep_list$protein
  # compute averages
  g_vline = protein %>% ungroup() %>% dplyr::group_by(pop_protein) %>% dplyr::summarise(obs_med = median(pep_av, na.rm = T))
  
  # draw histogram
  g = ggplot(protein, aes(x = pep_av)) +
    geom_histogram(bins = 100) +
    geom_vline(data = g_vline, aes(xintercept = obs_med)) +
    ggtitle("Observed Peptide Level Averages") +
    facet_wrap(pop_protein ~ ct) +
    xlab("Peptide Level Average") +
    theme(text = element_text(size = 35),
          axis.text = element_text(size = 22),
          legend.key.size = unit(2, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# plot protein observed variances (across cell types or proteins)
plot_protein_obs_var = function(prep_list, across_celltypes = T, peptide_level = T){
  protein = prep_list$protein
  # compute across cell types variance if true (at peptide or protein levels)
  if(across_celltypes == TRUE){
    if(peptide_level == TRUE){
      df = protein %>%
        dplyr::group_by(pep, pop_protein) %>%
        dplyr::summarise(obs_var = var(pep_av, na.rm = T))
    }
    if(peptide_level == FALSE){
      df = protein %>%
        dplyr::group_by(UNIPROT, pop_protein) %>%
        dplyr::summarise(obs_var = var(pep_av, na.rm = T))
    }
    g_title = "Observed Variance Across Cell Types" # plot title
  }
  
  # compute across protein variances if false
  if(across_celltypes == FALSE){
    df = protein %>%
      dplyr::group_by(ct, pop_protein) %>%
      dplyr::summarise(obs_var = var(pep_av, na.rm = T))
    
    g_title = "Observed Variance Across Genes" # plot title
  }
  
  # vertical line at median
  g_vline = df %>% ungroup() %>% dplyr::group_by(pop_protein) %>% dplyr::summarise(var_med = median(obs_var, na.rm = T))
  
  # draw histogram
  g = ggplot(df, aes(x = obs_var)) +
    geom_histogram(bins = 100) +
    geom_vline(data = g_vline, aes(xintercept = var_med)) +
    ggtitle(g_title) +
    facet_wrap(vars(pop_protein)) +
    xlab("Observed Variance") +
    theme(text = element_text(size = 35),
          axis.text = element_text(size = 22),
          legend.key.size = unit(2, "cm"),
          legend.position = "bottom",
          panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# plot z transformed averages of mrna, protein and posterior intervals of significant genes
plot_zscore_means = function(mrna_z, protein_z, plex_z, gene_res,
                             uni, clusters,
                             mrna_pop_label, protein_pop_label,
                             mrna_color, protein_color,
                             pass_text = 40, title_add = NULL, sig = TRUE){
  # number of mrna and protein data sets
  n_pop_m = length(mrna_z)
  n_pop_p = length(protein_z)
  
  # identify symbols and genenames associated with gene
  gene_rec = AnnotationDbi::select(org.Hs.eg.db, keys = uni, keytype = "UNIPROT", columns = c("SYMBOL", "GENENAME"))
  
  # label data sets for each modality and formatting
  mrna_z = lapply(mrna_z, function(x) merge(x, gene_rec))
  mrna_z = lapply(mrna_z, function(x) mutate(x, df_lab = NULL))
  for(i in 1:length(mrna_z)){
    mrna_z[[i]] = mutate(mrna_z[[i]], measure = mrna_pop_label[i])
  }
  mrna_z = rlist::list.rbind(mrna_z)
  
  protein_z = lapply(protein_z, function(x) mutate(x, df_lab = NULL))
  for(i in 1:length(protein_z)){
    protein_z[[i]] = mutate(protein_z[[i]], measure = protein_pop_label[i])
  }
  protein_z = rlist::list.rbind(protein_z)
  
  # filter to selected gene for each modality, arrange by cell type
  mrna_z = mrna_z %>%
    dplyr::filter(UNIPROT == uni & ct %in% clusters) %>%
    dplyr::mutate(ct_factor = factor(ct, levels = clusters)) %>%
    arrange(ct_factor)
  
  protein_z = protein_z %>%
    dplyr::filter(UNIPROT == uni & ct %in% clusters) %>%
    dplyr::mutate(ct_factor = factor(ct, levels = clusters)) %>%
    arrange(ct_factor)
  
  # same thing for fitted parameter data
  gene_res = gene_res %>%
    merge(gene_rec) %>%
    filter(UNIPROT == uni & ct %in% clusters) %>%
    mutate(ct_factor = factor(ct, levels = clusters)) %>%
    arrange(ct_factor)
  
  # use full gene names for selected uniprot id
  gene_name_select = gene_rec %>% pull(SYMBOL) %>% unique()
  if(length(gene_name_select) >= 1){
    gene_name_select = gene_name_select %>% sample(1) # if there is more than one gene name present, randomly select one
  }
  ylab_select = gene_name_select # axis label
  
    # draw observed z-scrore plot
    g1 = ggplot() +
      geom_hline(yintercept = 0, linetype = 3) +
      geom_point(data = mrna_z, aes(x = ct_factor, y = mrna_zscore, color = measure), shape = 19, size = 20) +
      geom_path(data = mrna_z, aes(x = ct_factor, y = mrna_zscore, color = measure, group = pop_mrna, linetype = measure), size = 5) +
      geom_point(data = protein_z, aes(x = ct_factor, y = protein_zscore, color = measure), shape = 19, size = 20) +
      geom_path(data = protein_z, aes(x = ct_factor, y = protein_zscore, color = measure, group = pop_protein, linetype = measure), size = 5) +
      ylab("Log2 Fold Change") +
      xlab("") +
      scale_color_manual(guide = "none", values = c(mrna_color, mrna_color, protein_color, protein_color)) +
      scale_linetype_manual(name = "", values = c(1, 2, 2, 1)) +
      scale_size_manual(guide = "none") +
      guides(linetype = guide_legend(nrow = 2, override.aes = list(color = c(mrna_color, mrna_color, protein_color, protein_color),
                                                                   size = 40))) +
      theme(panel.background = element_rect(fill = 'white', color = "slategray4"),
            panel.grid.major = element_line(color = 'slategray2'),
            panel.grid.minor = element_line(color = 'slategray1'),
            panel.border = element_blank(),
            legend.text = element_text(margin = margin(r = 1, unit="inch"), size = pass_text),
            axis.text.y = element_text(size = pass_text + 20),
            axis.text.x = element_text(size = pass_text + 20, vjust = 0),
            legend.key = element_rect(colour = NA, fill = NA),
            legend.box.background = element_blank(),
            legend.position = "bottom")
    
    # draw posterior interval for ratio
    g2 = ggplot() +
      geom_hline(yintercept = 0, linetype = 3) +
      geom_point(data = gene_res, mapping = aes(x = ct_factor, y = r_av), color = "gray47", size = 15) +
      geom_linerange(data = gene_res, mapping = aes(x = ct_factor, ymin = r_lwr, ymax = r_upr), color = "gray47", linewidth = 15) +
      geom_path(data = gene_res, mapping = aes(x = ct_factor, y = r_av, group = UNIPROT)) +
      geom_text(x = 1.5, y = -1, aes(label = ylab_select), size = 40) +
      xlab("") +
      ylab("rPTR Interval") +
      theme(panel.background = element_rect(fill = 'white', color = "slategray4"),
            panel.grid.major = element_line(color = 'slategray2'),
            panel.grid.minor = element_line(color = 'slategray1'),
            legend.text = element_text(size = pass_text),
            panel.border = element_blank(),
            axis.text.x = element_blank())
    
    # if any cell types are significant, color significant interval 
    if(sig == TRUE){
      g2 = ggplot() +
        geom_hline(yintercept = 0, linetype = 3) +
        geom_point(data = gene_res, mapping = aes(x = ct_factor, y = r_av, color = significant), size = 15) +
        geom_linerange(data = gene_res, mapping = aes(x = ct_factor, ymin = r_lwr, ymax = r_upr, color = significant), linewidth = 15) +
        geom_path(data = gene_res, mapping = aes(x = ct_factor, y = r_av, group = UNIPROT)) +
        geom_text(x = 1.5, y = 1.5, aes(label = ylab_select), size = 45, color = "limegreen", fontface = "bold") +
        xlab("") +
        ylab("rPTR Interval") +
        scale_color_manual(values = c("purple", "gray47"), name = "", limits = c(TRUE, FALSE),
                           labels = c("Significant rPTR"), breaks = c(TRUE)) +
        scale_size_manual(guide = "none") +
        theme(panel.background = element_rect(fill = 'white', color = "slategray4"),
              panel.grid.major = element_line(color = 'slategray2'),
              panel.grid.minor = element_line(color = 'slategray1'),
              legend.position = c(0.19, 0.75),
              legend.text = element_text(size = pass_text),
              legend.background = element_blank(),
              legend.box.background = element_blank(),
              axis.text.y = element_text(size = pass_text + 20),
              panel.border = element_blank(),
              axis.text.x = element_blank())
    }
    
    if(!is.null(title_add)){
      g2 = g2 + ggtitle(title_add)
    }
  return(g2 / g1)
}

# plot observed z-transformed data for two genes (high and low reliability) without posterior rptr intervals
plot_zscore_means_double = function(mrna_z, protein_z, plex_z,
                             uni, organism = "human", clusters,
                             mrna_pop_label, protein_pop_label,
                             mrna_color, protein_color,
                             pass_text = 40, title_add = NULL, sig = TRUE){
  # number of data sets
  n_pop_m = length(mrna_z)
  n_pop_p = length(protein_z)
  
  # link UNIPROT id to symbol and genename
  gene_rec = AnnotationDbi::select(org.Hs.eg.db, keys = uni$UNIPROT, keytype = "UNIPROT", columns = c("SYMBOL", "GENENAME"))
  
  # label and format mrna
  mrna_z = lapply(mrna_z, function(x) merge(x, gene_rec))
  mrna_z = lapply(mrna_z, function(x) mutate(x, df_lab = NULL))
  for(i in 1:length(mrna_z)){
    mrna_z[[i]] = mutate(mrna_z[[i]], measure = mrna_pop_label[i])
  }
  mrna_z = rlist::list.rbind(mrna_z)
  
  # label and format protein
  protein_z = lapply(protein_z, function(x) mutate(x, df_lab = NULL))
  for(i in 1:length(protein_z)){
    protein_z[[i]] = mutate(protein_z[[i]], measure = protein_pop_label[i])
  }
  protein_z = rlist::list.rbind(protein_z)
  
  # filter to selected gene for each modality, arrange by cell type
  mrna_z = mrna_z %>%
    merge(uni) %>%
    dplyr::filter(ct %in% clusters) %>%
    dplyr::mutate(ct_factor = factor(ct, levels = clusters)) %>%
    arrange(ct_factor) 
  
  protein_z = protein_z %>%
    merge(uni) %>%
    dplyr::filter(ct %in% clusters) %>%
    dplyr::mutate(ct_factor = factor(ct, levels = clusters)) %>%
    arrange(ct_factor) 
  
  # use full gene names for selected uniprot id
  gene_name_select = gene_rec %>% pull(SYMBOL) %>% unique()
  ylab_select = gene_name_select # axis label
  name_set = mrna_z %>% dplyr::group_by(var_group) %>% dplyr::summarise(SYMBOL = unique(SYMBOL))
  
    # draw plot
    g1 = ggplot() +
      geom_hline(yintercept = 0, linetype = 3) +
      geom_point(data = mrna_z, aes(x = ct_factor, y = mrna_zscore, color = measure), shape = 19, size = 20, show.legend = F) +
      geom_path(data = mrna_z, aes(x = ct_factor, y = mrna_zscore, color = measure, group = pop_mrna, linetype = measure), size = 5, show.legend = T) +
      geom_point(data = protein_z, aes(x = ct_factor, y = protein_zscore, color = measure), shape = 19, size = 20, show.legend = F) +
      geom_path(data = protein_z, aes(x = ct_factor, y = protein_zscore, color = measure, group = pop_protein, linetype = measure), size = 5, show.legend = F) +
      geom_text(data = name_set, aes(x = 5, y = 1.5, label = SYMBOL), size = 50) +
      ylab("Log2 Fold Change") +
      xlab("") +
      facet_grid(rows = vars(var_group)) +
      scale_color_manual(guide = "none", values = c(mrna_color, mrna_color, protein_color, protein_color)) +
      scale_linetype_manual(name = "", values = c(1, 2, 2, 1)) +
      scale_size_manual(guide = "none") +
      guides(linetype = guide_legend(nrow = 2, override.aes = list(color =  c(mrna_color, mrna_color, protein_color, protein_color),
                                                                size = 80))) +
      theme(panel.background = element_rect(fill = 'white', color = "slategray4"),
            panel.grid.major = element_line(color = 'slategray2'),
            panel.grid.minor = element_line(color = 'slategray1'),
            legend.text = element_text(margin = margin(r = 0.5, unit="inch"), size = pass_text),
            legend.position = "bottom")
  return(g1)
}

# function to return list of dfs with significant genes that can be outputted as csv
# input tested gene information or set of posterior draws
output_significant_genes = function(gene_res = NULL, posterior_draws = NULL){
  if(is.null(gene_res)){
    gene_res = test_genes(posterior_draws)
  }
  fit_sig = filter(gene_res, significant == TRUE) %>% # filter to significant genes
    dplyr::select(UNIPROT, ct, fdr, r_av, r_lwr, r_upr, prot_av, mu_av) %>%
    arrange(fdr)
  colnames(fit_sig) = c("UNIPROT", "Cell Type", "Expected Proportion of False Discoveries",
                        "rPTR Average", "rPTR Lower", "rPTR Upper", "Protein Average", "mRNA Average")
  return(fit_sig)
}

# volcano plot of transformed fdr and posterior mean of rptr
plot_volcano = function(gene_res, sig_color = "purple", test_name = "Significant Nonzero", pass_text = 40, ct_select = NULL, ncol_facet = 3,
                        special_gene = NULL){
  gene_res = gene_res %>%
    dplyr::mutate(p_transform = -1*log10(fdr + 0.0001)) # transform fdr
  
  # possibly select only specific cell types 
  if(!is.null(ct_select)){
    gene_res = gene_res %>%
      dplyr::filter(ct %in% ct_select)
  }
  
  # summarize number of significant gene products in each cell type and across cell types
  ns_summary = gene_res %>% dplyr::group_by(ct) %>% dplyr::summarise(ns = sum(significant_fdr))
  ns_summary_across = gene_res %>% filter(significant_fdr == T) %>% pull(UNIPROT) %>% unique() %>% length()
  
  # plot transformed p-vlaue by r_av
  g = ggplot() +
    geom_point(data = gene_res, mapping = aes(x = r_av, y = p_transform, color = significant_fdr, alpha = significant_fdr, size = significant_fdr)) +
    geom_label(data = ns_summary, aes(x = 2, y = 1, label = ns), color = sig_color, size = 30) +
    facet_wrap(~factor(ct, levels = c("EC", "PTM", "LC", "SPG", "SPC", "St")), ncol = ncol_facet) +
    xlab("rPTR (Posterior Mean)") +
    ylab("-log10(FDR)") +
    ggtitle(paste("Total # Extreme Gene Products:", ns_summary_across)) +
    scale_color_manual(values = c("gray47", sig_color), guide = "none") +
    scale_size_manual(values = c(2, 12), guide = "none") +
    scale_alpha_manual(values = c(0.45, 1), guide = "none") +
    theme(panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  
  # if there's a "special gene" of interest, show that point in limegreen
  if(!is.null(special_gene)){
    gene_res_sub = gene_res %>% filter(UNIPROT == special_gene)
    g = ggplot() +
      geom_point(data = gene_res, mapping = aes(x = r_av, y = p_transform, color = significant_fdr, alpha = significant_fdr, size = significant_fdr)) +
      geom_point(data = gene_res_sub, mapping = aes(x = r_av, y = p_transform), alpha = 1, size = 20, color = "limegreen") +
      geom_text(data = ns_summary, aes(x = 2, y = 1, label = ns), color = sig_color, size = 30, fontface = "bold") +
      facet_wrap(~factor(ct, levels = c("EC", "PTM", "LC", "SPG", "SPC", "St")), ncol = ncol_facet) +
      xlab("rPTR (Posterior Mean)") +
      ylab("-log10(FDR)") +
      ggtitle(paste("Total # Significant Gene Products:", ns_summary_across)) +
      scale_x_continuous(breaks = c(-2.5, 0, 2.5), labels = c("-2.5", "0", "2.5")) +
      scale_color_manual(values = c("gray47", sig_color), guide = "none") +
      scale_size_manual(values = c(2, 12), guide = "none") +
      scale_alpha_manual(values = c(0.45, 1), guide = "none") +
      theme(panel.background = element_rect(fill = 'white', color = "slategray4"),
            panel.grid.major = element_line(color = 'slategray2'),
            panel.grid.minor = element_line(color = 'slategray1'))
  }
  return(g)
}

# plot mu + r by mu, with color corresponding to significant test results
# one of "none", "mrna" or "protein" for posterior_interval (mrna or protein gives 95 pct mu and mu + r interval, respectively)
plot_posterior_means = function(gene_res, posterior_interval = "none", facet_col = 3, sig_color = "purple", test_name = "Significant Nonzero", pass_text = 40, small = 3, big = 4.5){
  # draw point plot with reference diagonal, line at zero
  g = ggplot() +
    geom_point(data = gene_res, mapping = aes(x = mu_av, y = prot_av, color = significant_fdr, size = significant_fdr, alpha = significant_fdr)) +
    geom_abline(color = "red") +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_vline(xintercept = 0, linetype = 3) +
    facet_wrap(~ct, ncol = facet_col) +
    xlab("mRNA Posterior Mean") +
    ylab("Protein Posterior Mean") +
    scale_color_manual(values = c("gray47", sig_color), guide = "none", name = "Significant rPTR") +
    scale_size_manual(values = c(small, big), guide = "none") +
    scale_alpha_manual(values = c(0.45, 0.75), guide = "none") +
    theme(panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  
  # can optionally add error bars to x or y directions
  if(posterior_interval == "protein"){
    g = g + geom_errorbar(data = gene_res, aes(x = mu_av, y = prot_av, ymin = prot_lwr, ymax = prot_upr, color = significant))
  }
  
  if(posterior_interval == "mrna"){
    g = g + geom_errorbar(data = gene_res, aes(x = mu_av, y = prot_av, xmin = mu_lwr, xmax = mu_upr, color = significant))
  }
  return(g)
}

# returns point plot displaying observed and model-fit correlations
plot_correlation_comparison_med = function(cor_comparison, gene_res, pass_text = 40){
  # identify gene products with significance in at least once cell type
  gene_res = gene_res %>% dplyr::group_by(UNIPROT) %>% 
    dplyr::summarise(significant = sum(significant) >= 1) %>% ungroup()
  
  # compute posterior means of previously computed correlation comparison
  cor_info = cor_comparison %>%
    dplyr::group_by(UNIPROT) %>%
    dplyr::summarise(n_ct = round(median(n_ct, na.rm = T)), # use across-data set medians
                     param_cor = median(param_cor, na.rm = T),
                     pop_cor = median(pop_cor, na.m = T),
                     param_cor_lwr = median(param_cor_lwr, na.rm = T),
                     param_cor_upr = median(param_cor_upr, na.rm = T)) %>%
    ungroup() %>%
    merge(gene_res)
  
  # draw main plot
  g1 = ggplot() +
    geom_point(data = cor_info, aes(y = param_cor, x = pop_cor, color = significant, 
                                    size = significant, alpha = significant)) +
    geom_abline(color = "red") +
    geom_smooth(method = "lm", se = FALSE, color = "blue", alpha = 1) +
    geom_text(aes(label = paste("Significant rPTR", "(One or More Cell Types)",
                                sep = "\n")), x = 0.5, y = -0.75, size = 20) +
    xlab(paste("Across Clusters Correlation (Empirical)", " ", sep = "\n")) +
    ylab("Across Clusters Correlation (Model Fit)") +
    scale_x_continuous(limits=c(-1, 1), expand = c(0, 0)) +
    scale_y_continuous(limits=c(-1, 1), expand = c(0, 0)) +
    scale_color_manual(values = c("gray47", "purple")) +
    scale_alpha_manual(values = c(0.75, 1), guide = "none") +
    scale_size_manual(values = c(4, 6), guide = "none") +
    theme(panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'),
          plot.margin=unit(c(3, 3, 3, 3),"cm"),
          text = element_text(size = pass_text),
          axis.text = element_text(size = pass_text, vjust = 0.25, hjust = 0.25),
          legend.position = "none")
  
  # prepare data for marginal plot
  cor_info_piv = cor_comparison %>%
    merge(gene_res)
  g2_lab = cor_info_piv %>% dplyr::group_by(significant) %>% dplyr::summarise(pop_med = median(pop_cor, na.rm = T),
                                                                              lab = paste0("Median = ", round(pop_med, 2)),
                                                                              x = pop_med,
                                                                              y = 0.1)
  
  g3_lab = cor_info %>% dplyr::group_by(significant) %>% dplyr::summarise(param_med = median(param_cor, na.rm = T),
                                                                          y = 0.5, x = 3)
  
  # draw marginal density plots
  g2 = ggplot(data = cor_info_piv, mapping = aes(x = pop_cor, color = significant)) +
    geom_density(data = cor_info_piv, mapping = aes(x = pop_cor, color = significant, fill = significant), alpha = 0.1) +
    geom_vline(data = g2_lab, aes(xintercept = pop_med, color = significant)) +
    scale_x_continuous(limits=c(-1, 1), expand = c(0, 0), labels = NULL, breaks = NULL) +
    scale_color_manual(values = c("gray47", "purple"), guide = "none") +
    scale_fill_manual(values = c("gray47", "purple"), guide = "none") +
    # xlab("") +
    ylab("") +
    labs(x = NULL) +
    theme_void()
  
  # marginal plot for fitted parameter
  g3 = ggplot(data = cor_info, mapping = aes(y = param_cor, color = significant)) +
    geom_density(data = cor_info, mapping = aes(y = param_cor, color = significant, fill = significant), alpha = 0.1) +
    geom_hline(data = g3_lab, aes(yintercept = param_med, color = significant)) +
    scale_y_continuous(limits=c(-1, 1), expand = c(0, 0)) +
    scale_color_manual(values = c("gray47", "purple"), guide = "none", name = "Significant rPTR") +
    scale_fill_manual(values = c("gray47", "purple"), guide = "none") +
    xlab("") +
    ylab("") +
    theme(panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          panel.margin = unit(c(0, 0, 0, 0), "cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          legend.position = "none",
          axis.ticks.length = unit(0, "cm"),
          axis.ticks.margin = unit(0, "cm"),
          legend.margin = unit(0, "cm"))
  
  design = "
  AAA#
  BBBC
  BBBC
  BBBC
  "
  g = g2 + g1 + g3 + plot_layout(design = design)
  return(g)
}

# table of posterior predictive statistics
plot_ppc_stats = function(mrna_ppc_variance, protein_ppc_variance,
                          mrna_ppc_z, protein_ppc_z, mrna_obs_z, protein_obs_z,
                          text_pass = 20){
  
  # link data sets to measurement techniques 
  mrna_df_rec = data.frame(pop_label = paste0("mRNA Pop", 1:5),
                           dataset = c(rep("10x", 2), rep("Drop-Seq", 3)))
  
  protein_df_rec = data.frame(pop_label = paste0("Protein Pop", 1:5),
                              dataset = c(rep("SCoPE2", 4), "plexDIA"))
  
  # labels related to pp variance
  mrna_ppc_variance = mrna_ppc_variance %>%
    dplyr::mutate(pop_label = paste("mRNA", pop_mrna),
                  modality = "mRNA") %>%
    dplyr::select(-pop_mrna) %>%
    merge(mrna_df_rec)
  
  protein_ppc_variance = protein_ppc_variance %>%
    dplyr::mutate(pop_label = paste("Protein", pop_protein),
                  modality = "Protein") %>%
    dplyr::select(-pop_protein) %>%
    merge(protein_df_rec)
  
  ppc_var = rbind(mrna_ppc_variance, protein_ppc_variance) # combine pp variance for each modality 
  variance_lab = ppc_var %>% dplyr::group_by(dataset) %>%
    # compute coverage for each data set
    dplyr::summarise(cover = round(mean(obs_var >= ppc_var_lwr & obs_var <= ppc_var_upr), 2),
                     cover_lab = paste("Coverage =", round(cover, 2))) %>%
    ungroup() %>%
    dplyr::mutate(Variance = ifelse(cover == 1, 0.99, cover)) %>%
    dplyr::select(dataset, Variance) %>%
    pivot_longer(cols = c(Variance), names_to = "Posterior_Predictive") %>%
    pivot_wider(names_from = dataset, values_from = value) # do some formatting
  
  # prepare posterior predictive fold changes
  mrna_ppc_z = mrna_ppc_z %>%
    dplyr::mutate(pop_label = paste("mRNA", pop_mrna),
                  ppc_z = mrna_ppc_zscore,
                  ppc_lwr = mrna_ppc_lwr,
                  ppc_upr = mrna_ppc_upr,
                  ppc_med = mrna_ppc_med,
                  modality = "mRNA",
                  ct = as.character(ct)) %>%
    dplyr::select(-c(mrna_ppc_zscore, mrna_ppc_lwr, mrna_ppc_upr, mrna_ppc_med, pop_mrna)) %>%
    merge(mrna_df_rec)
  
  protein_ppc_z = protein_ppc_z %>%
    dplyr::mutate(pop_label = paste("Protein", pop_protein),
                  ppc_z = protein_ppc_zscore,
                  ppc_lwr = protein_ppc_lwr,
                  ppc_upr = protein_ppc_upr,
                  ppc_med = protein_ppc_med,
                  modality = "Protein",
                  ct = as.character(ct)) %>%
    dplyr::select(-c(protein_ppc_zscore, protein_ppc_lwr, protein_ppc_upr, protein_ppc_med, pop_protein)) %>%
    merge(protein_df_rec)
  
  # prepare observed fold changes
  mrna_obs_z = mrna_obs_z %>%
    dplyr::mutate(obs_z = mrna_zscore,
                  ct = as.character(ct))%>%
    dplyr::select(-mrna_zscore) %>% merge(mrna_df_rec)
  
  protein_obs_z = protein_obs_z %>%
    dplyr::mutate(obs_z = protein_zscore,
                  ct = as.character(ct)) %>%
    dplyr::select(-protein_zscore) %>% merge(protein_df_rec)
  
  coverage_df = rbind(mrna_ppc_z, protein_ppc_z) %>%
    merge(rbind(mrna_obs_z, protein_obs_z))
  
  # compute coverage of fold changes
  z_lab = coverage_df %>% dplyr::group_by(dataset) %>%
    dplyr::summarise(cover = round(mean(obs_z >= ppc_lwr & obs_z <= ppc_upr, na.rm = T), 2)) %>%
    ungroup() %>%
    dplyr::mutate("Fold Change" = as.character(ifelse(cover == 1, 0.99, cover))) %>%
    dplyr::select(dataset, "Fold Change") %>%
    pivot_longer(cols = c("Fold Change"), names_to = "Posterior_Predictive") %>%
    pivot_wider(names_from = dataset, values_from = value)
  
  # combine variance and fold change posterior predictive statistics
  coverage_lab = rbind(z_lab, variance_lab)
  colnames(coverage_lab)[1] = "Posterior Predictive"
  
  # draw table
  g = kableExtra::kbl(coverage_lab, align = "c") %>%
    kableExtra::kable_classic(font_size = text_pass, full_width = F) %>%
    kableExtra::add_header_above(c(" ", "mRNA" = 2, "Protein" = 2))
  return(g)
}

