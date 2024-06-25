#########################################################################################################################################
#                                                        FUNCTIONS RELATED TO GO TESTING 
#########################################################################################################################################
# function to run GO test. Input posterior samples, and GO_info, GO_info2, GO_tab contain GO mapping from GO_prep.RData
run_go_test_hs = function(posterior_draws, GO_info, GO_info2, GO_tab, min_size = 1, max_size = 150, threshold = 0.01){
  UNIPROT_in = unique(posterior_draws$UNIPROT) 
  # combine previous information
  GO_info = GO_info %>%
    merge(GO_tab, by.x = "GO", by.y = "GOID") %>%
    merge(GO_info2)
  rm(list = c("GO_tab", "GO_info2"))
  
  # filter go groups by total number of genes and number of genes present in data
  GO_info = GO_info %>%
    filter(total_terms < max_size, nd > min_size)
  
  # combine go information with posterior draws for later testing
  posterior_draws_go = posterior_draws %>% 
    base::merge(GO_info, by = "UNIPROT") %>% 
    na.omit()
  
  # record size information           
  size_info = posterior_draws_go %>%
    dplyr::group_by(GO) %>%
    dplyr::summarise(total_terms = unique(total_terms),
                     nd = unique(nd))
  
  test_res = posterior_draws_go %>%
    dplyr::group_by(.iteration, .chain, GO) %>% # record average across genes not associated with each group
    dplyr::mutate(r_group = mean(r, na.rm = T)) %>%
    ungroup() %>%
    dplyr::mutate(r = r - r_group) %>%
    dplyr::group_by(.iteration, .chain, GO, ct) %>%
    dplyr::summarise(r_mean = mean(r, na.rm = T), # compute sample average across genes for each group, posterior draw
                     mu_mean = mean(mu, na.rm = T),
                     prot_mean = mean(prot, na.rm = T),
                     TERM = unique(TERM)) %>%
    ungroup() %>%
    dplyr::group_by(GO, ct) %>%
    dplyr::summarise(TERM = unique(TERM), 
                     r_av = mean(r_mean, na.rm = T), # average r in group
                     mu_av = mean(mu_mean, na.rm = T), # average mu in group
                     prot_av = mean(prot_mean, na.rm = T), # average mu + r in group
                     r_lwr = quantile(r_mean, 0.025, na.rm = T), # 95 pct interval r
                     r_med = median(r_mean, na.rm = T),
                     r_upr = quantile(r_mean, 0.975, na.rm = T),
                     p_r = ifelse(mean(r_mean < 0) <= 1 - mean(r_mean < 0), mean(r_mean < 0), 1 - mean(r_mean < 0)),
                     significant = !(r_lwr < 0 & r_upr > 0)) %>%
    ungroup()
  list(posterior_draws_go = posterior_draws_go, test_res = test_res)
}

# compute gene level information for significant go groups
compute_ticks_draws = function(posterior_draws_go, test_res, celltype){
  posterior_draws_go = posterior_draws_go %>%
    dplyr::group_by(.iteration, .chain, GO) %>% # record average across genes not associated with each group
    dplyr::mutate(r_group = mean(r, na.rm = T)) %>%
    ungroup() %>%
    dplyr::mutate(param = r - r_group) %>%
    dplyr::filter(ct == celltype)
  
  go = test_res %>%
    dplyr::group_by(ct) %>%
    arrange(p_r, .by_group = T) %>% # compute expected proportion of false discoveries for r, mu, mu + r
    dplyr::mutate(fdr = cummean(p_r[order(p_r)])) %>%
    ungroup() %>%
    dplyr::filter(significant == T) %>%
    dplyr::select(ct, GO, fdr, r_lwr, r_upr) %>%
    distinct(.keep_all = TRUE)
  
  # compute gene-level summaries for members in each GO group
  GO_dat = posterior_draws_go %>%
    merge(go) %>%
    dplyr::group_by(UNIPROT) %>%
    dplyr::summarise(param_mean = mean(param),
                     GO = unique(GO),
                     TERM = unique(TERM),
                     ct = celltype,
                     prt = mean(fdr),
                     r_lwr = mean(r_lwr),
                     r_upr = mean(r_upr)) %>%
    ungroup() %>%
    dplyr::group_by(GO) %>%
    dplyr::mutate(go_mean = mean(param_mean)) %>%
    ungroup()
  return(GO_dat)
}

# output list of data frames with summary of significant results
output_significant_groups = function(test_res, gene_ref){
  # get number of observed genes associated with each group
  gene_ref2 = gene_ref %>% dplyr::group_by(GO) %>% dplyr::summarise(nd = length(unique(UNIPROT)))
  # keep significant test results
  test_res = test_res %>% filter(significant_fdr == TRUE)
  go_sig = test_res %>% pull(GO) %>% unique()
  # get number of genes associated with each group (may not be observed)
  GO_rec = AnnotationDbi::select(org.Hs.eg.db, keys = go_sig, keytype = "GO", columns = c("UNIPROT")) %>%
    ungroup() %>%
    dplyr::group_by(GO) %>%
    dplyr::summarise(total_terms = length(unique(UNIPROT))) %>%
    merge(gene_ref2)
  
  go_info = test_res %>% # select relevant summary columns
    dplyr::select(GO, ct, TERM,
                  fdr, r_av, mu_av, prot_av,
                  significant_fdr, r_lwr, r_upr) %>%
    merge(GO_rec) %>%
    distinct(.keep_all = TRUE) %>%
    arrange(fdr)
  
  colnames(go_info) = c("GO", "CT", "TERM",
                        "Expected Proportion of False Discoveries", "Ratio GO", "mRNA GO", "Protein GO",
                        "Significant Ratio", "Ratio Lower", "Ratio Upper", "Total Genes", "Genes in Data")
  # recover observed genes in each group for output file
  gene_ref = gene_ref %>% dplyr::group_by(GO) %>% dplyr::summarise(UNIPROT = paste(UNIPROT, collapse = ", "))
  go_info = merge(go_info, gene_ref)
  return(go_info)
}

# return ggplot with heatmap showing sample average of posterior means for each group
plot_go_heatmap = function(test_res, go_select = NULL, TERM_abbrev = NULL, ct_select = "PTM"){
  if(is.null(go_select)){
    go_select = test_res %>% filter(significant == TRUE) %>% pull(GO) %>% sample(10) # unless list of groups are provided, random sample
  }
  test_res = test_res %>%
    filter(GO %in% go_select) # filter to select groups
  
  # arrange go groups in ascending order
  go_order = test_res %>% filter(ct == ct_select) %>%
    ungroup() %>%
    mutate(GO_fct = factor(GO, levels = GO[order(-1*r_av)])) %>%
    dplyr::select(GO, GO_fct)
  
  # if any terms are abbreviated, merge abbreviated names
  if(!is.null(TERM_abbrev)){
    test_res = test_res %>% merge(TERM_abbrev) %>%
      dplyr::mutate(TERM = TERM_abbrev)
  }
  
    # merge ordering and order cell types
    test_res = test_res %>% merge(go_order, by = "GO") %>%
    dplyr::mutate(param_mean = r_av) %>% 
    mutate(ct = factor(ct, levels = c("EC", "PTM", "LC", "SPG", "SPC", "St"))) %>% 
    arrange(GO_fct)
  
  # generate heatmap
  g = ggplot(data = test_res, mapping = aes(x = ct, y = forcats::fct_inorder(TERM), fill = param_mean)) +
    geom_tile() +
    ylab("") +
    ggtitle("Average rPTR (Posterior Mean)") +
    scale_y_discrete(labels = function(y) str_wrap(y, width = 15)) +
    scale_fill_gradientn(colors = c("darkblue", "blue", "white", "red", "darkred"),
                         limits = c(-1, 1),
                         oob = scales::squish,
                         na.value = NA, name = "rPTR", labels = c(-1, "", 0, "", 1)) +
    theme(panel.background = element_rect(fill = 'white', color = "white"),
          panel.grid.major = element_line(color = 'white'),
          panel.grid.minor = element_line(color = 'white'),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          panel.margin = unit(c(0, 0, 0, 0), "cm"))
  return(g)
}

# return point plot of posterior means
plot_go_means = function(test_res){
  # point plot of posterior means of mu and mu + R for each go group
  g = ggplot(test_res, aes(x = mu_av, y = prot_av)) +
    geom_point(aes(color = significant, size = significant, alpha = significant), size = 4) +
    facet_wrap(vars(ct), ncol = 3) +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_vline(xintercept = 0, linetype = 3) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    scale_color_manual(guide = "none", values = c("gray47", "purple")) +
    scale_size_manual(values = c(3, 4.5), guide = "none") +
    scale_alpha_manual(values = c(0.45, 0.75), guide = "none") +
    xlab("mRNA Posterior Mean (Among GO Groups)") +
    ylab("Protein Posterior Mean (Among GO Groups)")  +
    theme(panel.background = element_rect(fill = 'white', color = "slategray4"),
          panel.grid.major = element_line(color = 'slategray2'),
          panel.grid.minor = element_line(color = 'slategray1'))
  return(g)
}

# draw tick mark plot displaying posterior mean for genes associated with selected significant go groups
plot_ticks = function(GO_dat, gene_res,
                       n = 10, seed_id = 15,
                       GO_list = NULL, celltype, celltype_title, TERM_abbrev = NULL,
                       pass_text = 40, go_order = NULL){
  
  plot_param = "r" # use ratio as plotting parameter
  x_title = "rPTR (Posterior Mean)" # axis label
  
  # if a list of groups are not supplied, randomly select groups with some filtering
  if(is.null(GO_list)){
    set.seed(seed_id)
    GO_list = GO_dat %>%
      dplyr::group_by(GO) %>%
      dplyr::summarise(go_av = mean(param_mean, na.rm = T),
                       ng = length(unique(UNIPROT))) %>%
      ungroup() %>%
      dplyr::mutate(ng_group = ntile(ng, 2)) %>% # bin number of genes observed for each group into two
      filter(ng_group == 2) %>% # keep groups in upper half of observed genes
      dplyr::mutate(av_group = ntile(go_av, 4)) %>% # bin group-level average into 4 groups 
      dplyr::group_by(av_group) %>%
      dplyr::sample_n(size = 2) %>% # select two groups from each bin
      ungroup() %>%
      pull(GO)
  }
  
  # filter go groups to those selected, put in descending order
  GO_dat = GO_dat %>%
    dplyr::filter(GO %in% GO_list) %>%
    dplyr::mutate(TERM_factor = fct_reorder(TERM, go_mean, .desc = TRUE)) %>%
    arrange(TERM_factor)
  
  # if abbreviated terms are given, merge and order accordingly
  if(!is.null(TERM_abbrev)){
    GO_dat = GO_dat %>% merge(TERM_abbrev) %>%
      dplyr::mutate(TERM_factor = fct_reorder(TERM_abbrev, go_mean, .desc = TRUE))
  }
  
  # use min and max across genes associated with all groups to auto-set axis limits
  min_check = min(GO_dat$param_mean[is.finite(GO_dat$param_mean)], na.rm = T)
  max_check = max(GO_dat$param_mean[is.finite(GO_dat$param_mean)], na.rm = T)
  lim_val = max(abs(min_check), abs(max_check))
  
  # draw tick mark plot
  bb =  ggplot() +
    geom_point(data = GO_dat, mapping = aes(x = param_mean, y = fct_inorder(TERM_factor), color = param_mean, group = UNIPROT),
               size = 80, stroke = 80, shape="|", alpha = 1) +
    geom_point(data = GO_dat, mapping = aes(x = go_mean, y = fct_inorder(TERM_factor), color = go_mean),
               size = 30, stroke = 10, shape = 4) +
    geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
    geom_text(x = 0, y = length(GO_list) + 0.75, aes(label = celltype_title), size = 40) +
    xlab(paste(x_title)) +
    ylab("") +
    xlim(-1*(lim_val), (lim_val)) +
    scale_color_gradientn(colors = c("darkblue", "blue", "white", "red", "darkred"),
                          limits = c(-1, 1),
                          oob = scales::squish,
                          na.value = NA, name = "", guide = "none") +
    scale_y_discrete(labels = function(y) str_wrap(y, width = 15)) +
    theme(panel.background = element_rect(fill = 'white', color = "white"),
          plot.margin = unit(c(2, 2, 2, 2), "cm"),
          panel.margin = unit(c(2, 2, 2, 2), "cm"),
          legend.position = "bottom",
          axis.text.y = element_blank(), # angle = 90, hjust = 0.5, vjust = 0.5),
          axis.text.x = element_text(size = pass_text),
          axis.title = element_text(size = pass_text),
          panel.grid.major = element_line(color = 'slategray1'),
          panel.grid.minor = element_line(color = 'white'))
  
  # marginal table displaying posterior probability associated with each group
  tbl_df = gene_res %>%
    ungroup() %>%
    filter(ct == celltype) %>%
    dplyr::group_by(TERM) %>%
    dplyr::summarise(Posterior_Probability = unique(fdr)) %>%
    dplyr::mutate(Posterior_Probability = ifelse(round(Posterior_Probability, 5) == 0, 0.00001, Posterior_Probability)) %>%
    merge(select(GO_dat, TERM, TERM_factor), all.y = T) %>%
    distinct() %>%
    dplyr::select(Posterior_Probability, TERM_factor) %>%
    # formatC(unique(prt), format = "e", digits = 2)) %>%
    arrange(TERM_factor)
  
  # plot table 
  tbl = ggplot(tbl_df, aes(y = fct_inorder(TERM_factor), x = 1)) +
    geom_text(mapping = aes(x = 1, y = fct_inorder(TERM_factor), 
                            label = as.character(vaxedemic::scientific_10x(Posterior_Probability, digits = 0))), 
              size = 50, parse = T) +
    xlab("") +
    ylab("") +
    theme(panel.background = element_rect(fill = 'white', color = "white"),
          plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          panel.grid.major = element_line(color = 'white'),
          panel.grid.minor = element_line(color = 'white'),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          panel.margin = unit(c(0, 0, 0, 0), "cm"))
  
  # draw title panel for table
  tbl_fill = ggplot() +
    geom_text(mapping = aes(x = 1, y = 1, label = "FDR"), size = 40) +
    xlab("") +
    ylab("") +
    theme(panel.background = element_rect(fill = 'white', color = "white"),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          panel.grid.major = element_line(color = 'white'),
          panel.grid.minor = element_line(color = 'white'),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          panel.margin = unit(c(0, 0, 0, 0), "cm"))
  
  # marignal density plot across all genes
  bd <- ggplot(gene_res) +
    geom_density(aes(x = r_av), fill = "gray47", alpha = 0.1) +
    xlab("") +
    ylab("") +
    xlim(-1*(lim_val), (lim_val)) +
    theme(panel.background = element_rect(fill = 'white', color = "white"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major = element_line(color = 'white'),
          panel.grid.minor = element_line(color = 'white'),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          panel.margin = unit(c(0, 0, 0, 0), "cm"))
  
  b = (((bd / bb) + plot_layout(height = c(1, 4))) | ((tbl_fill / tbl) + plot_layout(height = c(1, 4)))) + plot_layout(width = c(4, 1))
  
  return(b)
}

