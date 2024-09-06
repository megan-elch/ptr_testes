
# compute protein consistency at peptide level
compute_protein_consistency = function(df, df_info, protein_label = "prot", peptide_label = "pep", seed_id = 100){
  df[["UNIPROT"]] = df[[protein_label]]
  df[["pep"]] = df[[peptide_label]]
  cells = intersect(df_info$id, colnames(df))
  
  set.seed(seed_id) 
  # assign peptides to group a
  df_group_a = df %>% 
    dplyr::group_by(UNIPROT) %>%
    dplyr::sample_n(size = floor(length(unique(pep))/2)) %>%
    ungroup()
  
  # all remaining peptides assigned to group b
  peptides_a = unique(df_group_a$pep)
  peptides_b = setdiff(unique(df$pep), peptides_a)
  
  # compute averages for group a
  df_group_a = df_group_a %>%
    pivot_longer(cols = cells, names_to = "id") %>%
    merge(df_info) %>%
    dplyr::group_by(UNIPROT, ct) %>%
    dplyr::summarise(prot_a = mean(value, na.rm = T)) %>%
    ungroup()
  
  # compute averages for group b and merge with group a averages, correlate
  reliabilities = df %>%
    filter(pep %in% peptides_b) %>%
    pivot_longer(cols = cells, names_to = "id") %>%
    merge(df_info) %>%
    dplyr::group_by(UNIPROT, ct) %>%
    dplyr::summarise(prot_b = mean(value, na.rm = T)) %>%
    ungroup() %>%
    merge(df_group_a) %>%
    dplyr::group_by(UNIPROT) %>%
    dplyr::summarise(reliability = cor(prot_a, prot_b, use = "pairwise.complete.obs"))
  return(reliabilities)
}
