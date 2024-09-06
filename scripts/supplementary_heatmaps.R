library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
library(colorspace)
library(viridis)
library(ggridges)
library(ggvenn)
library(ggpubr)
library(ggrepel)
library(tidyverse)
library(patchwork)
library(ggbeeswarm)
source("scripts/plot_functions.R")
source("scripts/go_functions.R")
clusters = c("EC", "PTM", "LC", "SPG", "SPC", "St")

mrna_color = "#2166AC"
protein_color = "#B2182B"
multi_colors = "slategray"
sig_colors = "magenta3"
text_size = 80

# format cell type titles
celltype_title = c("Endothelial Cells", "Peritubular Myoid Cells", "Leydig Cells", "Spermatagonia Cells", "Spermatocyte Cells", "Spermatid Cells")

####################################################################################################################################### 
#                                       Heatmaps, Tick Mark Plots for mRNA, Protein GO-level test results     
####################################################################################################################################### 

# load gene level test results
load("model_output/gene_res.RData")
gene_res = gene_res %>% 
  dplyr::group_by(ct) %>%
  arrange(p_r_scaled, .by_group = T) %>%
  dplyr::mutate(fdr = cummean(p_r_scaled[order(p_r_scaled)]),
                significant_fdr = fdr < 0.01,
                significant = significant_fdr)

# load go level test results
load("model_output/go_test_av.RData")

# map go terms to uniprot ids
gene_ref = AnnotationDbi::select(org.Hs.eg.db, keys = unique(test_res$GO), keytype = "GO", columns = c("UNIPROT")) %>%
  filter(UNIPROT %in% gene_res$UNIPROT) %>%
  dplyr::select(GO, UNIPROT)

# selected terms for heatmap and tick mark plot
TERMS = c("chaperone-mediated protein complex assembly", "histone methyltransferase activity (H3-K36 specific)",
          "glycolytic process", "citrate metabolic process", "fatty acid beta-oxidation using acyl-CoA dehydrogenase",
          "NADH dehydrogenase (ubiquinone) activity", "cytoplasmic translation", "double-strand break repair via break-induced replication")
GO_heatmap = test_res %>% filter(TERM %in% TERMS) %>% pull(GO) %>% unique()

# load complex test results
load("model_output/complex_test_av.RData")
# data frame mapping protein complexes to UNIPROT IDs
complex_uni =  read_delim("processed_data/coreComplexes.txt", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
complex_uni = complex_uni %>% 
  dplyr::group_by(ComplexName) %>% 
  dplyr::reframe(UNIPROT = unlist(strsplit(`subunits(UniProt IDs)`, ";", fixed = T)))

# select complexes of interest
complex_labels = c("Emerin architectural complex", "PA28gamma-20S proteasome", "Large Drosha complex", "HRD1 complex",
                   "TRBP containing complex (DICER, RPL7A, EIF6, MOV10 and subunits of the 60S ribosomal particle)", 
                   "Intraflagellar transport complex B", "Kinase maturation complex 1", "MCM complex")

####################################################################################################################################### 
#                                                     Summary Info: mRNA GO-level test results     
####################################################################################################################################### 

# fdr table of complex test 
complex_fdr = test_res_complex %>%
  dplyr::group_by(ct) %>%
  arrange(p_mu, .by_group = T) %>% # compute expected proportion of false discoveries for r, mu, mu + r
  dplyr::mutate(fdr = cummean(p_mu[order(p_mu)]), # expected proportion
                significant_fdr = fdr <= 0.05) %>%
  ungroup() %>%
  filter(ComplexName %in% complex_labels) %>%
  dplyr::select(ComplexName, fdr, ct)

# collect summary info for mrna complex heatmap and tick mark plot
complex_heat_mrna = test_res_complex %>%
  filter(ComplexName %in% complex_labels) %>%
  mutate(ComplexName = ifelse(ComplexName == "LAS1L-PELP1-TEX10-WDR18-NOL9-SENP3  complex",
                              "LAS1L-PELP1-TEX10- WDR18-NOL9-SENP3  complex",
                              ifelse(ComplexName == "Retromer complex (SNX1, SNX2, VPS35, VPS29, VPS26A)",
                                     "Retromer complex",
                                     ifelse(ComplexName == "Oligosaccharyltransferase complex (Stt3A variant)",
                                            "Oligosaccharyl transferase complex (Stt3A variant)", ComplexName)))) %>%
  dplyr::group_by(ComplexName) %>%
  dplyr::transmute(ct, mu_av, r_av_st = r_av[ct == "St"]) %>%
  ungroup() %>%
  mutate(ComplexName_factor = forcats::fct_reorder(ComplexName, r_av_st, .desc = T),
         type = "Protein Complex") %>%
  merge(complex_fdr) %>%
  arrange(ComplexName_factor)

# fdr table of go test
test_res = test_res %>%
  dplyr::group_by(ct) %>%
  arrange(p_mu, .by_group = T) %>% # compute expected proportion of false discoveries for r, mu, mu + r
  dplyr::mutate(fdr = cummean(p_mu[order(p_mu)]), # expected proportion
                significant_fdr = fdr <= 0.05)

# collect summary info for go heatmap and tick mark plot
go_heat_mrna = test_res %>%
  filter(GO %in% GO_heatmap) %>%
  dplyr::transmute(ComplexName = TERM, mu_av = mu_av, ct, fdr, r_av = r_av) %>%
  dplyr::group_by(ComplexName) %>%
  dplyr::mutate(r_av_st = r_av[ct == "St"]) %>%
  ungroup() %>%
  mutate(ComplexName_factor = forcats::fct_reorder(ComplexName, r_av_st, .desc = T),
         type = "GO Group") %>%
  arrange(ComplexName_factor) %>%
  dplyr::select(-r_av)

# combine complex and go level results
r_av_all = rbind(complex_heat_mrna, go_heat_mrna) %>%
  dplyr::mutate(ct_factor = factor(ct, levels = clusters))
complex_order = r_av_all %>%
  dplyr::group_by(ComplexName_factor) %>%
  dplyr::summarise(ComplexName = unique(ComplexName))

# mrna summary tick mark info for protein complexes st
complex_st_tick_mrna = gene_res %>%
  merge(complex_uni) %>%
  filter(ComplexName %in% complex_labels & ct == "St") %>%
  dplyr::group_by(ComplexName) %>%
  dplyr::mutate(complex_mean = mean(mu_av, na.rm = T)) %>%
  ungroup() %>%
  transmute(UNIPROT, ct, mu_av, ComplexName, complex_mean, type = "Protein Complex")

# mrna summary tick mark info for go groups st
go_st_tick_mrna = gene_ref %>%
  filter(GO %in% GO_heatmap) %>%
  merge(transmute(test_res, ct, GO, ComplexName = TERM, complex_mean = mu_av)) %>%
  merge(transmute(gene_res, ct, UNIPROT, ct, mu_av, type = "GO Group")) %>%
  filter(ct == clusters[6]) %>%
  dplyr::select(-GO)

# marginal table displaying posterior probability associated with each group
tbl_df = r_av_all %>%
  ungroup() %>%
  filter(ct == "St") %>%
  dplyr::group_by(ComplexName_factor, type) %>%
  dplyr::summarise(Posterior_Probability = unique(fdr)) %>%
  dplyr::mutate(Posterior_Probability = ifelse(round(Posterior_Probability, 5) == 0, 0.00001, Posterior_Probability)) %>%
  distinct() %>%
  dplyr::select(Posterior_Probability, ComplexName_factor, type) %>%
  arrange(ComplexName_factor, .by_group = T)

# separate complexes and go groups
complex_fdr_st_mrna = tbl_df %>% filter(type == "Protein Complex")
go_fdr_st_mrna = tbl_df %>% filter(type == "GO Group")

complex_spc_tick_mrna = gene_res %>%
  merge(complex_uni) %>%
  filter(ComplexName %in% complex_labels & ct == "SPC") %>%
  dplyr::group_by(ComplexName) %>%
  dplyr::mutate(complex_mean = mean(mu_av, na.rm = T)) %>%
  ungroup() %>%
  transmute(UNIPROT, ct, mu_av, ComplexName, complex_mean, type = "Protein Complex")

# plot tick marks: posterior mean of rptr for genes in each group
go_spc_tick_mrna = gene_ref %>%
  filter(GO %in% GO_heatmap) %>%
  merge(transmute(test_res, ct, GO, ComplexName = TERM, complex_mean = mu_av)) %>%
  merge(transmute(gene_res, ct, UNIPROT, ct, mu_av, type = "GO Group")) %>%
  filter(ct == clusters[5]) %>%
  dplyr::select(-GO)

####################################################################################################################################### 
#                                                     Summary Info: Protein GO-level test results     
####################################################################################################################################### 

# marginal table displaying posterior probability associated with each group
tbl_df = r_av_all %>%
  ungroup() %>%
  filter(ct == "SPC") %>%
  dplyr::group_by(ComplexName_factor, type) %>%
  dplyr::summarise(Posterior_Probability = unique(fdr)) %>%
  dplyr::mutate(Posterior_Probability = ifelse(round(Posterior_Probability, 5) == 0, 0.00001, Posterior_Probability)) %>%
  distinct() %>%
  dplyr::select(Posterior_Probability, ComplexName_factor, type) %>%
  arrange(ComplexName_factor, .by_group = T)

# separate complexes and go groups
complex_fdr_spc_mrna = tbl_df %>% filter(type == "Protein Complex")
go_fdr_spc_mrna = tbl_df %>% filter(type == "GO Group")

# selected terms for heatmap and tick mark plot
# format cell type title
complex_fdr = test_res_complex %>%
  dplyr::group_by(ct) %>%
  arrange(p_prot, .by_group = T) %>% # compute expected proportion of false discoveries for r, mu, mu + r
  dplyr::mutate(fdr = cummean(p_prot[order(p_prot)]), # expected proportion
                significant_fdr = fdr <= 0.05) %>%
  ungroup() %>%
  filter(ComplexName %in% complex_labels) %>%
  dplyr::select(ComplexName, fdr, ct)

# collect summary info for complex heatmap and tick mark plot
complex_heat_protein = test_res_complex %>%
  filter(ComplexName %in% complex_labels) %>%
  dplyr::group_by(ComplexName) %>%
  dplyr::transmute(ct, prot_av, r_av_st = r_av[ct == "St"]) %>%
  ungroup() %>%
  mutate(ComplexName_factor = forcats::fct_reorder(ComplexName, r_av_st, .desc = T),
         type = "Protein Complex") %>%
  merge(complex_fdr) %>%
  arrange(ComplexName_factor)

# go test results summary protein info 
test_res = test_res %>%
  dplyr::group_by(ct) %>%
  arrange(p_prot, .by_group = T) %>% # compute expected proportion of false discoveries for r, mu, mu + r
  dplyr::mutate(fdr = cummean(p_prot[order(p_prot)]), # expected proportion
                significant_fdr = fdr <= 0.05)

# collect summary info for go heatmap and tick mark plot
go_heat_protein = test_res %>%
  filter(GO %in% GO_heatmap) %>%
  dplyr::transmute(ComplexName = TERM, prot_av = prot_av, ct, fdr, r_av = r_av) %>%
  dplyr::group_by(ComplexName) %>%
  dplyr::mutate(r_av_st = r_av[ct == "St"]) %>%
  ungroup() %>%
  mutate(ComplexName_factor = forcats::fct_reorder(ComplexName, r_av_st, .desc = T),
         type = "GO Group") %>%
  arrange(ComplexName_factor) %>%
  dplyr::select(-r_av)

# combine go and complex results and track ordering
r_av_all = rbind(complex_r_av, go_r_av) %>%
  dplyr::mutate(ct_factor = factor(ct, levels = clusters))
complex_order = r_av_all %>%
  dplyr::group_by(ComplexName_factor) %>%
  dplyr::summarise(ComplexName = unique(ComplexName))

# summary tick mark info for protein complex test st
complex_st_tick_protein = gene_res %>%
  merge(complex_uni) %>%
  filter(ComplexName %in% complex_labels & ct == "St") %>%
  dplyr::group_by(ComplexName) %>%
  dplyr::mutate(complex_mean = mean(prot_av, na.rm = T)) %>%
  ungroup() %>%
  transmute(UNIPROT, ct, prot_av, ComplexName, complex_mean, type = "Protein Complex")

# summary tick mark info for protein go test st
go_st_tick_protein = gene_ref %>%
  filter(GO %in% GO_heatmap) %>%
  merge(transmute(test_res, ct, GO, ComplexName = TERM, complex_mean = prot_av)) %>%
  merge(transmute(gene_res, ct, UNIPROT, ct, prot_av, type = "GO Group")) %>%
  filter(ct == clusters[6]) %>%
  dplyr::select(-GO)

# marginal table displaying posterior probability associated with each group
tbl_df = r_av_all %>%
  ungroup() %>%
  filter(ct == "St") %>%
  dplyr::group_by(ComplexName_factor, type) %>%
  dplyr::summarise(Posterior_Probability = unique(fdr)) %>%
  dplyr::mutate(Posterior_Probability = ifelse(round(Posterior_Probability, 5) == 0, 0.00001, Posterior_Probability)) %>%
  distinct() %>%
  dplyr::select(Posterior_Probability, ComplexName_factor, type) %>%
  arrange(ComplexName_factor, .by_group = T)

# separate go and complex fdrs
complex_fdr_st_protein = tbl_df %>% filter(type == "Protein Complex")
go_fdr_st_protein = tbl_df %>% filter(type == "GO Group")

# tick mark summary info for protein go spc
complex_spc_tick_protein = gene_res %>%
  merge(complex_uni) %>%
  filter(ComplexName %in% complex_labels & ct == "SPC") %>%
  dplyr::group_by(ComplexName) %>%
  dplyr::mutate(complex_mean = mean(prot_av, na.rm = T)) %>%
  ungroup() %>%
  transmute(UNIPROT, ct, prot_av, ComplexName, complex_mean, type = "Protein Complex")

# tick mark summary info for protein complex spc
go_spc_tick_protein = gene_ref %>%
  filter(GO %in% GO_heatmap) %>%
  merge(transmute(test_res, ct, GO, ComplexName = TERM, complex_mean = prot_av)) %>%
  merge(transmute(gene_res, ct, UNIPROT, ct, prot_av, type = "GO Group")) %>%
  filter(ct == clusters[5]) %>%
  dplyr::select(-GO)

# marginal table displaying posterior probability associated with each group
tbl_df = r_av_all %>%
  ungroup() %>%
  filter(ct == "SPC") %>%
  dplyr::group_by(ComplexName_factor, type) %>%
  dplyr::summarise(Posterior_Probability = unique(fdr)) %>%
  dplyr::mutate(Posterior_Probability = ifelse(round(Posterior_Probability, 5) == 0, 0.00001, Posterior_Probability)) %>%
  distinct() %>%
  dplyr::select(Posterior_Probability, ComplexName_factor, type) %>%
  arrange(ComplexName_factor, .by_group = T)

# separate fdrs for complexes and go groups
complex_fdr_spc_protein = tbl_df %>% filter(type == "Protein Complex")
go_fdr_spc_protein = tbl_df %>% filter(type == "GO Group")


####################################################################################################################################### 
#                                                               Plot St GO heatmap    
####################################################################################################################################### 

#  go heatmap plot components (St)
r_av_all = go_heat_mrna %>% transmute(ComplexName, av = mu_av, ct, fdr, r_av_st, ComplexName_factor, type = "mRNA") %>%
  rbind(transmute(go_heat_protein, ComplexName, av = prot_av, ct, fdr, r_av_st, ComplexName_factor, type = "Protein")) %>%
  dplyr::mutate(ct_factor = factor(ct, levels = clusters))

# Plot go spermatid heatmap
# generate heatmap
g = ggplot(data = r_av_all, mapping = aes(x = ct_factor, y = ComplexName_factor, fill = av)) +
  geom_tile() +
  facet_wrap(vars(type), ncol = 1, scales = "free_y") +
  ylab("") +
  xlab("Cell Type") +
  ggtitle("GO Group Average (Posterior Mean)") +
  scale_y_discrete(labels = function(y) str_wrap(y, width = 25)) +
  scale_fill_gradientn(colors = c("darkblue", "blue", "white", "red", "darkred"),
                       limits = c(-1, 1),
                       oob = scales::squish,
                       na.value = NA, name = "rPTR", labels = c(-1, "", 0, "", 1)) +
  theme(panel.background = element_rect(fill = 'white', color = "white"),
        panel.grid.major = element_line(color = 'white'),
        panel.grid.minor = element_line(color = 'white'),
        panel.spacing = unit(2, "lines"),
        plot.title = element_text(size = text_size + 50),
        axis.text.y = element_text(size = text_size),
        axis.text.x = element_text(size = text_size + 50),
        axis.title = element_text(size = text_size + 50),
        strip.text = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"),
        legend.position = "none",
        legend.key.size = unit(4, "cm"))
g
# ggsave(plot = g, filename = "figs_results/heatmap_st_go.pdf", width = 45, height = 70, units = "in", limitsize = FALSE)

####################################################################################################################################### 
#                                                               Plot St Complex heatmap    
####################################################################################################################################### 

# combine mrna and protein complex info for st
r_av_all = complex_heat_mrna %>% transmute(ComplexName, av = mu_av, ct, fdr, r_av_st, ComplexName_factor, type = "mRNA") %>%
  rbind(transmute(complex_heat_protein, ComplexName, av = prot_av, ct, fdr, r_av_st, ComplexName_factor, type = "Protein")) %>%
  dplyr::mutate(ct_factor = factor(ct, levels = clusters))

# generate st complex heatmap
g = ggplot(data = r_av_all, mapping = aes(x = ct_factor, y = ComplexName_factor, fill = av)) +
  geom_tile() +
  facet_wrap(vars(type), ncol = 1, scales = "free_y") +
  ylab("") +
  xlab("Cell Type") +
  ggtitle("Complex Average (Posterior Mean)") +
  scale_y_discrete(labels = function(y) str_wrap(y, width = 25)) +
  scale_fill_gradientn(colors = c("darkblue", "blue", "white", "red", "darkred"),
                       limits = c(-1, 1),
                       oob = scales::squish,
                       na.value = NA, name = "rPTR", labels = c(-1, "", 0, "", 1)) +
  theme(panel.background = element_rect(fill = 'white', color = "white"),
        panel.grid.major = element_line(color = 'white'),
        panel.grid.minor = element_line(color = 'white'),
        panel.spacing = unit(2, "lines"),
        plot.title = element_text(size = text_size + 50),
        axis.text.y = element_text(size = text_size),
        axis.text.x = element_text(size = text_size + 50),
        axis.title = element_text(size = text_size + 50),
        strip.text = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"),
        legend.position = "none",
        legend.key.size = unit(4, "cm"))
g
# ggsave(plot = g, filename = "figs_results/heatmap_st_complex.pdf", width = 45, height = 70, units = "in", limitsize = FALSE)

####################################################################################################################################### 
#                                                               Plot St GO tick marks    
####################################################################################################################################### 

# combine mrna and protein tick mark go info st
uni_r_all = go_st_tick_mrna %>% transmute(UNIPROT, ct, ComplexName, complex_mean, av = mu_av, type = "mRNA") %>%
  rbind(transmute(go_st_tick_protein, UNIPROT, ct, ComplexName, complex_mean, av = prot_av, type = "Protein")) %>%
  merge(complex_order)

# use min and max across genes associated with all groups to auto-set axis limits
min_check = min(uni_r_all$av[is.finite(uni_r_all$av)], na.rm = T)
max_check = max(uni_r_all$av[is.finite(uni_r_all$av)], na.rm = T)
lim_val = max(abs(min_check), abs(max_check))

# draw tick mark plot
bb =  ggplot() +
  geom_point(data = uni_r_all,
             mapping = aes(x = av, y = ComplexName_factor, color = av, group = UNIPROT),
             size = 60, stroke = 60, shape="|", alpha = 1) +
  geom_point(data = uni_r_all, mapping = aes(x = complex_mean, y = ComplexName_factor, fill = complex_mean),
             size = 30, stroke = 10, shape = 24, color = "black") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed", size = 2) +
  facet_wrap(vars(type), ncol = 1, scales = "free_y") +
  xlab("Gene Product Posterior Mean") +
  ylab("") +
  xlim(-1*(lim_val), (lim_val)) +
  scale_color_gradientn(colors = c("darkblue", "blue", "white", "red", "darkred"),
                        limits = c(-1, 1),
                        oob = scales::squish,
                        na.value = NA, name = "", guide = "none") +
  scale_fill_gradientn(colors = c("darkblue", "blue", "white", "red", "darkred"),
                       limits = c(-1, 1),
                       oob = scales::squish,
                       na.value = NA, name = "", guide = "none") +
  scale_y_discrete(labels = function(y) str_wrap(y, width = 25)) +
  theme(panel.background = element_rect(fill = 'white', color = "white"),
        panel.spacing = unit(2, "lines"),
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = text_size + 50),
        axis.title = element_text(size = text_size + 50),
        strip.text = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = 'slategray1'),
        panel.grid.minor = element_line(color = 'white'))

# ggsave(plot = bb, filename = "figs_results/ticks_st_go.pdf", width = 30, height = 65, units = "in", limitsize = FALSE)
bb

####################################################################################################################################### 
#                                                             Plot St complex tick marks    
####################################################################################################################################### 

# combine mrna and protein tick mark complex info st
uni_r_all = complex_st_tick_mrna %>% transmute(UNIPROT, ct, ComplexName, complex_mean, av = mu_av, type = "mRNA") %>%
  rbind(transmute(complex_st_tick_protein, UNIPROT, ct, ComplexName, complex_mean, av = prot_av, type = "Protein")) %>%
  merge(complex_order)

# use min and max across genes associated with all groups to auto-set axis limits
min_check = min(uni_r_all$av[is.finite(uni_r_all$av)], na.rm = T)
max_check = max(uni_r_all$av[is.finite(uni_r_all$av)], na.rm = T)
lim_val = max(abs(min_check), abs(max_check))

# draw tick mark plot
bb =  ggplot() +
  geom_point(data = uni_r_all,
             mapping = aes(x = av, y = ComplexName_factor, color = av, group = UNIPROT),
             size = 60, stroke = 60, shape="|", alpha = 1) +
  geom_point(data = uni_r_all, mapping = aes(x = complex_mean, y = ComplexName_factor, fill = complex_mean),
             size = 30, stroke = 10, shape = 24, color = "black") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed", size = 2) +
  facet_wrap(vars(type), ncol = 1, scales = "free_y") +
  xlab("Gene Product Posterior Mean") +
  ylab("") +
  xlim(-1*(lim_val), (lim_val)) +
  scale_color_gradientn(colors = c("darkblue", "blue", "white", "red", "darkred"),
                        limits = c(-1, 1),
                        oob = scales::squish,
                        na.value = NA, name = "", guide = "none") +
  scale_fill_gradientn(colors = c("darkblue", "blue", "white", "red", "darkred"),
                       limits = c(-1, 1),
                       oob = scales::squish,
                       na.value = NA, name = "", guide = "none") +
  scale_y_discrete(labels = function(y) str_wrap(y, width = 25)) +
  theme(panel.background = element_rect(fill = 'white', color = "white"),
        panel.spacing = unit(2, "lines"),
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = text_size + 50),
        axis.title = element_text(size = text_size + 50),
        strip.text = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = 'slategray1'),
        panel.grid.minor = element_line(color = 'white'))

bb
# ggsave(plot = bb, filename = "figs_results/ticks_st_complex.pdf", width = 30, height = 65, units = "in", limitsize = FALSE)

####################################################################################################################################### 
#                                                               Plot spc GO heatmap    
####################################################################################################################################### 

# combine spc go tick mark info
uni_r_all = go_spc_tick_mrna %>% transmute(UNIPROT, ct, ComplexName, complex_mean, av = mu_av, type = "mRNA") %>%
  rbind(transmute(go_spc_tick_protein, UNIPROT, ct, ComplexName, complex_mean, av = prot_av, type = "Protein")) %>%
  merge(complex_order)

# use min and max across genes associated with all groups to auto-set axis limits
min_check = min(uni_r_all$av[is.finite(uni_r_all$av)], na.rm = T)
max_check = max(uni_r_all$av[is.finite(uni_r_all$av)], na.rm = T)
lim_val = max(abs(min_check), abs(max_check))

# draw tick mark plot
bb =  ggplot() +
  geom_point(data = uni_r_all,
             mapping = aes(x = av, y = ComplexName_factor, color = av, group = UNIPROT),
             size = 60, stroke = 60, shape="|", alpha = 1) +
  geom_point(data = uni_r_all, mapping = aes(x = complex_mean, y = ComplexName_factor, fill = complex_mean),
             size = 30, stroke = 10, shape = 24, color = "black") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed", size = 2) +
  facet_wrap(vars(type), ncol = 1, scales = "free_y") +
  xlab("Gene Product Posterior Mean") +
  ylab("") +
  xlim(-1*(lim_val), (lim_val)) +
  scale_color_gradientn(colors = c("darkblue", "blue", "white", "red", "darkred"),
                        limits = c(-1, 1),
                        oob = scales::squish,
                        na.value = NA, name = "", guide = "none") +
  scale_fill_gradientn(colors = c("darkblue", "blue", "white", "red", "darkred"),
                       limits = c(-1, 1),
                       oob = scales::squish,
                       na.value = NA, name = "", guide = "none") +
  scale_y_discrete(labels = function(y) str_wrap(y, width = 25)) +
  theme(panel.background = element_rect(fill = 'white', color = "white"),
        panel.spacing = unit(2, "lines"),
        legend.position = "bottom",
        # axis.text.y = element_text(size = text_size),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = text_size + 50),
        axis.title = element_text(size = text_size + 50),
        strip.text = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = 'slategray1'),
        panel.grid.minor = element_line(color = 'white'))
bb
# ggsave(plot = bb, filename = "figs_results/ticks_spc_go.pdf", width = 30, height = 65, units = "in", limitsize = FALSE)

####################################################################################################################################### 
#                                                               Plot spc complex heatmap    
####################################################################################################################################### 

# combine spc mrna and protein complex tick mark info
uni_r_all = complex_spc_tick_mrna %>% transmute(UNIPROT, ct, ComplexName, complex_mean, av = mu_av, type = "mRNA") %>%
  rbind(transmute(complex_spc_tick_protein, UNIPROT, ct, ComplexName, complex_mean, av = prot_av, type = "Protein")) %>%
  merge(complex_order)

# use min and max across genes associated with all groups to auto-set axis limits
min_check = min(uni_r_all$av[is.finite(uni_r_all$av)], na.rm = T)
max_check = max(uni_r_all$av[is.finite(uni_r_all$av)], na.rm = T)
lim_val = max(abs(min_check), abs(max_check))

# draw tick mark plot
bb =  ggplot() +
  geom_point(data = uni_r_all,
             mapping = aes(x = av, y = ComplexName_factor, color = av, group = UNIPROT),
             size = 60, stroke = 60, shape="|", alpha = 1) +
  geom_point(data = uni_r_all, mapping = aes(x = complex_mean, y = ComplexName_factor, fill = complex_mean),
             size = 30, stroke = 10, shape = 24, color = "black") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed", size = 2) +
  facet_wrap(vars(type), ncol = 1, scales = "free_y") +
  xlab("Gene Product Posterior Mean") +
  ylab("") +
  xlim(-1*(lim_val), (lim_val)) +
  scale_color_gradientn(colors = c("darkblue", "blue", "white", "red", "darkred"),
                        limits = c(-1, 1),
                        oob = scales::squish,
                        na.value = NA, name = "", guide = "none") +
  scale_fill_gradientn(colors = c("darkblue", "blue", "white", "red", "darkred"),
                       limits = c(-1, 1),
                       oob = scales::squish,
                       na.value = NA, name = "", guide = "none") +
  scale_y_discrete(labels = function(y) str_wrap(y, width = 25)) +
  theme(panel.background = element_rect(fill = 'white', color = "white"),
        panel.spacing = unit(2, "lines"),
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = text_size + 50),
        axis.title = element_text(size = text_size + 50),
        strip.text = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = 'slategray1'),
        panel.grid.minor = element_line(color = 'white'))

bb
# ggsave(plot = bb, filename = "figs_results/ticks_spc_complex.pdf", width = 30, height = 65, units = "in", limitsize = FALSE)

####################################################################################################################################### 
#                                                               Plot St fdr tables    
####################################################################################################################################### 

# st go test fdr info
tbl_df = go_fdr_st_mrna %>% dplyr::mutate(type = "mRNA") %>%
  rbind(mutate(go_fdr_st_protein, type = "Protein"))
tbl = ggplot(tbl_df, aes(y = fct_inorder(ComplexName_factor), x = 1)) +
  geom_text(mapping = aes(x = 1, y = fct_inorder(ComplexName_factor),
                          label = as.character(vaxedemic::scientific_10x(Posterior_Probability, digits = 0))),
            size = 50, parse = T) +
  facet_wrap(vars(type), ncol = 1, scales = "free_y") +
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = 'white', color = "white"),
        panel.spacing = unit(2, "lines"),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_line(color = 'white'),
        panel.grid.minor = element_line(color = 'white'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"))
# ggsave(plot = tbl, filename = "figs_results/fdr_st_go.pdf", width = 10, height = 60, units = "in", limitsize = FALSE)
tbl

# st complex fdr go info
tbl_df = complex_fdr_st_mrna %>% dplyr::mutate(type = "mRNA") %>%
  rbind(mutate(complex_fdr_st_protein, type = "Protein"))
tbl = ggplot(tbl_df, aes(y = fct_inorder(ComplexName_factor), x = 1)) +
  geom_text(mapping = aes(x = 1, y = fct_inorder(ComplexName_factor),
                          label = as.character(vaxedemic::scientific_10x(Posterior_Probability, digits = 0))),
            size = 50, parse = T) +
  facet_wrap(vars(type), ncol = 1, scales = "free_y") +
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = 'white', color = "white"),
        panel.spacing = unit(2, "lines"),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_line(color = 'white'),
        panel.grid.minor = element_line(color = 'white'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"))
# ggsave(plot = tbl, filename = "figs_results/fdr_st_complex.pdf", width = 10, height = 60, units = "in", limitsize = FALSE)
tbl

####################################################################################################################################### 
#                                                               Plot SPC fdr table    
####################################################################################################################################### 

# go test fdr info spc 
tbl_df = go_fdr_spc_mrna %>% dplyr::mutate(type = "mRNA") %>%
  rbind(mutate(go_fdr_spc_protein, type = "Protein"))
tbl = ggplot(tbl_df, aes(y = fct_inorder(ComplexName_factor), x = 1)) +
  geom_text(mapping = aes(x = 1, y = fct_inorder(ComplexName_factor),
                          label = as.character(vaxedemic::scientific_10x(Posterior_Probability, digits = 0))),
            size = 50, parse = T) +
  facet_wrap(vars(type), ncol = 1, scales = "free_y") +
  xlab("") +
  ylab("") +
  theme(panel.background = element_rect(fill = 'white', color = "white"),
        panel.spacing = unit(2, "lines"),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_line(color = 'white'),
        panel.grid.minor = element_line(color = 'white'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"))
# ggsave(plot = tbl, filename = "figs_results/fdr_spc_go.pdf", width = 10, height = 60, units = "in", limitsize = FALSE)
tbl

# complex test fdr info spc
tbl_df = complex_fdr_spc_mrna %>% dplyr::mutate(type = "mRNA") %>%
  rbind(mutate(complex_fdr_spc_protein, type = "Protein"))
tbl = ggplot(tbl_df, aes(y = fct_inorder(ComplexName_factor), x = 1)) +
  geom_text(mapping = aes(x = 1, y = fct_inorder(ComplexName_factor),
                          label = as.character(vaxedemic::scientific_10x(Posterior_Probability, digits = 0))),
            size = 50, parse = T) +
  facet_wrap(vars(type), ncol = 1, scales = "free_y") +
  xlab("") +
  ylab("") + 
  theme(panel.background = element_rect(fill = 'white', color = "white"),
        panel.spacing = unit(2, "lines"),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_line(color = 'white'),
        panel.grid.minor = element_line(color = 'white'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"))
# ggsave(plot = tbl, filename = "figs_results/fdr_spc_complex.pdf", width = 10, height = 60, units = "in", limitsize = FALSE)
tbl