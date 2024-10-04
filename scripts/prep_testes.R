library(tidyverse)
library(Matrix)
library(rlist)
library(readxl)
library(org.Hs.eg.db)  

# identify project
source("scripts/model_functions.R")
clusters = c("EC", "LC", "PTM", "SPC", "SPG", "St") # identify cell types order

# load previously filtered mrna data
load(paste0("data/scRNA_seq/mrna_prepped_align.RData"))

# load in peptide data sets 
peptides_pop1 = read.table(paste0("data/scProtein/Pop1/testis_nPoP.181122_peptideLevelMatrix_dartUpdate_medIntNorm_BCAtt3.noProteo.txt"), 
                           sep = '\t', header = TRUE)
peptides_pop2 = read.table(paste0("data/scProtein/Pop2/testis_nPoP.081221_peptideLevelMatrix_dartUpdate.medIntNorm.noProteo.txt"), 
                           sep = '\t', header = TRUE)
peptides_pop3 = read.table(paste0("data/scProtein/Pop3/testis_nPoP.130122_peptideLevelMatrix_dartUpdate_medIntNormNoBC.noProteo.txt"), 
                           sep = '\t', header = TRUE)
peptides_pop4 = read.table(paste0("data/scProtein/Pop4/testis_nPoP.030222_peptideLevelMatrix_dartUpdate_medIntNorm.txt"), 
                           sep = '\t', header = TRUE)
peptides_pop5 = read.table(paste0("data/scProtein/Pop5/testis_nPoP.130124_peptideLevelMatrix_BC.txt"), 
                           sep = '\t', header = TRUE)

# organize column names across data sets
colnames(peptides_pop1) = ifelse(colnames(peptides_pop1) == "prot", "prot", ifelse(colnames(peptides_pop1) == "pep", "pep", paste0("One_", colnames(peptides_pop1))))
colnames(peptides_pop2) = ifelse(colnames(peptides_pop2) == "prot", "prot", ifelse(colnames(peptides_pop2) == "pep", "pep", paste0("Two_", colnames(peptides_pop2))))
colnames(peptides_pop3) = ifelse(colnames(peptides_pop3) == "prot", "prot", ifelse(colnames(peptides_pop3) == "pep", "pep", paste0("Three_", colnames(peptides_pop3))))
colnames(peptides_pop4) = ifelse(colnames(peptides_pop4) == "prot", "prot", ifelse(colnames(peptides_pop4) == "pep", "pep", paste0("Four_", colnames(peptides_pop4))))
colnames(peptides_pop5) = ifelse(colnames(peptides_pop5) == "prot", "prot", ifelse(colnames(peptides_pop5) == "pep", "pep", paste0("Five_", colnames(peptides_pop5))))

# protein meta data
protein_meta = read.delim(paste0("data/scProtein/npop_cellIDs_pDIA.RNABC_kNNAll.txt"))
protein_info = data.frame(ct = protein_meta$cellType, id = protein_meta$id)

# filter for each data set
protein_pop1_info = protein_info %>% filter(id %in% colnames(peptides_pop1))
protein_pop2_info = protein_info %>% filter(id %in% colnames(peptides_pop2))
protein_pop3_info = protein_info %>% filter(id %in% colnames(peptides_pop3))
protein_pop4_info = protein_info %>% filter(id %in% colnames(peptides_pop4))
protein_pop5_info = protein_info %>% filter(id %in% colnames(peptides_pop5))

# print summary stats about data
n_transcripts = length(union(rownames(mrna_pop1_counts), rownames(mrna_pop2_counts)))
print(paste0("Number of Transcripts ", n_transcripts))

n_cells_m = length(mrna_pop1_info$id) + length(mrna_pop2_info$id)
print(paste0("Number of Cells (mRNA) ", n_cells_m))

n_proteins = union(peptides_pop1$prot, peptides_pop2$prot) %>% union(peptides_pop3$prot) %>% union(peptides_pop4$prot) %>% union(peptides_pop5$prot) %>% length()
print(paste0("Number of Proteins ", n_proteins))

n_peptides = union(peptides_pop1$pep, peptides_pop2$pep) %>% union(peptides_pop3$pep) %>% union(peptides_pop4$pep) %>% union(peptides_pop5$pep) %>% length()
print(paste0("Number of Peptides ", n_peptides))

n_cells_p = length(protein_info$id) 
print(paste0("Number of Cells (Proteins) ", n_cells_p))

# we input transcript count sums and log2 peptide intensity averages for each gene product, data set, cluster.
# since Stan does not accept NAs, we must flatten data into one long vector for each modality and remove NAs
# during prep steps we will also generate reference vectors that map our observations to corresponding data sets,
# cluster labels etc (see testes.R)

# compute summary statistics for protein data (peptide level averages for each cluster and data set)
protein_suff_pop1 = compute_suff_stats(peptides_pop1, protein_pop1_info, type = "protein", clusters = clusters, pop_label = "Pop1") %>%
    dplyr::mutate(measure = "scope")
print("protein1")
protein_suff_pop2 = compute_suff_stats(peptides_pop2, protein_pop2_info, type = "protein", clusters = clusters, pop_label = "Pop2") %>%
    dplyr::mutate(measure = "scope")
print("protein2")
protein_suff_pop3 = compute_suff_stats(peptides_pop3, protein_pop3_info, type = "protein", clusters = clusters, pop_label = "Pop3") %>%
    dplyr::mutate(measure = "scope")
print("protein3")
protein_suff_pop4 = compute_suff_stats(peptides_pop4, protein_pop4_info, type = "protein", clusters = clusters, pop_label = "Pop4") %>%
    dplyr::mutate(measure = "scope")
print("protein4")
protein_suff_pop5 = compute_suff_stats(peptides_pop5, protein_pop5_info, type = "protein", clusters = clusters, pop_label = "Pop5") %>%
    dplyr::mutate(measure = "plex")
print("protein5")
protein_suff = rbind(protein_suff_pop1, protein_suff_pop2) %>% rbind(protein_suff_pop3) %>% rbind(protein_suff_pop4) %>% rbind(protein_suff_pop5)
print(head(protein_suff))

save(protein_suff, file = paste0("hs/data/suff_stats_prot.RData"))

# compute sufficient statistics for mrna data (transcript count sums across cells for each gene, cell type, data set)
mrna_suff_pop1 = compute_suff_stats(mrna_pop1_counts, mrna_pop1_info, clusters = clusters, pop_label = "Pop1")
print("mrna1")
mrna_suff_pop2 = compute_suff_stats(mrna_pop2_counts, mrna_pop2_info, clusters = clusters, pop_label = "Pop2")
print("mrna2")
mrna_suff_pop3 = compute_suff_stats(mrna_pop3_counts, mrna_pop3_info, clusters = clusters, pop_label = "Pop3")
print("mrna3")
mrna_suff_pop4 = compute_suff_stats(mrna_pop4_counts, mrna_pop4_info, clusters = clusters, pop_label = "Pop4")
print("mrna4")
mrna_suff_pop5 = compute_suff_stats(mrna_pop5_counts, mrna_pop5_info, clusters = clusters, pop_label = "Pop5")
print("mrna5")
mrna_suff = rbind(mrna_suff_pop1, mrna_suff_pop2) %>%
    rbind(mrna_suff_pop3, mrna_suff_pop4, mrna_suff_pop5)
save(mrna_suff, file = paste0("data/suff_stats_mrna.RData"))
print(head(mrna_suff))

# link together labels for each modality
protein_genes_union = unique(protein_suff$UNIPROT)
mrna_genes_union = unique(mrna_suff$SYMBOL)

protein_genes = AnnotationDbi::select(org.Hs.eg.db, keys = protein_genes_union,
                                      keytype ='UNIPROT', columns = c("SYMBOL"))
protein_genes = protein_genes[!duplicated(protein_genes[,"UNIPROT"]),]
mrna_genes = AnnotationDbi::select(org.Hs.eg.db, keys= mrna_genes_union,
                                   keytype = 'SYMBOL', columns = c("UNIPROT"))

# combine with sufficient statistics
all_genes = merge(protein_genes, mrna_genes)
protein_suff = merge(protein_suff, all_genes)
mrna_suff = merge(mrna_suff, all_genes)

save(mrna_suff, protein_suff, file = paste0("data/suff_stats.RData"))