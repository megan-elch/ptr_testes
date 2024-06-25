library(tidyverse)
library(Matrix)
library(rlist)
library(readxl)
library(org.Hs.eg.db)  

# identify project
project_name = "testes"
source("scripts/model_functions.R")
source("scripts/reliability_functions.R")
clusters = c("EC", "LC", "PTM", "SPC", "SPG", "St") # identify cell types order
plot_path = paste0("figs_results/")
output_path = paste0("processed_data/")
print("loaded")

load(paste0("data/", project_name, "/scRNA_seq/mrna_prepped_align.RData"))

# load in peptide data sets 
peptides_pop1 = read.table(paste0("data/", project_name, "/scProtein/Pop1/testis_nPoP.181122_peptideLevelMatrix_dartUpdate_medIntNorm_BCAtt3.noProteo.txt"), 
                           sep = '\t', header = TRUE)
peptides_pop2 = read.table(paste0("data/", project_name, "/scProtein/Pop2/testis_nPoP.081221_peptideLevelMatrix_dartUpdate.medIntNorm.noProteo.txt"), 
                           sep = '\t', header = TRUE)
peptides_pop3 = read.table(paste0("data/", project_name, "/scProtein/Pop3/testis_nPoP.130122_peptideLevelMatrix_dartUpdate_medIntNormNoBC.noProteo.txt"), 
                           sep = '\t', header = TRUE)
peptides_pop4 = read.table(paste0("data/", project_name, "/scProtein/Pop4/testis_nPoP.030222_peptideLevelMatrix_dartUpdate_medIntNorm.txt"), 
                           sep = '\t', header = TRUE)
peptides_pop5 = read.table(paste0("data/", project_name, "/scProtein/Pop5/testis_nPoP.130124_peptideLevelMatrix_BC.txt"), 
                           sep = '\t', header = TRUE)

colnames(peptides_pop1) = ifelse(colnames(peptides_pop1) == "prot", "prot", ifelse(colnames(peptides_pop1) == "pep", "pep", paste0("One_", colnames(peptides_pop1))))
colnames(peptides_pop2) = ifelse(colnames(peptides_pop2) == "prot", "prot", ifelse(colnames(peptides_pop2) == "pep", "pep", paste0("Two_", colnames(peptides_pop2))))
colnames(peptides_pop3) = ifelse(colnames(peptides_pop3) == "prot", "prot", ifelse(colnames(peptides_pop3) == "pep", "pep", paste0("Three_", colnames(peptides_pop3))))
colnames(peptides_pop4) = ifelse(colnames(peptides_pop4) == "prot", "prot", ifelse(colnames(peptides_pop4) == "pep", "pep", paste0("Four_", colnames(peptides_pop4))))
colnames(peptides_pop5) = ifelse(colnames(peptides_pop5) == "prot", "prot", ifelse(colnames(peptides_pop5) == "pep", "pep", paste0("Five_", colnames(peptides_pop5))))

# protein meta data
protein_meta = read.delim(paste0("data/", project_name, "/scProtein/npop_cellIDs_pDIA.RNABC_kNNAll.txt"))
protein_info = data.frame(ct = protein_meta$cellType, id = protein_meta$id)
protein_pop1_info = protein_info %>% filter(id %in% colnames(peptides_pop1))
protein_pop2_info = protein_info %>% filter(id %in% colnames(peptides_pop2))
protein_pop3_info = protein_info %>% filter(id %in% colnames(peptides_pop3))
protein_pop4_info = protein_info %>% filter(id %in% colnames(peptides_pop4))
protein_pop5_info = protein_info %>% filter(id %in% colnames(peptides_pop5))

# compute protein consistency
protein_pop1_reliability = compute_protein_consistency(peptides_pop1, protein_pop1_info) %>%
   mutate(pop_protein = "Pop1")
protein_pop2_reliability = compute_protein_consistency(peptides_pop2, protein_pop2_info) %>%
   mutate(pop_protein = "Pop2")
protein_pop3_reliability = compute_protein_consistency(peptides_pop3, protein_pop3_info) %>%
   mutate(pop_protein = "Pop3")
protein_pop4_reliability = compute_protein_consistency(peptides_pop4, protein_pop4_info) %>%
   mutate(pop_protein = "Pop4") 
protein_pop5_reliability = compute_protein_consistency(peptides_pop5, protein_pop5_info) %>%
  mutate(pop_protein = "Pop5") 
protein_reliability = rbind(protein_pop1_reliability, protein_pop2_reliability) %>%
   rbind(protein_pop3_reliability) %>%
   rbind(protein_pop4_reliability) %>% 
   rbind(protein_pop5_reliability)
save(protein_reliability, file = paste0(output_path, "protein_consistency.RData"))

