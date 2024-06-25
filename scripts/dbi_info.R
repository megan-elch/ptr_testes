library(tidyverse)
library(AnnotationDbi)
library(GO.db)
library(org.Hs.eg.db)

load("model_output/gene_res.RData")

# recover go groups associated with proteins fit in model
GO_info = AnnotationDbi::select(org.Hs.eg.db, keys = unique(gene_res$UNIPROT), keytype = "UNIPROT", columns = "GO") %>%
    dplyr::select(-c(EVIDENCE, ONTOLOGY)) %>%
    distinct(.keep_all = TRUE) 
  
# recover all proteins associated with present go groups
GO_info2 = AnnotationDbi::select(org.Hs.eg.db, keys = unique(GO_info$GO), keytype = "GO", columns = "UNIPROT") %>%
    dplyr::group_by(GO) %>%
    dplyr::summarise(total_terms = length(unique(UNIPROT)),
                     nd = length(unique(intersect(gene_res$UNIPROT, UNIPROT)))) %>%
    ungroup() 

# recover go "terms" associated with groups (more detailed labels)
GO_tab = AnnotationDbi::select(GO.db, keys = unique(GO_info$GO), keytype = "GOID", columns = "TERM") %>%
  ungroup() %>%
  dplyr::select(GOID, TERM) %>%
  distinct(.keep_all = TRUE) 

save(GO_info, GO_info2, GO_tab, file = "processed_data/GO_prep.RData")
