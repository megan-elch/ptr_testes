# ptr_testes
Code for BayesPG Model, developed with Slavov and Franks Labs.

BayesPG allows for the quantification of post transcriptional regulation by disentangling biological and technical variation to compute the relative protein-to-RNA ratio (rPTR). 

* eda.Rmd: Exploration of protein consistency and data set agreement metrics, z-tranformed cluster-level fold changes.
* results.Rmd: Examination of model results, including posterior means, gene, GO, and complex test results. 
* ppc.Rmd: Comparison of posterior predictive statistics to observed.
* scripts folder: contains supporting R scripts for running BayesPG and exploring results. Includes:
  * prep_testes.R: Data preparation and formatting prior to model application
  * testes.R: Running BayesPG
  * hs.stan: Stan code used for BayesPG
  * results.R: Gene, GO, and complex-level significance testing
  * model_functions.R: Utility functions for running model and gene-level testing.
  * go_functions.R: Utility functions related to GO testing 
  * dbi_info.R: Organizes AnnotationDbi information related to GO testing
  * mrna_ppc.stan: Stan code for generating posterior predictive draws for transcript data
  * protein_ppc.stan: Stan code for generating posterior predictive draws for protein data
  * ppc.R: Collect posterior predictive draws and summary statistics.
  * ppc_functions.R: Utility functions related to posterior predictive checking
  * plot_functions.R: Utility functions for plotting
  * protein_reliability.R: Protein consistency computations.
  * reliability_functions.R: Utility functions for reliability computations.
