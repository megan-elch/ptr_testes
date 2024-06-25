data{
  int<lower=0> n_g; //number of gene products
  int<lower=0> total_obs_mrna; //total number of mrna obs
  int<lower=0> total_obs_protein; //total number of protein obs
  
  array[total_obs_mrna] int<lower=0> mrna_sum; //mrna sum of transcript counts across cells in each cluster
  vector[total_obs_protein] protein_av; //protein average across cells in each clusters

  vector<lower=0>[total_obs_mrna] n_mrna; //number of mrna data sets
  vector<lower=0>[total_obs_protein] n_protein; //number of protein data sets

  array[total_obs_mrna] int<lower=0> mrna_gene_rec; //map each mrna observation to its gene
  array[total_obs_protein] int<lower=0> protein_gene_rec; //map each protein observation to gene
  
  array[total_obs_mrna] int<lower=0> mrna_celltype_rec; //map each mrna obs to cell type
  array[total_obs_protein] int<lower=0> protein_celltype_rec; //map each potein obs to cell type

  array[total_obs_mrna] int<lower=0> mrna_pop_rec; //map each mrna obs to data set
  array[total_obs_protein] int<lower=0> protein_pop_rec; //map each protein obs to data set
  array[total_obs_protein] int<lower=0> protein_measure_rec; //map each protein obs to scope or plexdia
}

transformed data{
  int<lower=0> n_ct;
  int<lower=0> n_pop_mrna;
  int<lower=0> n_pop_protein;
  int<lower=0> n_measure_protein;

  n_ct = max(protein_celltype_rec); //total number of cell types
  n_pop_mrna = max(mrna_pop_rec); //total number of mrna data sets
  n_pop_protein = max(protein_pop_rec); //total number of protein data sets
  n_measure_protein = max(protein_measure_rec); //total number of protein measurement types 
}

parameters{
  matrix[n_g, n_ct] r; //parameter for protein to rna ratio rptr
  matrix[n_g, n_ct] mu; //parameter for consensus transcript abundance
  vector<lower=0>[n_g] sigma_r; //gene-level variance of rptr
  vector<lower=0>[n_g] sigma_mu; //gene-level variance of mrna
  matrix<lower=0>[n_measure_protein, n_g] scale_param; //scale linking protein to rna
  
  matrix[n_pop_mrna, n_g] gamma_m; //gene level mrna technical effect
  matrix[n_pop_mrna, n_ct] a_mrna; //cluster level mrna technical effect
  matrix[n_pop_protein, n_g] kappa; //gene level protein technical effect
  matrix[n_pop_protein, n_ct] b_protein; //cluster level protein technical effect
  
  matrix<lower=0>[n_pop_protein, n_g] sigma_p_sq; //data set, gene product error variance for protein
  vector<lower=0>[n_pop_protein] sigma_d; //data set level hierarchical mean for protein error variance
  vector<lower=0>[n_pop_protein] sd_var; //data set level hierarchical sd for protein error variance
  vector<lower=0>[n_g] phi_m; //gene level mrna overdispersion
}

 transformed parameters{
  vector<lower=0>[total_obs_mrna] mrna_mean;
  vector<lower=0>[total_obs_mrna] phi_transform;
  vector[total_obs_protein] protein_mean;
  vector<lower=0>[total_obs_protein] protein_sd;
  
  for(i in 1:total_obs_mrna){
    phi_transform[i] = pow(phi_m[mrna_gene_rec[i]], -2.0); //transform mrna overdisperison 
    mrna_mean[i] = exp2(mu[mrna_gene_rec[i], mrna_celltype_rec[i]] + gamma_m[mrna_pop_rec[i], mrna_gene_rec[i]] + a_mrna[mrna_pop_rec[i], mrna_celltype_rec[i]]); //define mrna by signal and technical effects on log scale
  }
  
  for(i in 1:total_obs_protein){
    protein_mean[i] = inv(scale_param[protein_measure_rec[i], protein_gene_rec[i]])*(mu[protein_gene_rec[i], protein_celltype_rec[i]] + r[protein_gene_rec[i], protein_celltype_rec[i]]) + kappa[protein_pop_rec[i], protein_gene_rec[i]] + b_protein[protein_pop_rec[i], protein_celltype_rec[i]]; //define protein mean
    protein_sd[i] = sqrt(sigma_p_sq[protein_pop_rec[i], protein_gene_rec[i]]); //protein standard deviation
  } 
}

model{
  //set priors
  to_vector(scale_param) ~ gamma(2, 1);
  sigma_mu ~ normal(0, 1);
  sigma_r ~ normal(0, 1);
  
  phi_m ~ normal(0, 1);
  sigma_d ~ normal(0, 1);
  sd_var ~ normal(0, 1);
  to_vector(sigma_p_sq) ~ normal(to_vector(rep_matrix(sigma_d, n_g)), to_vector(rep_matrix(sd_var, n_g)));
 
  to_vector(gamma_m) ~ normal(0, 10);
  to_vector(a_mrna) ~ normal(0, 10);
  to_vector(kappa) ~ normal(0, 1);
  to_vector(b_protein) ~ normal(0, 1);
  
  to_vector(mu) ~ normal(rep_vector(0, n_g*n_ct), to_vector(rep_matrix(sqrt(sigma_mu), n_ct)));
  to_vector(r) ~ normal(rep_vector(0, n_g*n_ct), to_vector(rep_matrix(sqrt(sigma_r), n_ct)));
    
  mrna_sum ~ neg_binomial_2(mrna_mean, phi_transform);
  protein_av ~ normal(protein_mean, protein_sd);
}