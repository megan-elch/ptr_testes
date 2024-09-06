data{
  int<lower=0> n_g; // number of gene product pairs
  int<lower=0> total_obs_mrna; // total number of mrna observations
  int<lower=0> total_obs_protein; // total number of protein observations
  
  array[total_obs_mrna] int<lower=0> mrna_sum; // sum of transcript counts
  vector[total_obs_protein] protein_av; // average log2 peptide intensity

  vector<lower=0>[total_obs_mrna] n_mrna; // number of mrna observations
  vector<lower=0>[total_obs_protein] n_protein; // number of protein observations

  array[total_obs_mrna] int<lower=0> mrna_gene_rec; // map mrna observations to gene
  array[total_obs_protein] int<lower=0> protein_gene_rec; // map protein observations to gene
  
  array[total_obs_mrna] int<lower=0> mrna_celltype_rec; // map mrna obs to cell type
  array[total_obs_protein] int<lower=0> protein_celltype_rec; // map protein obs to cell type

  array[total_obs_mrna] int<lower=0> mrna_pop_rec; // map mrna obs to data set
  array[total_obs_protein] int<lower=0> protein_pop_rec; // map protein obs to data set
  array[total_obs_protein] int<lower=0> protein_measure_rec; // map protein obs to measurment type
}

transformed data{
  int<lower=0> n_ct; // number of cell types 
  int<lower=0> n_pop_mrna; // number of mrna data sets
  int<lower=0> n_pop_protein; // number of protein data sets
  int<lower=0> n_measure_protein; // number of protein measurement types

  n_ct = max(protein_celltype_rec); // number of cell types
  n_pop_mrna = max(mrna_pop_rec); // number of mrna data sets
  n_pop_protein = max(protein_pop_rec); // number of protein data sets
  n_measure_protein = max(protein_measure_rec); // number of protein measurement types
}

parameters{
  matrix[n_g, n_ct] r; // protein to mrna ratio
  matrix[n_g, n_ct] mu; // consensus mrna
  matrix<lower=0>[n_measure_protein, n_g] scale_param; // scale parameter linking mrna, protein
  
  matrix[n_pop_mrna, n_g] gamma_m; // gene level normalization term
  matrix[n_pop_mrna, n_ct] a_mrna; // cluster level normalization term 
  matrix[n_pop_protein, n_g] kappa; // gene level normalization term 
  matrix[n_pop_protein, n_ct] b_protein; // cluster level normalization term
  
  vector<lower=0>[n_pop_protein] sigma_d; // hyperprior mean for protein level sampling sd
  matrix<lower=0>[n_pop_protein, n_g] sigma_p; // protein level sampling var
  vector<lower=0>[n_pop_protein] sd_var; // hyperprior sd for protein level sampling sd
  vector<lower=0>[n_g] phi_m; // mrna overdispersion
}

 transformed parameters{
  vector<lower=0>[total_obs_mrna] mrna_mean; // average transcript abundance param
  vector<lower=0>[total_obs_mrna] phi_transform; // vector form of overdispersion param
  vector[total_obs_protein] protein_mean; // average protein intensity param 
  vector<lower=0>[total_obs_protein] protein_sd; // vector form of protein sampling sd
  
  for(i in 1:total_obs_mrna){
    // transform overdispersion param to prevent extreme overdispersion a priori mass
    phi_transform[i] = pow(phi_m[mrna_gene_rec[i]], -2.0);
    // build average transcript abundance param
    mrna_mean[i] = exp2(mu[mrna_gene_rec[i], mrna_celltype_rec[i]] + gamma_m[mrna_pop_rec[i], mrna_gene_rec[i]] + a_mrna[mrna_pop_rec[i], mrna_celltype_rec[i]]);
  }
  
  for(i in 1:total_obs_protein){
    // build average protein intensity param 
    protein_mean[i] = inv(scale_param[protein_measure_rec[i], protein_gene_rec[i]])*(mu[protein_gene_rec[i], protein_celltype_rec[i]] + r[protein_gene_rec[i], protein_celltype_rec[i]]) + kappa[protein_pop_rec[i], protein_gene_rec[i]] + b_protein[protein_pop_rec[i], protein_celltype_rec[i]];
    // protein sampling standard dev 
    protein_sd[i] = sqrt(sigma_p[protein_pop_rec[i], protein_gene_rec[i]]);
  } 
}

model{
  // model inverse scale param
  to_vector(scale_param) ~ gamma(2, 1); 

  phi_m ~ normal(0, 1);
  sigma_d ~ normal(0, 1);
  sd_var ~ normal(0, 1);
  // protein sampling var distributed across data set level mean
  to_vector(sigma_p) ~ normal(to_vector(rep_matrix(sigma_d, n_g)), to_vector(rep_matrix(sd_var, n_g))); 
 
  to_vector(gamma_m) ~ normal(0, 10);
  to_vector(a_mrna) ~ normal(0, 10);
  to_vector(kappa) ~ normal(0, 1);
  to_vector(b_protein) ~ normal(0, 1);
  
  to_vector(mu) ~ normal(0, 1); 
  to_vector(r) ~ normal(0, 1); 
    
  mrna_sum ~ neg_binomial_2(mrna_mean, phi_transform);
  protein_av ~ normal(protein_mean, protein_sd);
}