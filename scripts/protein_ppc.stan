data{
  int<lower=0> n_g;
  int<lower=0> total_obs_mrna;
  int<lower=0> total_obs_protein; 
  
  array[total_obs_mrna] int<lower=0> mrna_sum;
  vector[total_obs_protein] protein_av;

  vector<lower=0>[total_obs_mrna] n_mrna;
  vector<lower=0>[total_obs_protein] n_protein;

  array[total_obs_mrna] int<lower=0> mrna_gene_rec;
  array[total_obs_protein] int<lower=0> protein_gene_rec;
  
  array[total_obs_mrna] int<lower=0> mrna_celltype_rec;
  array[total_obs_protein] int<lower=0> protein_celltype_rec;

  array[total_obs_mrna] int<lower=0> mrna_pop_rec;
  array[total_obs_protein] int<lower=0> protein_pop_rec;
  array[total_obs_protein] int<lower=0> protein_measure_rec;
}

transformed data{
  int<lower=0> n_ct;
  int<lower=0> n_pop_mrna;
  int<lower=0> n_pop_protein;
  int<lower=0> n_measure_protein;
  int<lower=0> total_obs;

  n_ct = max(protein_celltype_rec);
  n_pop_mrna = max(mrna_pop_rec);
  n_pop_protein = max(protein_pop_rec);
  n_measure_protein = max(protein_measure_rec);
  total_obs = total_obs_protein + total_obs_mrna;
}

parameters{
  matrix[n_g, n_ct] r;
  matrix[n_g, n_ct] mu;
  vector<lower=0>[n_g] sigma_r;
  vector<lower=0>[n_g] sigma_mu;
  matrix<lower=0>[n_measure_protein, n_g] scale_param;
  
  matrix[n_pop_mrna, n_g] gamma_m;
  matrix[n_pop_mrna, n_ct] a_mrna;
  matrix[n_pop_protein, n_g] kappa;
  matrix[n_pop_protein, n_ct] b_protein;
  
  vector<lower=0>[n_pop_protein] sigma_d;
  matrix<lower=0>[n_pop_protein, n_g] sigma_p;
  vector<lower=0>[n_pop_protein] sd_var;
  vector<lower=0>[n_g] phi_m;
}

 transformed parameters{
  vector<lower=0>[total_obs_mrna] mrna_mean;
  vector<lower=0>[total_obs_mrna] phi_transform;
  vector[total_obs_protein] protein_mean;
  vector<lower=0>[total_obs_protein] protein_sd;
  
  for(i in 1:total_obs_mrna){
    phi_transform[i] = pow(phi_m[mrna_gene_rec[i]], -2.0);
    mrna_mean[i] = exp2(mu[mrna_gene_rec[i], mrna_celltype_rec[i]] + gamma_m[mrna_pop_rec[i], mrna_gene_rec[i]] + a_mrna[mrna_pop_rec[i], mrna_celltype_rec[i]]);
  }
  
  for(i in 1:total_obs_protein){
    protein_mean[i] = inv(scale_param[protein_measure_rec[i], protein_gene_rec[i]])*(mu[protein_gene_rec[i], protein_celltype_rec[i]] + r[protein_gene_rec[i], protein_celltype_rec[i]]) + kappa[protein_pop_rec[i], protein_gene_rec[i]] + b_protein[protein_pop_rec[i], protein_celltype_rec[i]];
    protein_sd[i] = sqrt(sigma_p[protein_pop_rec[i], protein_gene_rec[i]]);
  } 
}

generated quantities{
  array[total_obs_protein] real protein_bar_rep;
  vector[total_obs_protein] protein_log_lik;
  vector[total_obs_mrna] mrna_log_lik;
  vector[total_obs] log_lik;

for(n in 1:total_obs_protein){
  protein_log_lik[n] = normal_lpdf(protein_av[n] | protein_mean[n], protein_sd[n]);
}
for(n in 1:total_obs_mrna){
  mrna_log_lik[n] = neg_binomial_2_lpmf(mrna_sum[n] | mrna_mean[n], phi_transform[n]);
}

  log_lik = append_row(protein_log_lik, mrna_log_lik);
  protein_bar_rep = normal_rng(protein_mean, protein_sd);
} 