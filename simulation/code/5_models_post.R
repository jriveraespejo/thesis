# # preliminar ####
# 
# # erase all objects
# rm(list=ls())
# 
# # load libraries
# librerias <- c('stringr','dplyr','ggplot2','ggpubr','knitr','tidyverse',
#                'reshape2','tinytex','gt','haven',
#                'dagitty','ellipse','mvtnorm','MASS','splines','gtools',
#                'rethinking','rstan','coda','runjags','rjags', #'loo',
#                'cmdstanr','posterior','bayesplot')
# sapply(librerias, require, character.only=T)
# # sapply(librerias, install.packages, character.only=T)


# set working directory
setwd('~/Desktop/simulations/#final')


# # load set of extra functions
# source("3_functions_extra.R")

# # 0. data ####
# model_dir = 'ListFormat_J100_l0.4_S1.RData'
# load(file.path('data', model_dir) )
# str(data_post)



# 1. FOLV ####

## 1.1 centered mod ####
mcmc_code <- "
data{
    // evaluation
    int N;
    int J;
    int K;
    int L;
    int D;
    int IDj[N];
    int IDk[N];
    int IDl[N];
    int IDd[N];
    int GE[N];
    real AG[N];
    int ED[N];
    int XP[N];
    int y[N];
    
    // individuals
    int IDind[J];
    int G[J];
    real A[J];
    int E[J];
    int X[J];
 
    // items
    int IDitem[K];
    int IDtext[K];
    int IDdim[K];
}
parameters{
    // items
    real m_b[L];
    real<lower=0> s_b[L];
    real b_k[K];
    
    // abilities
    corr_matrix[D] Rho_theta_sub; // sub-dimensions: lit, inf, ref
    vector[D] theta_sub[J]; 
}
model{
    // declare
    vector[D] m_mult[J];
    real v;
    real p;
              
    // items
    m_b ~ normal(0, 1);
    s_b ~ exponential(2);
    for(k in 1:K){            // priors
      b_k[k] ~ normal( m_b[ IDtext[k] ], s_b[ IDtext[k] ]); 
    }
    
    // abilities
    Rho_theta_sub ~ lkj_corr(2);
    theta_sub ~ multi_normal( rep_vector(0, D), Rho_theta_sub );
    
    // model
    for( i in 1:N ) {
      v = theta_sub[ IDj[i], IDd[i] ] - b_k[ IDk[i] ];
      p = inv_logit(v);
      y[i] ~ bernoulli(p);
    }
}
//# generated quantities{
//#     vector[N] log_lik;
//#     real v;
//#     real p;
//#     
//#     // likelihood
//#     for( i in 1:N ) {
//#       v = theta_sub[ IDj[i], IDd[i] ] - b_k[ IDk[i] ];
//#       p = inv_logit(v);
//#       log_lik[i] = bernoulli_lpmf( y[i] | p);
//#     }
//# }
"
save_code = "FOLV_CE_mod.stan"
writeLines(mcmc_code, con=file.path(getwd(), 'models_post', save_code) )

# # (just to check)
# set_cmdstan_path('~/cmdstan')
# mod = cmdstan_model( file.path(getwd(), 'models_post', save_code) )
# mod$sample(data=data_post, output_dir=getwd(),
#            chains=1, parallel_chains=1, init=0, adapt_delta=0.99 )




## 1.2 centered ####

mcmc_code <- "
data{
    // evaluation
    int N;
    int J;
    int K;
    int L;
    int D;
    int IDj[N];
    int IDk[N];
    int IDl[N];
    int IDd[N];
    int GE[N];
    real AG[N];
    int ED[N];
    int XP[N];
    int y[N];
    
    // individuals
    int IDind[J];
    int G[J];
    real A[J];
    int E[J];
    int X[J];
 
    // items
    int IDitem[K];
    int IDtext[K];
    int IDdim[K];
}
parameters{
    // items
    real m_b[L];
    real<lower=0> s_b[L];
    real b_k[K];

    // betas
    real a;
    real b_G[2];
    real b_A;
    real b_E[3];
    real b_X[4];
    
    // abilities
    corr_matrix[D] Rho_theta_sub; // sub-dimensions: lit, inf, ref
    vector[D] theta_sub[J]; 
}
model{
    // declare
    vector[D] m_mult[J];
    real v;
    real p;
              
    // items
    m_b ~ normal(0, 1);
    s_b ~ exponential(2);
    for(k in 1:K){            // priors
      b_k[k] ~ normal( m_b[ IDtext[k] ], s_b[ IDtext[k] ]); 
    }
    
    // abilities
    a ~ normal(0, 0.5);
    b_G ~ normal(0, 0.5);
    b_A ~ normal(0, 0.5);
    b_E ~ normal(0, 1);
    b_X ~ normal(0, 0.5);
    Rho_theta_sub ~ lkj_corr(2);
    for(j in 1:J){
      m_mult[j,] = rep_vector( a + b_G[ G[j] ] +
                                b_A * ( A[j] - min(A) ) +
                                b_E[ E[j] ] + 
                                b_X[ X[j] ], D);
      // using the same predictor for the three latents
    }
    theta_sub ~ multi_normal( m_mult, Rho_theta_sub );
    
    // model
    for( i in 1:N ) {
      v = theta_sub[ IDj[i], IDd[i] ] - b_k[ IDk[i] ];
      p = inv_logit(v);
      y[i] ~ bernoulli(p);
    }
}
//# generated quantities{
//#     vector[N] log_lik;
//#     real v;
//#     real p;
//#     
//#     // likelihood
//#     for( i in 1:N ) {
//#       v = theta_sub[ IDj[i], IDd[i] ] - b_k[ IDk[i] ];
//#       p = inv_logit(v);
//#       log_lik[i] = bernoulli_lpmf( y[i] | p);
//#     }
//# }
"
save_code = "FOLV_CE.stan"
writeLines(mcmc_code, con=file.path(getwd(), 'models_post', save_code) )

# # (just to check)
# set_cmdstan_path('~/cmdstan')
# mod = cmdstan_model( file.path(getwd(), 'models_post', save_code) )
# mod$sample(data=data_post, output_dir=getwd(),
#            chains=1, parallel_chains=1, init=0, adapt_delta=0.99 )





## 1.3 non-centered mod ####

mcmc_code <- "
data{
   // evaluation
    int N;
    int J;
    int K;
    int L;
    int D;
    int IDj[N];
    int IDk[N];
    int IDl[N];
    int IDd[N];
    int GE[N];
    real AG[N];
    int ED[N];
    int XP[N];
    int y[N];
    
    // individuals
    int IDind[J];
    int G[J];
    real A[J];
    int E[J];
    int X[J];
 
    // items
    int IDitem[K];
    int IDtext[K];
    int IDdim[K];
}
parameters{
    // items
    real m_b[L];
    real<lower=0> s_b[L];
    real zb_k[K];
    
    // abilities
    cholesky_factor_corr[D] L_Rho_theta_sub; // sub-dimensions: lit, inf, ref
    matrix[D, J] ztheta_sub; 
}
transformed parameters{
    real b_k[K];
    matrix[D, D] Rho_theta_sub;
    matrix[J, D] m_mult;
    matrix[J, D] theta_sub;
    
    // items
    for(k in 1:K){
      b_k[k] = m_b[ IDtext[k] ] + s_b[ IDtext[k] ] * zb_k[k]; 
    }
    
    // abilities
    Rho_theta_sub = multiply_lower_tri_self_transpose(L_Rho_theta_sub);
    theta_sub = (L_Rho_theta_sub * ztheta_sub)';
}
model{
    // declare
    real v;
    real p;
              
    // items
    m_b ~ normal(0, 1);
    s_b ~ exponential(2);
    zb_k ~ normal(0, 1);
    
    // abilities
    L_Rho_theta_sub ~ lkj_corr_cholesky(2);
    to_vector(ztheta_sub) ~ normal(0,1);
    
    // model
    for( i in 1:N ) {
      v = theta_sub[ IDj[i], IDd[i] ] - b_k[ IDk[i] ];
      p = inv_logit(v);
      y[i] ~ bernoulli(p);
    }
}
//# generated quantities{
//#     vector[N] log_lik;
//#     real v;
//#     real p;
//#     
//#     // likelihood
//#     for( i in 1:N ) {
//#       v = theta_sub[ IDj[i], IDd[i] ] - b_k[ IDk[i] ];
//#       p = inv_logit(v);
//#       log_lik[i] = bernoulli_lpmf( y[i] | p);
//#     }
//# }
"
save_code = "FOLV_NC_mod.stan"
writeLines(mcmc_code, con=file.path(getwd(), 'models_post', save_code) )

# # (just to check)
# set_cmdstan_path('~/cmdstan')
# mod = cmdstan_model( file.path(getwd(), 'models_post', save_code) )
# mod$sample(data=data_post, output_dir=getwd(),
#            chains=1, parallel_chains=1, init=0, adapt_delta=0.99 )




## 1.2 non-centered ####

mcmc_code <- "
data{
   // evaluation
    int N;
    int J;
    int K;
    int L;
    int D;
    int IDj[N];
    int IDk[N];
    int IDl[N];
    int IDd[N];
    int GE[N];
    real AG[N];
    int ED[N];
    int XP[N];
    int y[N];
    
    // individuals
    int IDind[J];
    int G[J];
    real A[J];
    int E[J];
    int X[J];
 
    // items
    int IDitem[K];
    int IDtext[K];
    int IDdim[K];
}
parameters{
    // items
    real m_b[L];
    real<lower=0> s_b[L];
    real zb_k[K];

    // betas
    real a;
    real b_G[2];
    real b_A;
    real b_E[3];
    real b_X[4];
    
    // abilities
    cholesky_factor_corr[D] L_Rho_theta_sub; // sub-dimensions: lit, inf, ref
    matrix[D, J] ztheta_sub; 
}
transformed parameters{
    real b_k[K];
    matrix[D, D] Rho_theta_sub;
    matrix[J, D] m_mult;
    matrix[J, D] theta_sub;
    
    // items
    for(k in 1:K){
      b_k[k] = m_b[ IDtext[k] ] + s_b[ IDtext[k] ] * zb_k[k]; 
    }
    
    // abilities
    Rho_theta_sub = multiply_lower_tri_self_transpose(L_Rho_theta_sub);
    for(j in 1:J){
      m_mult[j,] = rep_row_vector( a + b_G[ G[j] ] + 
                                    b_A * (A[j] - min(A) ) +
                                    b_E[ E[j] ] + 
                                    b_X[ X[j] ], D);
    }
    theta_sub = (L_Rho_theta_sub * ztheta_sub)';
    theta_sub = theta_sub + m_mult;
}
model{
    // declare
    real v;
    real p;
              
    // items
    m_b ~ normal(0, 1);
    s_b ~ exponential(2);
    zb_k ~ normal(0, 1);
    
    // abilities
    a ~ normal(0, 0.5);
    b_G ~ normal(0, 0.5);
    b_A ~ normal(0, 0.5);
    b_E ~ normal(0, 1);
    b_X ~ normal(0, 0.5);
    L_Rho_theta_sub ~ lkj_corr_cholesky(2);
    to_vector(ztheta_sub) ~ normal(0,1);
    
    // model
    for( i in 1:N ) {
      v = theta_sub[ IDj[i], IDd[i] ] - b_k[ IDk[i] ];
      p = inv_logit(v);
      y[i] ~ bernoulli(p);
    }
}
//# generated quantities{
//#     vector[N] log_lik;
//#     real v;
//#     real p;
//#     
//#     // likelihood
//#     for( i in 1:N ) {
//#       v = theta_sub[ IDj[i], IDd[i] ] - b_k[ IDk[i] ];
//#       p = inv_logit(v);
//#       log_lik[i] = bernoulli_lpmf( y[i] | p);
//#     }
//# }
"
save_code = "FOLV_NC.stan"
writeLines(mcmc_code, con=file.path(getwd(), 'models_post', save_code) )

# # (just to check)
# set_cmdstan_path('~/cmdstan')
# mod = cmdstan_model( file.path(getwd(), 'models_post', save_code) )
# mod$sample(data=data_post, output_dir=getwd(),
#            chains=1, parallel_chains=1, init=0, adapt_delta=0.99 )





# 2. SOLV ####

## 2.1 centered ####
mcmc_code <- "
data{
    // evaluation
    int N;
    int J;
    int K;
    int L;
    int D;
    int IDj[N];
    int IDk[N];
    int IDl[N];
    int IDd[N];
    int GE[N];
    real AG[N];
    int ED[N];
    int XP[N];
    int y[N];
    
    // individuals
    int IDind[J];
    int G[J];
    real A[J];
    int E[J];
    int X[J];
 
    // items
    int IDitem[K];
    int IDtext[K];
    int IDdim[K];
}
parameters{
    // items
    real m_b[L];
    real<lower=0> s_b[L];
    real b_k[K];

    // betas
    real a;
    real b_G[2];
    real b_A;
    real b_E[3];
    real b_X[4];
    
    // abilities
    vector[J] theta;              // reading comprehension
    real<lower=0> loads[D];       // loadings
    corr_matrix[D] Rho_theta_sub; // sub-dimensions: lit, inf, ref
    vector[D] theta_sub[J]; 
}
model{
    // declare
    vector[J] m_theta;
    vector[D] m_mult[J];
    real v;
    real p;
              
    // items
    m_b ~ normal(0, 1);
    s_b ~ exponential(2);
    for(k in 1:K){            // priors
      b_k[k] ~ normal( m_b[ IDtext[k] ], s_b[ IDtext[k] ]); 
    }
    
    // SOLV
    a ~ normal(0, 0.5);
    b_G ~ normal(0, 0.5);
    b_A ~ normal(0, 0.5);
    b_E ~ normal(0, 1);
    b_X ~ normal(0, 0.5);
    for(j in 1:J){
      m_theta[j] = a + b_G[ G[j] ] + 
                    b_A * ( A[j] - min(A) ) +
                    b_E[ E[j] ] + 
                    b_X[ X[j] ];
    }
    theta ~ normal(m_theta, 1); 
    // setting the scale at this level
    
    // FOLV
    loads ~ lognormal(0, 0.5);
    for(j in 1:J){
      m_mult[j,] = [ loads[1]*theta[j], 
                      loads[2]*theta[j], 
                      loads[3]*theta[j] ]';
    }
    Rho_theta_sub ~ lkj_corr(2);
    theta_sub ~ multi_normal( m_mult, Rho_theta_sub ); 
    // also setting scale here
    
    // model
    for( i in 1:N ) {
      v = theta_sub[ IDj[i], IDd[i] ] - b_k[ IDk[i] ];
      p = inv_logit(v);
      y[i] ~ bernoulli(p);
    }
}
//# generated quantities{
//#     vector[N] log_lik;
//#     real v;
//#     real p;
//#     
//#     // likelihood
//#     for( i in 1:N ) {
//#       v = theta_sub[ IDj[i], IDd[i] ] - b_k[ IDk[i] ];
//#       p = inv_logit(v);
//#       log_lik[i] = bernoulli_lpmf( y[i] | p);
//#     }
//# }
"
save_code = "SOLV_CE.stan"
writeLines(mcmc_code, con=file.path(getwd(), 'models_post', save_code) )

# # (just to check)
# set_cmdstan_path('~/cmdstan')
# mod = cmdstan_model( file.path(getwd(), 'models_post', save_code) )
# mod$sample(data=data_post, output_dir=getwd(),
#            chains=1, parallel_chains=1, init=0, adapt_delta=0.99)




## 2.1 non-centered ####
mcmc_code <- "
data{
    // evaluation
    int N;
    int J;
    int K;
    int L;
    int D;
    int IDj[N];
    int IDk[N];
    int IDl[N];
    int IDd[N];
    int GE[N];
    real AG[N];
    int ED[N];
    int XP[N];
    int y[N];
    
    // individuals
    int IDind[J];
    int G[J];
    real A[J];
    int E[J];
    int X[J];
 
    // items
    int IDitem[K];
    int IDtext[K];
    int IDdim[K];
}
parameters{
    // items
    real m_b[L];
    real<lower=0> s_b[L];
    real zb_k[K];

    // betas
    real a;
    real b_G[2];
    real b_A;
    real b_E[3];
    real b_X[4];
    
    // abilities
    vector[J] ztheta;             // reading comprehension
    real<lower=0> loads[D];       // loadings
    cholesky_factor_corr[D] L_Rho_theta_sub; // sub-dimensions: lit, inf, ref
    matrix[D, J] ztheta_sub; 
}
transformed parameters{
    real b_k[K];
    vector[J] m_theta;
    vector[J] theta;             // reading comprehension
    matrix[D, D] Rho_theta_sub;
    matrix[J, D] m_mult;
    matrix[J, D] theta_sub;
    
    // items
    for(k in 1:K){
      b_k[k] = m_b[ IDtext[k] ] + s_b[ IDtext[k] ] * zb_k[k]; 
    }
    
    // SOLV
    for(j in 1:J){
      m_theta[j] = a + b_G[ G[j] ] + 
                    b_A * ( A[j] - min(A) ) +
                    b_E[ E[j] ] + 
                    b_X[ X[j] ];
    }
    theta = m_theta + 1 * ztheta;
    // setting the scale at this level
    
    // FOLV
    Rho_theta_sub = multiply_lower_tri_self_transpose(L_Rho_theta_sub);
    for(j in 1:J){
      m_mult[j,] = to_row_vector( [ loads[1]*theta[j], 
                                    loads[2]*theta[j], 
                                    loads[3]*theta[j] ] );
    }
    theta_sub = m_mult + (L_Rho_theta_sub * ztheta_sub)';
    // also setting scale here
}
model{
    // declare
    real v;
    real p;
              
    // items
    m_b ~ normal(0, 1);
    s_b ~ exponential(2);
    zb_k ~ normal(0, 1);
    
    // abilities
    a ~ normal(0, 0.5);
    b_G ~ normal(0, 0.5);
    b_A ~ normal(0, 0.5);
    b_E ~ normal(0, 1);
    b_X ~ normal(0, 0.5);
    ztheta ~ normal(0,1);
    L_Rho_theta_sub ~ lkj_corr_cholesky(2);
    loads ~ lognormal(0, 0.5);
    to_vector(ztheta_sub) ~ normal(0,1);
    
    // model
    for( i in 1:N ) {
      v = theta_sub[ IDj[i], IDd[i] ] - b_k[ IDk[i] ];
      p = inv_logit(v);
      y[i] ~ bernoulli(p);
    }
}
//# generated quantities{
//#     vector[N] log_lik;
//#     real v;
//#     real p;
//#     
//#     // likelihood
//#     for( i in 1:N ) {
//#       v = theta_sub[ IDj[i], IDd[i] ] - b_k[ IDk[i] ];
//#       p = inv_logit(v);
//#       log_lik[i] = bernoulli_lpmf( y[i] | p);
//#     }
//# }
"
save_code = "SOLV_NC.stan"
writeLines(mcmc_code, con=file.path(getwd(), 'models_post', save_code) )

# # (just to check)
# set_cmdstan_path('~/cmdstan')
# mod = cmdstan_model( file.path(getwd(), 'models_post', save_code) )
# mod$sample(data=data_post, output_dir=getwd(),
#            chains=1, parallel_chains=1, init=0, adapt_delta=0.99)



# 3. multilevel ####

## 3.1 centered ####
mcmc_code <- "
data{
    // evaluation
    int N;
    int J;
    int K;
    int L;
    int D;
    int IDj[N];
    int IDk[N];
    int IDl[N];
    int IDd[N];
    int GE[N];
    real AG[N];
    int ED[N];
    int XP[N];
    int y[N];
    
    // individuals
    int IDind[J];
    int G[J];
    real A[J];
    int E[J];
    int X[J];
 
    // items
    int IDitem[K];
    int IDtext[K];
    int IDdim[K];
}
parameters{
    // individuals
    real<lower=0> s_ind;
    vector[J] ind_j;
    real<lower=0> s_dim;
    vector[D] dim_d;
    
    // items and text
    real<lower=0> s_item;   
    vector[K] item_k;
    real<lower=0> s_text;
    vector[L] text_l;

    // betas
    real a;
    real b_G[2];
    real b_A;
    real b_E[3];
    real b_X[4];
}
model{
    // declare
    real v;
    real p; 
    
    // abilities
    s_ind ~ exponential(2);
    ind_j ~ normal(0, s_ind);
    s_dim ~ exponential(2);
    dim_d ~ normal(0, s_dim);
    
    a ~ normal(0, 0.5);
    b_G ~ normal(0, 0.5);
    b_A ~ normal(0, 0.5);
    b_E ~ normal(0, 1);
    b_X ~ normal(0, 0.5);
    
    // items and texts
    s_item ~ exponential(2);
    item_k ~ normal(0, s_item);
    s_text ~ exponential(2);
    text_l ~ normal(0, s_text);

    // model
    for( i in 1:N ) {
      v = ( ind_j[ IDj[i] ] + dim_d[ IDd[i] ] ) -
          ( item_k[ IDk[i] ] + text_l[ IDl[i] ] ) +
          a + b_G[ GE[i] ] +
                b_A * ( AG[i] - min(A) ) +
                b_E[ ED[i] ] +
                b_X[ XP[i] ];
      p = inv_logit(v);
      p = inv_logit(v);
      y[i] ~ bernoulli(p);
    }
}
//# generated quantities{
//#     vector[N] log_lik;
//#     real v;
//#     real p;
//# 
//#     // likelihood
//#     for( i in 1:N ) {
//#       v = ( ind_j[ IDj[i] ] + dim_d[ IDd[i] ] ) -
//#           ( item_k[ IDk[i] ] + text_l[ IDl[i] ] ) +
//#           a + b_G[ GE[i] ] +
//#                 b_A * ( AG[i] - min(A) ) +
//#                 b_E[ ED[i] ] +
//#                 b_X[ XP[i] ];
//#       p = inv_logit(v);
//#       log_lik[i] = bernoulli_lpmf( y[i] | p);
//#     }
//# }
"
save_code = "MULT_CE.stan"
writeLines(mcmc_code, con=file.path(getwd(), 'models_post', save_code) )

# # (just to check)
# set_cmdstan_path('~/cmdstan')
# mod = cmdstan_model( file.path(getwd(), 'models_post', save_code) )
# mod$sample(data=data_post, output_dir=getwd(),
#            chains=1, parallel_chains=1, init=0, adapt_delta=0.99)




## 3.2 non-centered ####
mcmc_code <- "
data{
    // evaluation
    int N;
    int J;
    int K;
    int L;
    int D;
    int IDj[N];
    int IDk[N];
    int IDl[N];
    int IDd[N];
    int GE[N];
    real AG[N];
    int ED[N];
    int XP[N];
    int y[N];
    
    // individuals
    int IDind[J];
    int G[J];
    real A[J];
    int E[J];
    int X[J];
 
    // items
    int IDitem[K];
    int IDtext[K];
    int IDdim[K];
}
parameters{
    // individuals
    real<lower=0> s_ind;
    vector[J] zind_j;
    real<lower=0> s_dim;
    vector[D] zdim_d;
    
    // items and text
    real<lower=0> s_item;   
    vector[K] zitem_k;
    real<lower=0> s_text;
    vector[L] ztext_l;

    // betas
    real a;
    real b_G[2];
    real b_A;
    real b_E[3];
    real b_X[4];
}
transformed parameters{
    vector[J] m_ind;
    vector[J] ind_j;
    vector[D] dim_d;
    vector[K] item_k;
    vector[L] text_l;
    
    // individuals
    ind_j = s_ind * zind_j;
    dim_d = s_dim * zdim_d;
    
    // items and texts
    item_k = s_item * zitem_k;
    text_l = s_text * ztext_l;
}
model{
    // declare
    real v;
    real p; 
    
    // abilities
    s_ind ~ exponential(2);
    zind_j ~ normal(0, 1);
    s_dim ~ exponential(2);
    zdim_d ~ normal(0, 1);
    
    a ~ normal(0, 0.5);
    b_G ~ normal(0, 0.5);
    b_A ~ normal(0, 0.5);
    b_E ~ normal(0, 1);
    b_X ~ normal(0, 0.5);
    
    // items and texts
    s_item ~ exponential(2);
    zitem_k ~ normal(0, 1);
    s_text ~ exponential(2);
    ztext_l ~ normal(0, 1);

    // model
    for( i in 1:N ) {
      v = ( ind_j[ IDj[i] ] + dim_d[ IDd[i] ] ) -
          ( item_k[ IDk[i] ] + text_l[ IDl[i] ] ) +
          a + b_G[ GE[i] ] +
                b_A * ( AG[i] - min(A) ) +
                b_E[ ED[i] ] +
                b_X[ XP[i] ];
      p = inv_logit(v);
      y[i] ~ bernoulli(p);
    }
}
//# generated quantities{
//#     vector[N] log_lik;
//#     real v;
//#     real p;
//# 
//#     // likelihood
//#     for( i in 1:N ) {
//#       v = ( ind_j[ IDj[i] ] + dim_d[ IDd[i] ] ) -
//#           ( item_k[ IDk[i] ] + text_l[ IDl[i] ] ) +
//#           a + b_G[ GE[i] ] +
//#                 b_A * ( AG[i] - min(A) ) +
//#                 b_E[ ED[i] ] +
//#                 b_X[ XP[i] ];
//#       p = inv_logit(v);
//#       log_lik[i] = bernoulli_lpmf( y[i] | p);
//#     }
//# }
"
save_code = "MULT_NC.stan"
writeLines(mcmc_code, con=file.path(getwd(), 'models_post', save_code) )

# # (just to check)
# set_cmdstan_path('~/cmdstan')
# mod = cmdstan_model( file.path(getwd(), 'models_post', save_code) )
# fit = mod$sample(data=data_post, output_dir=getwd(),
#                  chains=1, parallel_chains=1, init=0, adapt_delta=0.99)




