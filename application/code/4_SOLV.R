# preliminar ####

# erase all objects
rm(list=ls())

# load libraries
librerias <- c('stringr','dplyr','ggplot2','ggpubr','knitr','tidyverse',
               'reshape2','tinytex','gt','haven',
               'dagitty','ellipse','mvtnorm','MASS','splines','gtools',
               'rethinking','rstan','coda','runjags','rjags', #'loo',
               'cmdstanr','posterior','bayesplot')
sapply(librerias, require, character.only=T)
# sapply(librerias, install.packages, character.only=T)

# set working directory
setwd('~/Desktop/application')

# load set of extra functions
source(file.path(getwd(),'# code','2_functions_extra.R'))

# saving par
opar = par()



# 1. data ####
load( file=file.path(getwd(), 'data', '4_final', 'data_application.RData') )
# str(data_app)



# 2. models ####

## 2.1 centered version ####
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
    int G[N];
    real A[N];
    int E[N];
    int S[N];
    int Xpu[N];
    int Xpr[N];
    int Di[N];
    int y[N];
    
    // individuals
    int IDind[J];
    int GE[J];
    real AG[J];
    int ED[J];
    int SP[J];
    int XPpu[J];
    int XPpr[J];
    int DI[J];
 
    // items
    int IDitem[K];
    int IDtext[K];
    int IDdim[K];
} 
transformed data {
  real min_AG;
  min_AG = min(AG);
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
    real b_S[3];
    real b_Xpu[4];
    real b_Xpr[4];
    real b_Di[4];
    
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
    b_A~ normal(0, 0.5);
    b_E ~ normal(0, 1);
    b_S ~ normal(0, 0.5);
    b_Xpu ~ normal(0, 0.5);
    b_Xpr ~ normal(0, 0.5);
    b_Di ~ normal(0, 0.5);
    for(j in 1:J){
      m_theta[j] = a + b_G[ GE[j] ] + b_A*(AG[j] - min_AG) +
                    b_E[ ED[j] ] + b_S[ SP[j] ] + b_Di[ DI[j] ] +
                    b_Xpu[ XPpu[j] ] + b_Xpr[ XPpr[j] ];
    }
    theta ~ normal(m_theta, 1); 
    // setting the scale at this level
    
    // FOLV
    loads ~ lognormal(0, 0.5);
    for(j in 1:J){
      m_mult[j,] = [ loads[1]*theta[j], loads[2]*theta[j], loads[3]*theta[j] ]';
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
generated quantities{
    vector[N] log_lik;
    real v;
    real p;

    // likelihood
    for( i in 1:N ) {
      v = theta_sub[ IDj[i], IDd[i] ] - b_k[ IDk[i] ];
      p = inv_logit(v);
      log_lik[i] = bernoulli_lpmf( y[i] | p);
    }
}
"
save_code = file.path(getwd(), 'results', "SOLV_center.stan")
writeLines(mcmc_code, con=save_code)


# # (you need to run it the first time, then hide the following code)
# set_cmdstan_path('~/cmdstan')
# mod = cmdstan_model(save_code)
# 
# t0 = proc.time()
# mod$sample(data=data_app, output_dir=file.path(getwd(),'results'),
#            output_basename = 'SOLV_center',
#            chains=3, parallel_chains=3, init=0, adapt_delta=0.99 )
# t1 = proc.time()
# 
# elapsed = (t1-t0)
# save(elapsed, file=file.path(getwd(), 'results', 'SOLV_center_time.RData'))




## 2.2 non-centered version ####
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
    int G[N];
    real A[N];
    int E[N];
    int S[N];
    int Xpu[N];
    int Xpr[N];
    int Di[N];
    int y[N];
    
    // individuals
    int IDind[J];
    int GE[J];
    real AG[J];
    int ED[J];
    int SP[J];
    int XPpu[J];
    int XPpr[J];
    int DI[J];
 
    // items
    int IDitem[K];
    int IDtext[K];
    int IDdim[K];
} 
transformed data {
  real min_AG;
  min_AG = min(AG);
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
    real b_S[3];
    real b_Xpu[4];
    real b_Xpr[4];
    real b_Di[4];
    
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
      m_theta[j] = a + b_G[ GE[j] ] + b_A*(AG[j] - min_AG) +
                    b_E[ ED[j] ] + b_S[ SP[j] ] + b_Di[ DI[j] ] +
                    b_Xpu[ XPpu[j] ] + b_Xpr[ XPpr[j] ];
    }
    theta = m_theta + 1 * ztheta;
    // setting the scale at this level
    
    // FOLV
    Rho_theta_sub = multiply_lower_tri_self_transpose(L_Rho_theta_sub);
    for(j in 1:J){
      m_mult[j,] = to_row_vector( [ loads[1]*theta[j], loads[2]*theta[j], loads[3]*theta[j] ] );
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
    b_A~ normal(0, 0.5);
    b_E ~ normal(0, 1);
    b_S ~ normal(0, 0.5);
    b_Xpu ~ normal(0, 0.5);
    b_Xpr ~ normal(0, 0.5);
    b_Di ~ normal(0, 0.5);
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
generated quantities{
    vector[N] log_lik;
    real v;
    real p;

    // likelihood
    for( i in 1:N ) {
      v = theta_sub[ IDj[i], IDd[i] ] - b_k[ IDk[i] ];
      p = inv_logit(v);
      log_lik[i] = bernoulli_lpmf( y[i] | p);
    }
}
"
save_code = file.path(getwd(), 'results', "SOLV_noncenter.stan")
writeLines(mcmc_code, con=save_code)

# # (you need to run it the first time, then hide the following code)
# set_cmdstan_path('~/cmdstan')
# mod = cmdstan_model(save_code)
# 
# t0 = proc.time()
# mod$sample(data=data_app, output_dir=file.path(getwd(),'results'),
#            output_basename = 'SOLV_noncenter',
#            chains=3, parallel_chains=3, init=0, adapt_delta=0.99 )
# t1 = proc.time()
# 
# elapsed = (t1-t0)
# save(elapsed, file=file.path(getwd(), 'results', 'SOLV_noncenter_time.RData'))




# 3. stats ####

## 3.1 centered ####

fit_files = file.path(getwd(), 'results', paste0('SOLV_center-',1:3,'.csv'))
stan_model_c = rstan::read_stan_csv(fit_files)
# precis(stan_model_c, depth=3)

file_save = file.path(getwd(), 'tables')

# stats
result = precis(stan_model_c, depth=4)
save( result, file=file.path(file_save, 'SOLV_CE_stat.RData' ) )

# posterior samples
set.seed(45589)
result = extract.samples( stan_model_c )
save( result, file=file.path(file_save, 'SOLV_CE_post.RData' ) )


# WAIC and PSIS
result = WAIC(stan_model_c, pointwise=TRUE)
save( result, file=file.path(file_save, 'SOLV_CE_WAIC.RData' ) )


fit_files = file.path(getwd(), 'results', paste0('SOLV_center-',1,'.csv'))
stan_model_c = rstan::read_stan_csv(fit_files)
result = PSIS(stan_model_c, pointwise=TRUE)
save( result, file=file.path(file_save, 'SOLV_CE_PSIS_OneChain.RData' ) )




## 3.2 non-centered ####

fit_files = file.path(getwd(), 'results', paste0('SOLV_noncenter-',1:3,'.csv'))
stan_model_nc = rstan::read_stan_csv(fit_files)
# precis(stan_model_nc, depth=3)

file_save = file.path(getwd(), 'tables')

# stats
result = precis(stan_model_nc, depth=4)
save( result, file=file.path(file_save, 'SOLV_NC_stat.RData' ) )

# posterior samples
set.seed(45589)
result = extract.samples( stan_model_nc )
save( result, file=file.path(file_save, 'SOLV_NC_post.RData' ) )


# WAIC and PSIS
result = WAIC(stan_model_nc, pointwise=TRUE)
save( result, file=file.path(file_save, 'SOLV_NC_WAIC.RData' ) )


fit_files = file.path(getwd(), 'results', paste0('SOLV_noncenter-',1,'.csv'))
stan_model_nc = rstan::read_stan_csv(fit_files)
result = PSIS(stan_model_nc, pointwise=TRUE)
save( result, file=file.path(file_save, 'SOLV_NC_PSIS_OneChain.RData' ) )








# 4. assessing chains ####

## 4.1 centered ####

fit_files = file.path(getwd(), 'results', paste0('SOLV_center-',1:3,'.csv'))
stan_model_c = rstan::read_stan_csv(fit_files)
# precis(stan_model_c, depth=3)

file_save = file.path(getwd(), 'figures')


png(file.path(file_save, 'SOLV_CE_mb.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('m_b[',1:5,']')
tri_plot(stan_model_c, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_CE_sb.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('s_b[',1:5,']')
tri_plot(stan_model_c, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_CE_bk_1_5.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('b_k[',1:5,']')
tri_plot(stan_model_c, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_CE_bk_6_10.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('b_k[',6:10,']')
tri_plot(stan_model_c, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_CE_bk_11_15.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('b_k[',11:15,']')
tri_plot(stan_model_c, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_CE_bk_16_20.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('b_k[',16:20,']')
tri_plot(stan_model_c, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_CE_bk_21_25.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('b_k[',21:25,']')
tri_plot(stan_model_c, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_CE_reg_1.png'), 
    units='cm', width=30, height=30, res=100)
idx = c('a', paste0('b_G[',1:2,']'),'b_A')
tri_plot(stan_model_c, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_CE_reg_2.png'), 
    units='cm', width=30, height=30, res=100)
idx = c(paste0('b_E[',1:3,']'))
tri_plot(stan_model_c, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_CE_reg_3.png'), 
    units='cm', width=30, height=30, res=100)
idx = c(paste0('b_S[',1:3,']'))
tri_plot(stan_model_c, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_CE_reg_4.png'), 
    units='cm', width=30, height=30, res=100)
idx = c(paste0('b_Xpu[',1:4,']'))
tri_plot(stan_model_c, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_CE_reg_5.png'), 
    units='cm', width=30, height=30, res=100)
idx = c(paste0('b_Xpr[',1:4,']'))
tri_plot(stan_model_c, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_CE_reg_6.png'), 
    units='cm', width=30, height=30, res=100)
idx = c(paste0('b_Di[',1:4,']'))
tri_plot(stan_model_c, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_CE_rho.png'), 
    units='cm', width=30, height=30, res=100)
idx = c('Rho_theta_sub[1,2]','Rho_theta_sub[1,3]','Rho_theta_sub[2,3]')
tri_plot(stan_model_c, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_CE_loads.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('loads[',1:3,']')
tri_plot(stan_model_c, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_CE_theta1.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('theta_sub[',1:5,',1]')
tri_plot(stan_model_c, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_CE_theta2.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('theta_sub[',1:5,',2]')
tri_plot(stan_model_c, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_CE_theta3.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('theta_sub[',1:5,',3]')
tri_plot(stan_model_c, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_CE_theta.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('theta[',1:5,']')
tri_plot(stan_model_c, pars=idx)
dev.off()



## 4.1 non-centered ####

fit_files = file.path(getwd(), 'results', paste0('SOLV_noncenter-',1:3,'.csv'))
stan_model_nc = rstan::read_stan_csv(fit_files)
# precis(stan_model_nc, depth=3)

file_save = file.path(getwd(), 'figures')


png(file.path(file_save, 'SOLV_NC_mb.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('m_b[',1:5,']')
tri_plot(stan_model_nc, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_NC_sb.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('s_b[',1:5,']')
tri_plot(stan_model_nc, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_NC_bk_1_5.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('b_k[',1:5,']')
tri_plot(stan_model_nc, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_NC_bk_6_10.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('b_k[',6:10,']')
tri_plot(stan_model_nc, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_NC_bk_11_15.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('b_k[',11:15,']')
tri_plot(stan_model_nc, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_NC_bk_16_20.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('b_k[',16:20,']')
tri_plot(stan_model_nc, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_NC_bk_21_25.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('b_k[',21:25,']')
tri_plot(stan_model_nc, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_NC_reg_1.png'), 
    units='cm', width=30, height=30, res=100)
idx = c('a', paste0('b_G[',1:2,']'),'b_A')
tri_plot(stan_model_nc, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_NC_reg_2.png'), 
    units='cm', width=30, height=30, res=100)
idx = c(paste0('b_E[',1:3,']'))
tri_plot(stan_model_nc, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_NC_reg_3.png'), 
    units='cm', width=30, height=30, res=100)
idx = c(paste0('b_S[',1:3,']'))
tri_plot(stan_model_nc, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_NC_reg_4.png'), 
    units='cm', width=30, height=30, res=100)
idx = c(paste0('b_Xpu[',1:4,']'))
tri_plot(stan_model_nc, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_NC_reg_5.png'), 
    units='cm', width=30, height=30, res=100)
idx = c(paste0('b_Xpr[',1:4,']'))
tri_plot(stan_model_nc, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_NC_reg_6.png'), 
    units='cm', width=30, height=30, res=100)
idx = c(paste0('b_Di[',1:4,']'))
tri_plot(stan_model_nc, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_NC_rho.png'), 
    units='cm', width=30, height=30, res=100)
idx = c('Rho_theta_sub[1,2]','Rho_theta_sub[1,3]','Rho_theta_sub[2,3]')
tri_plot(stan_model_nc, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_NC_loads.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('loads[',1:3,']')
tri_plot(stan_model_nc, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_NC_theta1.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('theta_sub[',1:5,',1]')
tri_plot(stan_model_nc, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_NC_theta2.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('theta_sub[',1:5,',2]')
tri_plot(stan_model_nc, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_NC_theta3.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('theta_sub[',1:5,',3]')
tri_plot(stan_model_nc, pars=idx)
dev.off()

png(file.path(file_save, 'SOLV_NC_theta.png'), 
    units='cm', width=30, height=30, res=100)
idx = paste0('theta[',1:5,']')
tri_plot(stan_model_nc, pars=idx)
dev.off()



# 5. recovery tables and plots ####

file_save = file.path(getwd(), 'tables')

## 5.1 load ####

# load( file.path(getwd(), 'tables', 'SOLV_CE_stat.RData') )
# result_CE = result
# result_CE$parameter = row.names(result_CE)
# 
# load( file.path(getwd(), 'tables', 'SOLV_NC_stat.RData') )
# result_NC = result
# result_NC$parameter = row.names(result_NC)
# 
# stan_result = merge(result_CE, result_NC, by='parameter', all.x=T, all.y=F)
# save( stan_result, file=file.path(file_save, 'SOLV_stat.RData' ) )
# rm(list=c('result','result_CE','result_NC'))

load( file.path(file_save, 'SOLV_stat.RData' ) )


# statistic of interest
stat_table = stan_result[,c("parameter",
                            'n_eff.x','n_eff.y',
                            'Rhat4.x','Rhat4.y')]



## 5.2 stat plots ####
file_save = file.path(getwd(), 'figures')

### items  and texts ####
png(file.path(file_save, 'SOLV_stat_mb.png'), 
    units='cm', width=30, height=12, res=100)
plot_stat_both(stat_table=stat_table, par_int=paste0('m_b[',1:5,']'))
dev.off()

png(file.path(file_save, 'SOLV_stat_sb.png'), 
    units='cm', width=30, height=12, res=100)
plot_stat_both(stat_table=stat_table, par_int=paste0('m_b[',1:5,']'))
dev.off()

png(file.path(file_save, 'SOLV_stat_bk.png'), 
    units='cm', width=30, height=12, res=100)
plot_stat_both(stat_table=stat_table, par_int=paste0('b_k[',1:25,']'))
dev.off()


### regression ####
png(file.path(file_save, 'SOLV_stat_a.png'), 
    units='cm', width=30, height=12, res=100)
plot_stat_both(stat_table=stat_table, par_int='a' )
dev.off()

png(file.path(file_save, 'SOLV_stat_bG.png'), 
    units='cm', width=30, height=12, res=100)
plot_stat_both(stat_table=stat_table, par_int=paste0('b_G[',1:2,']'))
dev.off()

png(file.path(file_save, 'SOLV_stat_bA.png'), 
    units='cm', width=30, height=12, res=100)
plot_stat_both(stat_table=stat_table, par_int='b_A' )
dev.off()

png(file.path(file_save, 'SOLV_stat_bE.png'), 
    units='cm', width=30, height=12, res=100)
plot_stat_both(stat_table=stat_table, par_int=paste0('b_E[',1:3,']'))
dev.off()

png(file.path(file_save, 'SOLV_stat_bS.png'), 
    units='cm', width=30, height=12, res=100)
plot_stat_both(stat_table=stat_table, par_int=paste0('b_S[',1:4,']'))
dev.off()

png(file.path(file_save, 'SOLV_stat_bX.png'), 
    units='cm', width=30, height=12, res=100)
idx_var = c( paste0('b_Xpu[',1:4,']'), paste0('b_Xpr[',1:4,']'))
plot_stat_both(stat_table=stat_table, par_int=idx_var)
dev.off()

png(file.path(file_save, 'SOLV_stat_bDi.png'), 
    units='cm', width=30, height=12, res=100)
plot_stat_both(stat_table=stat_table, par_int=paste0('b_Di[',1:4,']'))
dev.off()

png(file.path(file_save, 'SOLV_stat_rho.png'), 
    units='cm', width=30, height=12, res=100)
idx_var = c('Rho_theta_sub[1,2]','Rho_theta_sub[1,3]','Rho_theta_sub[2,3]')
plot_stat_both(stat_table=stat_table, par_int=idx_var)
dev.off()

png(file.path(file_save, 'SOLV_stat_loads.png'), 
    units='cm', width=30, height=12, res=100)
idx_var = paste0('loads[',1:3,']')
plot_stat_both(stat_table=stat_table, par_int=idx_var)
dev.off()

### abilities ####

png(file.path(file_save, 'SOLV_stat_theta_sub.png'), 
    units='cm', width=30, height=12, res=100)
idx_var = c(paste0('theta_sub[',1:2000,',1]'), 
            paste0('theta_sub[',1:2000,',2]'),
            paste0('theta_sub[',1:2000,',3]'))
plot_stat_both(stat_table=stat_table, par_int=idx_var)
dev.off()




## 5.3 recovery plots ####

file_save = file.path(getwd(), 'figures')


### items and texts ####
png(file.path(file_save, 'SOLV_recovery_texts.png'), 
    units='cm', width=30, height=12, res=100)
par(mfrow=c(1, 2))
recover_plot(result_object = stan_result, title='Texts difficulties',
             par_int = c( paste0('m_b[',1:5,']')) )
recover_plot(result_object = stan_result, title='Texts deviations',
             par_int = c( paste0('s_b[',1:5,']')) )
par(mfrow=c(1, 1))
dev.off()


png(file.path(file_save, 'SOLV_recovery_items.png'), 
    units='cm', width=30, height=12, res=100)
par(mfrow=c(1, 1))
recover_plot(result_object = stan_result, title='Items',
             par_int = paste0('b_k[',1:25,']'))
par(mfrow=c(1, 1))
dev.off()


### regression ####

png(file.path(file_save, 'SOLV_recovery_reg1.png'), 
    units='cm', width=15, height=12, res=100)
par(mfrow=c(1, 1))
recover_plot(result_object = stan_result, title='Intercept and age', 
             par_int=c('a', 'b_A'))
par(mfrow=c(1, 1))
dev.off()

png(file.path(file_save, 'SOLV_recovery_reg2.png'), 
    units='cm', width=15, height=12, res=100)
par(mfrow=c(1, 1))
recover_plot(result_object = stan_result, title='Gender', 
             par_int = paste0('b_G[',1:2,']'))
par(mfrow=c(1, 1))
dev.off()

png(file.path(file_save, 'SOLV_recovery_reg3.png'), 
    units='cm', width=15, height=12, res=100)
par(mfrow=c(1, 1))
recover_plot(result_object = stan_result, title='Education', 
             par_int = paste0('b_E[',1:3,']'))
par(mfrow=c(1, 1))
dev.off()

png(file.path(file_save, 'SOLV_recovery_reg4.png'), 
    units='cm', width=15, height=12, res=100)
par(mfrow=c(1, 1))
recover_plot(result_object = stan_result, title='Specialty', 
             par_int = paste0('b_S[',1:4,']'))
par(mfrow=c(1, 1))
dev.off()

png(file.path(file_save, 'SOLV_recovery_reg5.png'), 
    units='cm', width=15, height=12, res=100)
par(mfrow=c(1, 1))
recover_plot(result_object = stan_result, title='Public experience', 
             par_int = paste0('b_Xpu[',1:4,']'))
par(mfrow=c(1, 1))
dev.off()

png(file.path(file_save, 'SOLV_recovery_reg6.png'), 
    units='cm', width=15, height=12, res=100)
par(mfrow=c(1, 1))
recover_plot(result_object = stan_result, title='Private experience', 
             par_int = paste0('b_Xpr[',1:4,']'))
par(mfrow=c(1, 1))
dev.off()

png(file.path(file_save, 'SOLV_recovery_reg7.png'), 
    units='cm', width=15, height=12, res=100)
par(mfrow=c(1, 1))
recover_plot(result_object = stan_result, title='Disability', 
             par_int = paste0('b_Di[',1:4,']'))
par(mfrow=c(1, 1))
dev.off()

png(file.path(file_save, 'SOLV_recovery_rho.png'), 
    units='cm', width=17, height=12, res=100)
par(mfrow=c(1, 1), mar=c(9,4.1,4.1,2.1))
idx_var = c('Rho_theta_sub[1,2]','Rho_theta_sub[1,3]','Rho_theta_sub[2,3]')
recover_plot(result_object = stan_result, title='Correlation', 
             par_int = idx_var)
par(mfrow=c(1, 1), mar=opar$mar)
dev.off()

png(file.path(file_save, 'SOLV_recovery_loads.png'), 
    units='cm', width=17, height=12, res=100)
par(mfrow=c(1, 1))
idx_var = paste0('loads[',1:3,']')
recover_plot(result_object = stan_result, title='Loadings', 
             par_int = idx_var)
par(mfrow=c(1, 1))
dev.off()


### abilities ####
set.seed(4668)
ind = sample(1:2000, size=50)

png(file.path(file_save, 'SOLV_recovery_ability1.png'), 
    units='cm', width=30, height=15, res=100)
par(mfrow=c(1, 1), mar=c(8,4.1,4.1,2.1))
recover_plot(result_object = stan_result, title='Literal',
             par_int = paste0('theta_sub[',ind,',1]'))
par(mfrow=c(1, 1), mar=opar$mar)
dev.off()

png(file.path(file_save, 'SOLV_recovery_ability2.png'), 
    units='cm', width=30, height=15, res=100)
par(mfrow=c(1, 1), mar=c(8,4.1,4.1,2.1))
recover_plot(result_object = stan_result, title='Inferential',
             par_int = paste0('theta_sub[',ind,',2]'))
par(mfrow=c(1, 1), mar=opar$mar)
dev.off()

png(file.path(file_save, 'SOLV_recovery_ability3.png'), 
    units='cm', width=30, height=15, res=100)
par(mfrow=c(1, 1), mar=c(8,4.1,4.1,2.1))
recover_plot(result_object = stan_result, title='Reflective',
             par_int = paste0('theta_sub[',ind,',3]'))
par(mfrow=c(1, 1), mar=opar$mar)
dev.off()

png(file.path(file_save, 'SOLV_recovery_ability.png'), 
    units='cm', width=30, height=15, res=100)
par(mfrow=c(1, 1), mar=c(6,4.1,4.1,2.1))
recover_plot(result_object = stan_result, title='Reading comprehension',
             par_int = paste0('theta[',ind,']'))
par(mfrow=c(1, 1), mar=opar$mar)
dev.off()



### contrasts ####

# calculations
file_save = file.path(getwd(), 'tables')

# load( file.path(file_save, 'SOLV_stat.RData' ) )
# 
# # centered
# load( file.path(file_save, 'SOLV_CE_post.RData' ) )
# contrast_CE = contrast_stat(results_object=stan_result,
#                             post=result,
#                             contr_pars=c('b_G','b_E','b_S','b_Xpu','b_Xpr','b_Di'))
# save( contrast_CE, file=file.path(file_save, 'SOLV_CE_contrast.RData' ) )
# 
# # non-centered
# load( file.path(file_save, 'SOLV_NC_post.RData' ) )
# contrast_NC = contrast_stat(results_object=stan_result,
#                             post=result,
#                             contr_pars=c('b_G','b_E','b_S','b_Xpu','b_Xpr','b_Di'))
# save( contrast_NC, file=file.path(file_save, 'SOLV_NC_contrast.RData' ) )
# 
# # merge
# contrast_res = merge(contrast_CE, contrast_NC, by='parameter', all.x=T)
# save( contrast_res, file=file.path(file_save, 'SOLV_contrasts.RData' ) )

load(file.path(file_save, 'SOLV_NC_contrast.RData' ))
# str(contrast_res)


file_save = file.path(getwd(), 'figures')

# plot
png(file.path(file_save, 'SOLV_recovery_contrast.png'), 
    units='cm', width=30, height=15, res=100)

par(mfrow=c(1, 1), mar=c(9,4.1,4.1,2.1))

y_lim = range( with(contrast_res, c(X5.5..x, X94.5..x, X5.5..y, X94.5..y)), na.rm=T )

plot( (1:nrow(contrast_res))-0.1, contrast_res$mean.x, xlim=c(0, nrow(contrast_res)+1), 
      ylim=y_lim, xaxt='n', yaxt='n', col=col.alpha('black', 0.3), pch=19, 
      xlab='', ylab='estimates', main='Contrasts')
abline(h=0, col=col.alpha('black', 0.3), lty=2)
axis(side=1, at=1:nrow(contrast_res), labels=contrast_res$parameter, las=2)
axis(side=2, at=round( seq(y_lim[1], y_lim[2], by=0.1), 2), las=1)
for(i in 1:nrow(contrast_res)){
  lines(x=rep(i,2)-0.1, y=with( contrast_res[i,], c(X5.5..x, X94.5..x) ) ,
        col=col.alpha('black', 0.3))
}
points( (1:nrow(contrast_res))+0.1, contrast_res$mean.y, col=col.alpha('blue', 0.3), pch=19)
for(i in 1:nrow(contrast_res)){
  lines(x=rep(i,2)+0.1, y=with( contrast_res[i,], c(X5.5..y, X94.5..y) ) ,
        col=col.alpha('blue', 0.3))
}
legend('topleft', legend=c('centered','non-centered'), pch=c(19,19), cex=.8, 
       bty='n', col=c(col.alpha('black', 0.3), col.alpha('blue', 0.3)) )

par(mfrow=c(1, 1), mar=opar$mar)

dev.off()






# 6. model fit ####


## 6.1 centered ####

# load
file_save = file.path(getwd(), 'tables')
load( file.path(file_save, 'SOLV_CE_WAIC.RData' ) )
WAIC_CE_res = result

load( file.path(file_save, 'SOLV_CE_PSIS_OneChain.RData' ) )
PSIS_CE_res = result


# fit
apply(WAIC_CE_res[,-4], 2, sum)
apply(PSIS_CE_res[,-c(4:5)], 2, sum)


# outlier detection
file_save = file.path(getwd(), 'figures')

png(file.path(file_save, 'SOLV_CE_outlier.png'), 
    units='cm', width=30, height=15, res=100)

plot( PSIS_CE_res$k , WAIC_CE_res$penalty, 
      col=col.alpha('blue', 0.1), pch=19, lwd=2 , 
      xlab="PSIS Pareto k", ylab="WAIC penalty"  )
abline(v=0.5, lty=2)

dev.off()




## 6.2 non-centered ####

# load
file_save = file.path(getwd(), 'tables')
load( file.path(file_save, 'SOLV_NC_WAIC.RData' ) )
WAIC_NC_res = result

load( file.path(file_save, 'SOLV_NC_PSIS_OneChain.RData' ) )
PSIS_NC_res = result


# fit
apply(WAIC_NC_res[,-4], 2, sum)
apply(PSIS_NC_res[,-c(4:5)], 2, sum)

# outlier detection
file_save = file.path(getwd(), 'figures')

png(file.path(file_save, 'SOLV_NC_outlier.png'), 
    units='cm', width=30, height=15, res=100)

plot( PSIS_NC_res$k , WAIC_NC_res$penalty, 
      col=col.alpha('blue', 0.1), pch=19, lwd=2 , 
      xlab="PSIS Pareto k", ylab="WAIC penalty"  )
abline(v=0.5, lty=2)

dev.off()







# 7. posterior predictive ####

## 7.1 centered ####

### 7.1.1 ICC's and IIC's ####

# load
file_save = file.path(getwd(), 'tables')
load( file.path(file_save, 'SOLV_CE_post.RData' ) )


file_save = file.path(getwd(), 'figures')

# ICC
png(file.path(file_save, 'SOLV_CE_ICC_posterior.png'), 
    units='cm', width=55, height=80, res=100)

par(mfrow=c(7,4))
for(k in 1:25){
  
  # simulated, marginal, average, and true ICC
  ds_ICC = sim_ICC(sim_b=result$b_k, item_id=k)
  ds_ICC_mar = mar_ICC(sim_ICC_object=ds_ICC)
  ds_ICC_ave = ave_ICC(sim_b=result$b_k, item_id=k)
  
  # plot
  plot(NULL, xlim=c(-3,3), ylim=c(0,1), main=paste0('item ', k),
       xlab='individual ability', ylab='probability of correct item')
  abline(h=0.5, v=0, lty=2, col=col.alpha('black', 0.5))
  for(i in 1:500){ 
    lines(ds_ICC$theta, ds_ICC[,i+1], col=col.alpha('black',0.03)) 
  }
  lines( ds_ICC_mar, col='blue', lwd=3, lty=2 )
  lines( ds_ICC_ave, col='red', lwd=2, lty=3 ) 
  
}
par(mfrow=c(1,1))

dev.off()


# IIF
png(file.path(file_save, 'SOLV_CE_IIF_posterior.png'), 
    units='cm', width=55, height=80, res=100)

par(mfrow=c(7,4))
for(k in 1:25){
  
  # simulated, marginal, average, and true IIF
  ds_IIF = sim_IIF(sim_b=result$b_k, item_id=k)
  ds_IIF_mar = mar_IIF(sim_IIF_object=ds_IIF)
  ds_IIF_ave = ave_IIF(sim_b=result$b_k, item_id=k)
  
  # plot
  plot(NULL, xlim=c(-3,3), ylim=c(0, 0.3), main=paste0('item ', k),
       xlab='individual ability', ylab='probability of correct item')
  abline(h=0.5, v=0, lty=2, col=col.alpha('black', 0.5))
  for(i in 1:500){ 
    lines(ds_IIF$theta, ds_IIF[,i+1], col=col.alpha('black',0.03)) 
  }
  lines( ds_IIF_mar, col='blue', lwd=2, lty=2 )
  lines( ds_IIF_ave, col='red', lwd=2, lty=2 ) 
  
}
par(mfrow=c(1,1))

dev.off()





### 7.1.2 hit rates #####

# # load requirements
# file_save = file.path(getwd(), 'tables')
# load( file.path(file_save, 'SOLV_CE_post.RData' ) )
# # str(result)
# 
# file_load = file.path(getwd(), 'data', '4_final')
# load( file.path(file_load, 'data_application.RData' ) )
# # str(data_app)
# 
# 
# # calculations
# ds_prob_post = sim_hit_rate(d_true=data_app, d_sim=result, prob=T)
# save(ds_prob_post, file=file.path(file_save, 'SOLV_CE_prob_post.RData'))
# 
# ds_out_post = sim_hit_rate(d_true=data_app, d_sim=result, prob=F)
# save(ds_out_post, file=file.path(file_save, 'SOLV_CE_out_post.RData'))
# # str(post_sim)


# load sim rates
file_save = file.path(getwd(), 'tables')
load( file.path(file_save, 'SOLV_CE_prob_post.RData') )
load( file.path(file_save, 'SOLV_CE_out_post.RData') )


file_save = file.path(getwd(), 'figures')


#### a. individuals ####

# average and marginal hit rate
ds_ave = ave_hit_rate(sim_hit_object=ds_prob_post$IDind, prob=0.95)
ds_mar = ave_hit_rate(sim_hit_object=ds_out_post$IDind, prob=0.95)

table(ds_prob_post$IDind$ED)

# random sample
n = 50
set.seed(2546) # to make all comparable
idx = sample(ds_ave$IDj, size=n)
idx = idx[order(idx)]

ds_s = ds_prob_post$IDind[ ds_prob_post$IDind$IDj %in% idx,]
ds_ave_s = ds_ave[ds_ave$IDj %in% idx,]
ds_mar_s = ds_mar[ds_mar$IDj %in% idx,]


# plot
png(file.path(file_save, 'SOLV_CE_HitRate_ind.png'), 
    units='cm', width=30, height=25, res=100)

par(mfrow=c(2,1))

# shrinkage
plot(1:length(idx), ds_s$true, ylim=c(0,1), xaxt='n', col='black', 
     main='Shrinkage',
     xlab='Individual ID', ylab='Proportion of correct items')
axis(side=1, at=1:length(idx), labels=idx)
with(ds_ave_s, points((1:length(idx))-0.15, p, col=col.alpha('blue', 0.3), pch=19 ) )
for(j in 1:length(idx)){
  lines(x=rep( j-0.15, 2 ), y=c( ds_s$true[j], ds_ave_s$p[j] ), col=col.alpha('blue', 0.5))  
}
with(ds_mar_s, points((1:length(idx))+0.15, p, col=col.alpha('red', 0.3), pch=19 ) )
legend('topleft', legend=c('observed','average','marginal'), pch=c(1,19,19), 
       cex=.8, bty='n', col=c('black', col.alpha('blue',0.3), col.alpha('red',0.3) ) )


# average compatibility intervals
plot(1:length(idx), ds_s$true, ylim=c(0,1), xaxt='n', col='black', 
     main='Compatibility Intervals',
     xlab='Individual ID', ylab='Proportion of correct items')
axis(side=1, at=1:length(idx), labels=idx)
with(ds_ave_s, points((1:length(idx))-0.15, p, col=col.alpha('blue', 0.3), pch=19 ) )
for(j in 1:length(idx)){
  with(ds_ave_s[j,], lines(x=rep(j-0.15, 2), y=c(lower, upper), 
                           col=col.alpha('blue', 0.3)) )
}
with(ds_mar_s, points((1:length(idx))+0.15, p, col=col.alpha('red', 0.3), pch=19 ) )
for(j in 1:length(idx)){
  with(ds_mar_s[j,], lines(x=rep(j+0.15, 2), y=c(lower, upper), 
                           col=col.alpha('red', 0.3)) )
}

par(mfrow=c(1,1))

dev.off()



#### b. individuals per variable ####

png(file.path(file_save, 'SOLV_CE_HitRate_var.png'), 
    units='cm', width=30, height=45, res=100)

par(mfrow=c(4,2))

# gender
plot(ds_s$GE-0.1, ds_s$true, xlim=with(ds_s, c(min(GE)-1, max(GE)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.5), pch=1,
     xlab='', ylab='% correct items', main='gender')
axis(side=1, at=with(ds_s, min(GE):max(GE)), labels=c('male','female') ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(j in (1:200)+1 ){
  points(ds_s[,c('GE', paste0('X',j) )], col=col.alpha('blue', 0.002), pch=19 )
}
legend('topleft', legend=c('observed','predicted'), pch=c(1,19,19), 
       cex=.8, bty='n', col=c('black', col.alpha('blue',0.3), col.alpha('red',0.3) ) )


# age
plot(ds_s$AG-0.2, ds_s$true, xlim=with(ds_s, c(min(AG)-1, max(AG)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.5), pch=1,
     xlab='', ylab='% correct items', main='age')
axis(side=1, at=with(ds_s, min(AG):max(AG)) ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(j in (1:500)+1 ){
  points(ds_s[,c('AG', paste0('X',j) )], col=col.alpha('blue', 0.002), pch=19 )
}

# edu
plot(ds_s$ED-0.1, ds_s$true, xlim=with(ds_s, c(min(ED)-1, max(ED)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.5), pch=1,
     xlab='', ylab='% correct items', main='education')
axis(side=1, at=with(ds_s, min(ED):max(ED)), labels=c('inst','uni','both') ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(j in (1:200)+1 ){
  points(ds_s[,c('ED', paste0('X',j) )], col=col.alpha('blue', 0.002), pch=19 )
}

# specialty
plot(ds_s$SP-0.1, ds_s$true, xlim=with(ds_s, c(min(SP)-1, max(SP)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.5), pch=1,
     xlab='', ylab='% correct items', main='Specialty')
axis(side=1, at=with(ds_s, min(SP):max(SP)), 
     labels=c('early child.','primary','secondary') ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(j in (1:200)+1 ){
  points(ds_s[,c('SP', paste0('X',j) )], col=col.alpha('blue', 0.002), pch=19 )
}

# public experience
plot(ds_s$XPpu-0.1, ds_s$true, xlim=with(ds_s, c(min(XPpu)-1, max(XPpu)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.5), pch=1,
     xlab='', ylab='% correct items', main='Public experience')
axis(side=1, at=with(ds_s, min(XPpu):max(XPpu)), 
     labels=c('0y','<5y','6y-10y','10+y') ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(j in (1:200)+1 ){
  points(ds_s[,c('XPpu', paste0('X',j) )], col=col.alpha('blue', 0.002), pch=19 )
}

# private experience
plot(ds_s$XPpr-0.1, ds_s$true, xlim=with(ds_s, c(min(XPpr)-1, max(XPpr)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.5), pch=1,
     xlab='', ylab='% correct items', main='Public experience')
axis(side=1, at=with(ds_s, min(XPpr):max(XPpr)), 
     labels=c('0y','<5y','6y-10y','10+y') ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(j in (1:200)+1 ){
  points(ds_s[,c('XPpr', paste0('X',j) )], col=col.alpha('blue', 0.002), pch=19 )
}

# disability
plot(ds_s$DI-0.1, ds_s$true, xlim=with(ds_s, c(min(DI)-1, max(DI)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.5), pch=1,
     xlab='', ylab='% correct items', main='Public experience')
axis(side=1, at=with(ds_s, min(DI):max(DI)), 
     labels=c('none','low vision','motor') ) #,'auditory' 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(j in (1:200)+1 ){
  points(ds_s[,c('DI', paste0('X',j) )], col=col.alpha('blue', 0.002), pch=19 )
}

par(mfrow=c(1,1))

dev.off()




#### c. items ####

# average and marginal hit rate
ds_s = ds_prob_post$IDitems
ds_ave_s = ave_hit_rate(sim_hit_object=ds_prob_post$IDitems, prob=0.95)
ds_mar_s = ave_hit_rate(sim_hit_object=ds_out_post$IDitems, prob=0.95)


# plot
png(file.path(file_save, 'SOLV_CE_HitRate_items.png'), 
    units='cm', width=30, height=25, res=100)

par(mfrow=c(2,1))

# shrinkage
plot(1:nrow(ds_s), ds_s$true, xlim=c(0, nrow(ds_s)+1), ylim=c(0,1), 
     xaxt='n', col='black', main='Shrinkage',
     xlab='Item ID', ylab='Proportion of correct individuals')
axis(side=1, at=1:nrow(ds_s), labels=ds_s$IDk)
with(ds_ave_s, points((1:nrow(ds_s))-0.15, p, col=col.alpha('blue', 0.3), pch=19 ) )
for(j in 1:nrow(ds_s)){
  lines(x=rep( j-0.15, 2 ), y=c( ds_s$true[j], ds_ave_s$p[j] ), col=col.alpha('blue', 0.5))  
}
with(ds_mar_s, points((1:nrow(ds_s))+0.15, p, col=col.alpha('red', 0.3), pch=19 ) )
legend('topleft', legend=c('observed','average','marginal'), pch=c(1,19,19), 
       cex=.8, bty='n', col=c('black', col.alpha('blue',0.3), col.alpha('red',0.3) ) )


# average compatibility intervals
plot(1:nrow(ds_s), ds_s$true, xlim=c(0, nrow(ds_s)+1), ylim=c(0,1), 
     xaxt='n', col='black', main='Compatibility Intervals',
     xlab='Item ID', ylab='Proportion of correct individuals')
axis(side=1, at=1:nrow(ds_s), labels=ds_s$IDk)
with(ds_ave_s, points((1:nrow(ds_s))-0.15, p, col=col.alpha('blue', 0.3), pch=19 ) )
for(j in 1:nrow(ds_s)){
  with(ds_ave_s[j,], lines(x=rep(j-0.15, 2), y=c(lower, upper), 
                           col=col.alpha('blue', 0.3)) )
}
with(ds_mar_s, points((1:nrow(ds_s))+0.15, p, col=col.alpha('red', 0.3), pch=19 ) )
for(j in 1:nrow(ds_s)){
  with(ds_mar_s[j,], lines(x=rep(j+0.15, 2), y=c(lower, upper), 
                           col=col.alpha('red', 0.3)) )
}

par(mfrow=c(1,1))

dev.off()




#### d. texts ####

# average and marginal hit rate
ds_s = ds_prob_post$IDtext
ds_ave_s = ave_hit_rate(sim_hit_object=ds_prob_post$IDtext, prob=0.95)
ds_mar_s = ave_hit_rate(sim_hit_object=ds_out_post$IDtext, prob=0.95)


# plot
png(file.path(file_save, 'SOLV_CE_HitRate_text.png'), 
    units='cm', width=30, height=25, res=100)

par(mfrow=c(2,1))

# shrinkage
plot(1:nrow(ds_s), ds_s$true, xlim=c(0, nrow(ds_s)+1), ylim=c(0,1), 
     xaxt='n', col='black', main='Shrinkage',
     xlab='Text ID', ylab='Proportion of correct individuals')
axis(side=1, at=1:nrow(ds_s), labels=ds_s$IDl)
with(ds_ave_s, points((1:nrow(ds_s))-0.15, p, col=col.alpha('blue', 0.3), pch=19 ) )
for(j in 1:nrow(ds_s)){
  lines(x=rep( j-0.15, 2 ), y=c( ds_s$true[j], ds_ave_s$p[j] ), col=col.alpha('blue', 0.5))  
}
with(ds_mar_s, points((1:nrow(ds_s))+0.15, p, col=col.alpha('red', 0.3), pch=19 ) )
legend('topleft', legend=c('observed','average','marginal'), pch=c(1,19,19), 
       cex=.8, bty='n', col=c('black', col.alpha('blue',0.3), col.alpha('red',0.3) ) )


# average compatibility intervals
plot(1:nrow(ds_s), ds_s$true, xlim=c(0, nrow(ds_s)+1), ylim=c(0,1), 
     xaxt='n', col='black', main='Compatibility Intervals',
     xlab='Text ID', ylab='Proportion of correct individuals')
axis(side=1, at=1:nrow(ds_s), labels=ds_s$IDl)
with(ds_ave_s, points((1:nrow(ds_s))-0.15, p, col=col.alpha('blue', 0.3), pch=19 ) )
for(j in 1:nrow(ds_s)){
  with(ds_ave_s[j,], lines(x=rep(j-0.15, 2), y=c(lower, upper), 
                           col=col.alpha('blue', 0.3)) )
}
with(ds_mar_s, points((1:nrow(ds_s))+0.15, p, col=col.alpha('red', 0.3), pch=19 ) )
for(j in 1:nrow(ds_s)){
  with(ds_mar_s[j,], lines(x=rep(j+0.15, 2), y=c(lower, upper), 
                           col=col.alpha('red', 0.3)) )
}

par(mfrow=c(1,1))

dev.off()




#### e. dimensions ####

# average and marginal hit rate
ds_s = ds_prob_post$IDdim
ds_ave_s = ave_hit_rate(sim_hit_object=ds_prob_post$IDdim, prob=0.95)
ds_mar_s = ave_hit_rate(sim_hit_object=ds_out_post$IDdim, prob=0.95)


# plot
png(file.path(file_save, 'SOLV_CE_HitRate_dim.png'), 
    units='cm', width=30, height=25, res=100)

par(mfrow=c(2,1))

# shrinkage
plot(1:nrow(ds_s), ds_s$true, xlim=c(0, nrow(ds_s)+1), ylim=c(0,1), 
     xaxt='n', col='black', main='Shrinkage',
     xlab='Dimension ID', ylab='Proportion of correct individuals')
axis(side=1, at=1:nrow(ds_s), labels=ds_s$IDd)
with(ds_ave_s, points((1:nrow(ds_s))-0.15, p, col=col.alpha('blue', 0.3), pch=19 ) )
for(j in 1:nrow(ds_s)){
  lines(x=rep( j-0.15, 2 ), y=c( ds_s$true[j], ds_ave_s$p[j] ), col=col.alpha('blue', 0.5))  
}
with(ds_mar_s, points((1:nrow(ds_s))+0.15, p, col=col.alpha('red', 0.3), pch=19 ) )
legend('topleft', legend=c('observed','average','marginal'), pch=c(1,19,19), 
       cex=.8, bty='n', col=c('black', col.alpha('blue',0.3), col.alpha('red',0.3) ) )


# average compatibility intervals
plot(1:nrow(ds_s), ds_s$true, xlim=c(0, nrow(ds_s)+1), ylim=c(0,1), 
     xaxt='n', col='black', main='Compatibility Intervals',
     xlab='Dimension ID', ylab='Proportion of correct individuals')
axis(side=1, at=1:nrow(ds_s), labels=ds_s$IDd)
with(ds_ave_s, points((1:nrow(ds_s))-0.15, p, col=col.alpha('blue', 0.3), pch=19 ) )
for(j in 1:nrow(ds_s)){
  with(ds_ave_s[j,], lines(x=rep(j-0.15, 2), y=c(lower, upper), 
                           col=col.alpha('blue', 0.3)) )
}
with(ds_mar_s, points((1:nrow(ds_s))+0.15, p, col=col.alpha('red', 0.3), pch=19 ) )
for(j in 1:nrow(ds_s)){
  with(ds_mar_s[j,], lines(x=rep(j+0.15, 2), y=c(lower, upper), 
                           col=col.alpha('red', 0.3)) )
}

par(mfrow=c(1,1))

dev.off()





## 7.2 non-centered ####

### 7.2.1 ICC's and IIC's ####

# load
file_save = file.path(getwd(), 'tables')
load( file.path(file_save, 'SOLV_NC_post.RData' ) )


file_save = file.path(getwd(), 'figures')

# ICC
png(file.path(file_save, 'SOLV_NC_ICC_posterior.png'), 
    units='cm', width=55, height=80, res=100)

par(mfrow=c(7,4))
for(k in 1:25){
  
  # simulated, marginal, average, and true ICC
  ds_ICC = sim_ICC(sim_b=result$b_k, item_id=k)
  ds_ICC_mar = mar_ICC(sim_ICC_object=ds_ICC)
  ds_ICC_ave = ave_ICC(sim_b=result$b_k, item_id=k)
  
  # plot
  plot(NULL, xlim=c(-3,3), ylim=c(0,1), main=paste0('item ', k),
       xlab='individual ability', ylab='probability of correct item')
  abline(h=0.5, v=0, lty=2, col=col.alpha('black', 0.5))
  for(i in 1:500){ 
    lines(ds_ICC$theta, ds_ICC[,i+1], col=col.alpha('black',0.03)) 
  }
  lines( ds_ICC_mar, col='blue', lwd=3, lty=2 )
  lines( ds_ICC_ave, col='red', lwd=2, lty=3 ) 
  
}
par(mfrow=c(1,1))

dev.off()


# IIF
png(file.path(file_save, 'SOLV_NC_IIF_posterior.png'), 
    units='cm', width=55, height=80, res=100)

par(mfrow=c(7,4))
for(k in 1:25){
  
  # simulated, marginal, average, and true IIF
  ds_IIF = sim_IIF(sim_b=result$b_k, item_id=k)
  ds_IIF_mar = mar_IIF(sim_IIF_object=ds_IIF)
  ds_IIF_ave = ave_IIF(sim_b=result$b_k, item_id=k)
  
  # plot
  plot(NULL, xlim=c(-3,3), ylim=c(0, 0.3), main=paste0('item ', k),
       xlab='individual ability', ylab='probability of correct item')
  abline(h=0.5, v=0, lty=2, col=col.alpha('black', 0.5))
  for(i in 1:500){ 
    lines(ds_IIF$theta, ds_IIF[,i+1], col=col.alpha('black',0.03)) 
  }
  lines( ds_IIF_mar, col='blue', lwd=2, lty=2 )
  lines( ds_IIF_ave, col='red', lwd=2, lty=2 ) 
  
}
par(mfrow=c(1,1))

dev.off()





### 7.2.2 hit rates #####

# # load requirements
# file_save = file.path(getwd(), 'tables')
# load( file.path(file_save, 'SOLV_NC_post.RData' ) )
# # str(result)
# 
# file_load = file.path(getwd(), 'data', '4_final')
# load( file.path(file_load, 'data_application.RData' ) )
# # str(data_app)
# 
# 
# # calculations
# ds_prob_post = sim_hit_rate(d_true=data_app, d_sim=result, prob=T)
# save(ds_prob_post, file=file.path(file_save, 'SOLV_NC_prob_post.RData'))
# 
# ds_out_post = sim_hit_rate(d_true=data_app, d_sim=result, prob=F)
# save(ds_out_post, file=file.path(file_save, 'SOLV_NC_out_post.RData'))
# # str(post_sim)


# load sim rates
file_save = file.path(getwd(), 'tables')
load( file.path(file_save, 'SOLV_NC_prob_post.RData') )
load( file.path(file_save, 'SOLV_NC_out_post.RData') )


file_save = file.path(getwd(), 'figures')


#### a. individuals ####

# average and marginal hit rate
ds_ave = ave_hit_rate(sim_hit_object=ds_prob_post$IDind, prob=0.95)
ds_mar = ave_hit_rate(sim_hit_object=ds_out_post$IDind, prob=0.95)

table(ds_prob_post$IDind$ED)

# random sample
n = 50
set.seed(2546) # to make all comparable
idx = sample(ds_ave$IDj, size=n)
idx = idx[order(idx)]

ds_s = ds_prob_post$IDind[ ds_prob_post$IDind$IDj %in% idx,]
ds_ave_s = ds_ave[ds_ave$IDj %in% idx,]
ds_mar_s = ds_mar[ds_mar$IDj %in% idx,]


# plot
png(file.path(file_save, 'SOLV_NC_HitRate_ind.png'), 
    units='cm', width=30, height=25, res=100)

par(mfrow=c(2,1))

# shrinkage
plot(1:length(idx), ds_s$true, ylim=c(0,1), xaxt='n', col='black', 
     main='Shrinkage',
     xlab='Individual ID', ylab='Proportion of correct items')
axis(side=1, at=1:length(idx), labels=idx)
with(ds_ave_s, points((1:length(idx))-0.15, p, col=col.alpha('blue', 0.3), pch=19 ) )
for(j in 1:length(idx)){
  lines(x=rep( j-0.15, 2 ), y=c( ds_s$true[j], ds_ave_s$p[j] ), col=col.alpha('blue', 0.5))  
}
with(ds_mar_s, points((1:length(idx))+0.15, p, col=col.alpha('red', 0.3), pch=19 ) )
legend('topleft', legend=c('observed','average','marginal'), pch=c(1,19,19), 
       cex=.8, bty='n', col=c('black', col.alpha('blue',0.3), col.alpha('red',0.3) ) )


# average compatibility intervals
plot(1:length(idx), ds_s$true, ylim=c(0,1), xaxt='n', col='black', 
     main='Compatibility Intervals',
     xlab='Individual ID', ylab='Proportion of correct items')
axis(side=1, at=1:length(idx), labels=idx)
with(ds_ave_s, points((1:length(idx))-0.15, p, col=col.alpha('blue', 0.3), pch=19 ) )
for(j in 1:length(idx)){
  with(ds_ave_s[j,], lines(x=rep(j-0.15, 2), y=c(lower, upper), 
                           col=col.alpha('blue', 0.3)) )
}
with(ds_mar_s, points((1:length(idx))+0.15, p, col=col.alpha('red', 0.3), pch=19 ) )
for(j in 1:length(idx)){
  with(ds_mar_s[j,], lines(x=rep(j+0.15, 2), y=c(lower, upper), 
                           col=col.alpha('red', 0.3)) )
}

par(mfrow=c(1,1))

dev.off()



#### b. individuals per variable ####

png(file.path(file_save, 'SOLV_NC_HitRate_var.png'), 
    units='cm', width=30, height=45, res=100)

par(mfrow=c(4,2))

# gender
plot(ds_s$GE-0.1, ds_s$true, xlim=with(ds_s, c(min(GE)-1, max(GE)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.5), pch=1,
     xlab='', ylab='% correct items', main='gender')
axis(side=1, at=with(ds_s, min(GE):max(GE)), labels=c('male','female') ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(j in (1:200)+1 ){
  points(ds_s[,c('GE', paste0('X',j) )], col=col.alpha('blue', 0.002), pch=19 )
}
legend('topleft', legend=c('observed','predicted'), pch=c(1,19,19), 
       cex=.8, bty='n', col=c('black', col.alpha('blue',0.3), col.alpha('red',0.3) ) )


# age
plot(ds_s$AG-0.2, ds_s$true, xlim=with(ds_s, c(min(AG)-1, max(AG)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.5), pch=1,
     xlab='', ylab='% correct items', main='age')
axis(side=1, at=with(ds_s, min(AG):max(AG)) ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(j in (1:500)+1 ){
  points(ds_s[,c('AG', paste0('X',j) )], col=col.alpha('blue', 0.002), pch=19 )
}

# edu
plot(ds_s$ED-0.1, ds_s$true, xlim=with(ds_s, c(min(ED)-1, max(ED)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.5), pch=1,
     xlab='', ylab='% correct items', main='education')
axis(side=1, at=with(ds_s, min(ED):max(ED)), labels=c('inst','uni','both') ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(j in (1:200)+1 ){
  points(ds_s[,c('ED', paste0('X',j) )], col=col.alpha('blue', 0.002), pch=19 )
}

# specialty
plot(ds_s$SP-0.1, ds_s$true, xlim=with(ds_s, c(min(SP)-1, max(SP)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.5), pch=1,
     xlab='', ylab='% correct items', main='Specialty')
axis(side=1, at=with(ds_s, min(SP):max(SP)), 
     labels=c('early child.','primary','secondary') ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(j in (1:200)+1 ){
  points(ds_s[,c('SP', paste0('X',j) )], col=col.alpha('blue', 0.002), pch=19 )
}

# public experience
plot(ds_s$XPpu-0.1, ds_s$true, xlim=with(ds_s, c(min(XPpu)-1, max(XPpu)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.5), pch=1,
     xlab='', ylab='% correct items', main='Public experience')
axis(side=1, at=with(ds_s, min(XPpu):max(XPpu)), 
     labels=c('0y','<5y','6y-10y','10+y') ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(j in (1:200)+1 ){
  points(ds_s[,c('XPpu', paste0('X',j) )], col=col.alpha('blue', 0.002), pch=19 )
}

# private experience
plot(ds_s$XPpr-0.1, ds_s$true, xlim=with(ds_s, c(min(XPpr)-1, max(XPpr)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.5), pch=1,
     xlab='', ylab='% correct items', main='Public experience')
axis(side=1, at=with(ds_s, min(XPpr):max(XPpr)), 
     labels=c('0y','<5y','6y-10y','10+y') ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(j in (1:200)+1 ){
  points(ds_s[,c('XPpr', paste0('X',j) )], col=col.alpha('blue', 0.002), pch=19 )
}

# disability
plot(ds_s$DI-0.1, ds_s$true, xlim=with(ds_s, c(min(DI)-1, max(DI)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.5), pch=1,
     xlab='', ylab='% correct items', main='Public experience')
axis(side=1, at=with(ds_s, min(DI):max(DI)), 
     labels=c('none','low vision','motor') ) #,'auditory' 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(j in (1:200)+1 ){
  points(ds_s[,c('DI', paste0('X',j) )], col=col.alpha('blue', 0.002), pch=19 )
}

par(mfrow=c(1,1))

dev.off()




#### c. items ####

# average and marginal hit rate
ds_s = ds_prob_post$IDitems
ds_ave_s = ave_hit_rate(sim_hit_object=ds_prob_post$IDitems, prob=0.95)
ds_mar_s = ave_hit_rate(sim_hit_object=ds_out_post$IDitems, prob=0.95)


# plot
png(file.path(file_save, 'SOLV_NC_HitRate_items.png'), 
    units='cm', width=30, height=25, res=100)

par(mfrow=c(2,1))

# shrinkage
plot(1:nrow(ds_s), ds_s$true, xlim=c(0, nrow(ds_s)+1), ylim=c(0,1), 
     xaxt='n', col='black', main='Shrinkage',
     xlab='Item ID', ylab='Proportion of correct individuals')
axis(side=1, at=1:nrow(ds_s), labels=ds_s$IDk)
with(ds_ave_s, points((1:nrow(ds_s))-0.15, p, col=col.alpha('blue', 0.3), pch=19 ) )
for(j in 1:nrow(ds_s)){
  lines(x=rep( j-0.15, 2 ), y=c( ds_s$true[j], ds_ave_s$p[j] ), col=col.alpha('blue', 0.5))  
}
with(ds_mar_s, points((1:nrow(ds_s))+0.15, p, col=col.alpha('red', 0.3), pch=19 ) )
legend('topleft', legend=c('observed','average','marginal'), pch=c(1,19,19), 
       cex=.8, bty='n', col=c('black', col.alpha('blue',0.3), col.alpha('red',0.3) ) )


# average compatibility intervals
plot(1:nrow(ds_s), ds_s$true, xlim=c(0, nrow(ds_s)+1), ylim=c(0,1), 
     xaxt='n', col='black', main='Compatibility Intervals',
     xlab='Item ID', ylab='Proportion of correct individuals')
axis(side=1, at=1:nrow(ds_s), labels=ds_s$IDk)
with(ds_ave_s, points((1:nrow(ds_s))-0.15, p, col=col.alpha('blue', 0.3), pch=19 ) )
for(j in 1:nrow(ds_s)){
  with(ds_ave_s[j,], lines(x=rep(j-0.15, 2), y=c(lower, upper), 
                           col=col.alpha('blue', 0.3)) )
}
with(ds_mar_s, points((1:nrow(ds_s))+0.15, p, col=col.alpha('red', 0.3), pch=19 ) )
for(j in 1:nrow(ds_s)){
  with(ds_mar_s[j,], lines(x=rep(j+0.15, 2), y=c(lower, upper), 
                           col=col.alpha('red', 0.3)) )
}

par(mfrow=c(1,1))

dev.off()




#### d. texts ####

# average and marginal hit rate
ds_s = ds_prob_post$IDtext
ds_ave_s = ave_hit_rate(sim_hit_object=ds_prob_post$IDtext, prob=0.95)
ds_mar_s = ave_hit_rate(sim_hit_object=ds_out_post$IDtext, prob=0.95)


# plot
png(file.path(file_save, 'SOLV_NC_HitRate_text.png'), 
    units='cm', width=30, height=25, res=100)

par(mfrow=c(2,1))

# shrinkage
plot(1:nrow(ds_s), ds_s$true, xlim=c(0, nrow(ds_s)+1), ylim=c(0,1), 
     xaxt='n', col='black', main='Shrinkage',
     xlab='Text ID', ylab='Proportion of correct individuals')
axis(side=1, at=1:nrow(ds_s), labels=ds_s$IDl)
with(ds_ave_s, points((1:nrow(ds_s))-0.15, p, col=col.alpha('blue', 0.3), pch=19 ) )
for(j in 1:nrow(ds_s)){
  lines(x=rep( j-0.15, 2 ), y=c( ds_s$true[j], ds_ave_s$p[j] ), col=col.alpha('blue', 0.5))  
}
with(ds_mar_s, points((1:nrow(ds_s))+0.15, p, col=col.alpha('red', 0.3), pch=19 ) )
legend('topleft', legend=c('observed','average','marginal'), pch=c(1,19,19), 
       cex=.8, bty='n', col=c('black', col.alpha('blue',0.3), col.alpha('red',0.3) ) )


# average compatibility intervals
plot(1:nrow(ds_s), ds_s$true, xlim=c(0, nrow(ds_s)+1), ylim=c(0,1), 
     xaxt='n', col='black', main='Compatibility Intervals',
     xlab='Text ID', ylab='Proportion of correct individuals')
axis(side=1, at=1:nrow(ds_s), labels=ds_s$IDl)
with(ds_ave_s, points((1:nrow(ds_s))-0.15, p, col=col.alpha('blue', 0.3), pch=19 ) )
for(j in 1:nrow(ds_s)){
  with(ds_ave_s[j,], lines(x=rep(j-0.15, 2), y=c(lower, upper), 
                           col=col.alpha('blue', 0.3)) )
}
with(ds_mar_s, points((1:nrow(ds_s))+0.15, p, col=col.alpha('red', 0.3), pch=19 ) )
for(j in 1:nrow(ds_s)){
  with(ds_mar_s[j,], lines(x=rep(j+0.15, 2), y=c(lower, upper), 
                           col=col.alpha('red', 0.3)) )
}

par(mfrow=c(1,1))

dev.off()




#### e. dimensions ####

# average and marginal hit rate
ds_s = ds_prob_post$IDdim
ds_ave_s = ave_hit_rate(sim_hit_object=ds_prob_post$IDdim, prob=0.95)
ds_mar_s = ave_hit_rate(sim_hit_object=ds_out_post$IDdim, prob=0.95)


# plot
png(file.path(file_save, 'SOLV_NC_HitRate_dim.png'), 
    units='cm', width=30, height=25, res=100)

par(mfrow=c(2,1))

# shrinkage
plot(1:nrow(ds_s), ds_s$true, xlim=c(0, nrow(ds_s)+1), ylim=c(0,1), 
     xaxt='n', col='black', main='Shrinkage',
     xlab='Dimension ID', ylab='Proportion of correct individuals')
axis(side=1, at=1:nrow(ds_s), labels=ds_s$IDd)
with(ds_ave_s, points((1:nrow(ds_s))-0.15, p, col=col.alpha('blue', 0.3), pch=19 ) )
for(j in 1:nrow(ds_s)){
  lines(x=rep( j-0.15, 2 ), y=c( ds_s$true[j], ds_ave_s$p[j] ), col=col.alpha('blue', 0.5))  
}
with(ds_mar_s, points((1:nrow(ds_s))+0.15, p, col=col.alpha('red', 0.3), pch=19 ) )
legend('topleft', legend=c('observed','average','marginal'), pch=c(1,19,19), 
       cex=.8, bty='n', col=c('black', col.alpha('blue',0.3), col.alpha('red',0.3) ) )


# average compatibility intervals
plot(1:nrow(ds_s), ds_s$true, xlim=c(0, nrow(ds_s)+1), ylim=c(0,1), 
     xaxt='n', col='black', main='Compatibility Intervals',
     xlab='Dimension ID', ylab='Proportion of correct individuals')
axis(side=1, at=1:nrow(ds_s), labels=ds_s$IDd)
with(ds_ave_s, points((1:nrow(ds_s))-0.15, p, col=col.alpha('blue', 0.3), pch=19 ) )
for(j in 1:nrow(ds_s)){
  with(ds_ave_s[j,], lines(x=rep(j-0.15, 2), y=c(lower, upper), 
                           col=col.alpha('blue', 0.3)) )
}
with(ds_mar_s, points((1:nrow(ds_s))+0.15, p, col=col.alpha('red', 0.3), pch=19 ) )
for(j in 1:nrow(ds_s)){
  with(ds_mar_s[j,], lines(x=rep(j+0.15, 2), y=c(lower, upper), 
                           col=col.alpha('red', 0.3)) )
}

par(mfrow=c(1,1))

dev.off()
