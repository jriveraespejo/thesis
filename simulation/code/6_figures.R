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
setwd('~/Desktop/simulations/#final')


# extras
source('3_functions_extra.R')
file_save = file.path(getwd(), 'figures')
opar = par()


# Chapter 3 ####

## 1. Devil's funnel ####

### 1.1 model ####
mcmc_code <- "
transformed data {
  int<lower=0> J;
  J = 1;
}
parameters {
  real theta[J];
  real v;
}
model {
  v ~ normal(0, 3);
  theta ~ normal(0, exp(v));
}
"
save_code = '1_stan_CE_simple.stan'
writeLines(mcmc_code, con=file.path(file_save, save_code) )

# set_cmdstan_path('~/cmdstan')
# mod = cmdstan_model( file.path(file_save, save_code) )
# mod$sample(output_dir=file.path(file_save), output_basename='1_stan_CE_simple',
#            chains=3, parallel_chains=3, seed=1 )


# load and transform simulations
fit_files = file.path(file_save, 
                      c('1_stan_CE_simple-1.csv',
                        '1_stan_CE_simple-2.csv',
                        '1_stan_CE_simple-3.csv') )
stan_CE = rstan::read_stan_csv(fit_files)
precis(stan_CE, depth=4)
# notice the small n-eff



### 1.2 figures ####

png(file.path(file_save, '1_trace_CE_simple.png'), units='cm', width=25, height=10, res=100)
idx = c('theta[1]', 'v')
par(mfrow=c(1,2))
traceplot_ulam(stan_CE, pars=idx, n_cols=2) 
par(mfrow=c(1,1))
dev.off()

png(file.path(file_save, '1_trank_CE_simple.png'), units='cm', width=25, height=10, res=100)
idx = c('theta[1]', 'v')
par(mfrow=c(1,3))
trankplot(stan_CE, pars=idx, n_cols=2) 
par(mfrow=c(1,1))
dev.off()

png(file.path(file_save, '1_acf_CE_simple.png'), units='cm', width=25, height=10, res=100)
idx = c('theta.1', 'v')
par(mfrow=c(1,2), mar=c(3,2.1,2.1,2.1))
acf_plot(stan_CE, pars=idx) 
par(mfrow=c(1,1), mar=opar$mar)
dev.off()



samples_CE = stan_CE@sim$samples[[1]]
example = data.frame(samples_CE[,1:2])
example$divergent = attr(samples_CE, 'sampler_params')$divergent__

sd=3
v = seq(-sd*qnorm(0.975), sd*qnorm(0.975), by=0.01)
z = seq(-2, 2, by=0.2)
theta = sapply(z, function(i){ exp(v)*i } )

png(file.path(file_save, '1_funnel_CE_simple.png'), units='cm', width=25, height=15, res=100)
plot(NULL, xlim=range(theta), ylim=range(v), xlab='theta', ylab='v')
for(i in 1:ncol(theta)){
  lines(x=theta[,i], y=v, col=col.alpha('black', 0.3), pch=19)
}
points(example[ example$divergent==0, 1:2], col=col.alpha('blue', 0.3), pch=16) 
points(example[ example$divergent==1, 1:2], col=col.alpha('red', 0.3), pch=16) 
mtext( paste0('# divergent = ', sum(example$divergent)), side=3, adj=0, cex=1)
dev.off()



  
## 2. Devil's funnel vs priors ####

### 2.1 model ####
mcmc_code <- "
transformed data {
  int<lower=0> J;
  J = 1;
}
parameters {
  real theta[J];
  real v;
}
model {
  v ~ normal(0, 1);
  theta ~ normal(0, exp(v));
}
"
save_code = '2_stan_CE_priors.stan'
writeLines(mcmc_code, con=file.path(file_save, save_code) )

# set_cmdstan_path('~/cmdstan')
# mod = cmdstan_model( file.path(file_save, save_code) )
# mod$sample(output_dir=file.path(file_save), output_basename='2_stan_CE_priors',
#            chains=3, parallel_chains=3, seed=1 )


# load and transform simulations
fit_files = file.path(file_save, 
                      c('2_stan_CE_priors-1.csv',
                        '2_stan_CE_priors-2.csv',
                        '2_stan_CE_priors-3.csv') )
stan_CE = rstan::read_stan_csv(fit_files)
precis(stan_CE, depth=4)
# notice the small n-eff



### 2.2 figures ####

png(file.path(file_save, '2_trace_CE_priors.png'), units='cm', width=25, height=10, res=100)
idx = c('theta[1]', 'v')
par(mfrow=c(1,2))
traceplot_ulam(stan_CE, pars=idx, n_cols=2) 
par(mfrow=c(1,1))
dev.off()

png(file.path(file_save, '2_trank_CE_priors.png'), units='cm', width=25, height=10, res=100)
idx = c('theta[1]', 'v')
par(mfrow=c(1,3))
trankplot(stan_CE, pars=idx, n_cols=2) 
par(mfrow=c(1,1))
dev.off()

png(file.path(file_save, '2_acf_CE_priors.png'), units='cm', width=25, height=10, res=100)
idx = c('theta.1', 'v')
par(mfrow=c(1,2), mar=c(3,2.1,2.1,2.1))
acf_plot(stan_CE, pars=idx) 
par(mfrow=c(1,1), mar=opar$mar)
dev.off()



samples_CE = stan_CE@sim$samples[[1]]
example = data.frame(samples_CE[,1:2])
example$divergent = attr(samples_CE, 'sampler_params')$divergent__

sd=1
v = seq(-sd*qnorm(0.975), sd*qnorm(0.975), by=0.01)
z = seq(-2, 2, by=0.2)
theta = sapply(z, function(i){ exp(v)*i } )

png(file.path(file_save, '2_funnel_CE_priors.png'), units='cm', width=25, height=15, res=100)
plot(NULL, xlim=range(theta), ylim=range(v), xlab='theta', ylab='v')
for(i in 1:ncol(theta)){
  lines(x=theta[,i], y=v, col=col.alpha('black', 0.3), pch=19)
}
points(example[ example$divergent==0, 1:2], col=col.alpha('blue', 0.3), pch=16) 
points(example[ example$divergent==1, 1:2], col=col.alpha('red', 0.3), pch=16)
mtext( paste0('# divergent = ', sum(example$divergent)), side=3, adj=0, cex=1)
dev.off()





## 3. Devil's funnel vs NC ####

### 3.1 model ####
mcmc_code <- "
transformed data {
  int<lower=0> J;
  J = 1;
}
parameters {
  real ztheta[J];
  real v;
}
transformed parameters{
  vector[J] theta;
  theta = exp(v) * to_vector( ztheta );
}
model {
  v ~ normal(0, 3);
  ztheta ~ normal(0, 1);
}
"
save_code = '3_stan_NC.stan'
writeLines(mcmc_code, con=file.path(file_save, save_code) )

# set_cmdstan_path('~/cmdstan')
# mod = cmdstan_model( file.path(file_save, save_code) )
# mod$sample(output_dir=file.path(file_save), output_basename='3_stan_NC',
#            chains=3, parallel_chains=3, seed=1 )


# load and transform simulations
fit_files = file.path(file_save, 
                      c('3_stan_NC-1.csv',
                        '3_stan_NC-2.csv',
                        '3_stan_NC-3.csv') )
stan_NC = rstan::read_stan_csv(fit_files)
precis(stan_NC, depth=4)
# notice the small n-eff



### 3.2 figures ####

png(file.path(file_save, '3_trace_NC.png'), units='cm', width=35, height=10, res=100)
idx = c('theta[1]','ztheta[1]', 'v')
par(mfrow=c(1,3))
traceplot_ulam(stan_NC, pars=idx, n_cols=3) 
par(mfrow=c(1,1))
dev.off()

png(file.path(file_save, '3_trank_NC.png'), units='cm', width=35, height=10, res=100)
idx = c('theta[1]','ztheta[1]', 'v')
par(mfrow=c(1,3))
trankplot(stan_NC, pars=idx, n_cols=3) 
par(mfrow=c(1,1))
dev.off()

png(file.path(file_save, '3_acf_NC.png'), units='cm', width=35, height=10, res=100)
idx = c('theta.1', 'ztheta.1', 'v')
par(mfrow=c(1,3), mar=c(3,2.1,2.1,2.1))
acf_plot(stan_NC, pars=idx) 
par(mfrow=c(1,1), mar=opar$mar)
dev.off()



samples_CE = stan_NC@sim$samples[[1]]
example = data.frame(samples_CE[,1:3])
example$divergent = attr(samples_CE, 'sampler_params')$divergent__

sd=3
v = seq(-sd*qnorm(0.975), sd*qnorm(0.975), by=0.01)
z = seq(-2, 2, by=0.2)
theta = sapply(z, function(i){ exp(v)*i } )

png(file.path(file_save, '3_funnel_NC.png'), units='cm', width=35, height=15, res=100)

par(mfrow=c(1, 2))

# actual sample space
plot(NULL, xlim=range(z), ylim=range(v), xlab='ztheta', ylab='v')
for ( l in c(0.1,0.3,0.5,0.8,0.99) ){
  lines( ellipse( diag( c(1, 3) ), centre=c(0,0), level=l), 
         col=col.alpha('black', 0.3))
}
points(example[ example$divergent==0, 1:2], col=col.alpha('blue', 0.3), pch=16) 
points(example[ example$divergent==1, 1:2], col=col.alpha('red', 0.3), pch=16)
mtext( paste0('# divergent = ', sum(example$divergent)), side=3, adj=0, cex=1)

# final sample space
plot(NULL, xlim=range(theta), ylim=range(v), xlab='theta', ylab='v')
for(i in 1:ncol(theta)){
  lines(x=theta[,i], y=v, col=col.alpha('black', 0.3))
}
points(example[ example$divergent==0, 3:2], col=col.alpha('blue', 0.3), pch=16) 
points(example[ example$divergent==1, 3:2], col=col.alpha('red', 0.3), pch=16)
mtext( paste0('# divergent = ', sum(example$divergent)), side=3, adj=0, cex=1)

par(mfrow=c(1, 1))

dev.off()





# Chapter 4 ####