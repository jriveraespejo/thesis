# preliminar ####

# erase all objects
rm(list=ls())

# load libraries
librerias <- c('stringr','dplyr','ggplot2','ggpubr','knitr','tidyverse',
               'reshape2','tinytex','gt','haven','xtable',
               'dagitty','ellipse','mvtnorm','MASS','splines','gtools',
               'rethinking','rstan','coda','runjags','rjags', #'loo',
               'cmdstanr','posterior','bayesplot')
sapply(librerias, require, character.only=T)
# sapply(librerias, install.packages, character.only=T)


# set working directory
setwd('~/Desktop/simulations/#final')


# extras
source('3_functions_extra.R')
opar = par()




# Chapter 3 ####

file_save = file.path(getwd(), 'figures1')



## 1. informative priors ####
n = 1000
set.seed(1235)

# example 1
theta = rnorm(n, 0, 100)
p1 = inv_logit(theta)

# example 2
sigma = rlnorm(n, 0, 3)
theta = rnorm(n, 0, sigma)
p2 = inv_logit(theta)

# example 3
theta = rnorm(n, 0, 1)
p3 = inv_logit(theta)

# example 4
sigma = rlnorm(n, 0, 0.5)
theta = rnorm(n, 0, sigma)
p4 = inv_logit(theta)


png(file.path(file_save, 'prior_elicitation.png'), units='cm', width=25, height=12, res=100)

par(mfrow=c(1,2))
dens( p1 , adj=0.2, xlab='Probability', main='(A)' )
dens( p3 , adj=0.2, lty=3, lwd=1.5, col='blue', add=T )

dens( p2 , adj=0.2, xlab='Probability', main='(B)')
dens( p4 , adj=0.2, lty=3, lwd=1.5, col='blue', add=T )
par(mfrow=c(1,1))

dev.off()





## 2. Devil's funnel ####

### 2.1 model ####

#### stan ####
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
# mod$sample(data=list(y=c(1,0,1,0,0)), 
#            output_dir=file.path(file_save), output_basename='1_stan_CE_simple',
#            chains=3, parallel_chains=3, seed=1 )

# load and transform simulations
fit_files = file.path(file_save, 
                      c('1_stan_CE_simple-1.csv',
                        '1_stan_CE_simple-2.csv',
                        '1_stan_CE_simple-3.csv') )
stan_CE = rstan::read_stan_csv(fit_files)
precis(stan_CE, depth=4)
# notice the small n-eff



#### jags ####
mcmc_code <- "
model{
    v ~ dnorm(0,3)
    theta ~ dnorm(0, exp(v))
}
"
save_code = '1_jags_CE_simple.txt'
writeLines(mcmc_code, con=file.path(file_save, save_code) )

nchains <- 3
nadapt <- 5000
nburn <- 5000
niter <- 5000
parameters <- c('v','theta')

jags_CE = run.jags( model = file.path(file_save, save_code) ,
                    monitor = parameters ,
                    data = list(y=1) ,
                    n.chains = nchains ,
                    adapt = nadapt ,
                    burnin = nburn ,
                    sample = niter)
summary(jags_CE)[,c('Mean','SD','Lower95','Upper95','SSeff','psrf')]



### 2.2 figures ####

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

png(file.path(file_save, '1_jags_CE_simple.png'), units='cm', width=25, height=20, res=100)
plot(jags_CE$mcmc)
dev.off()



stan_samples_CE = stan_CE@sim$samples[[1]]
stan_example = data.frame(stan_samples_CE[,1:2])
stan_example$divergent = attr(stan_samples_CE, 'sampler_params')$divergent__

jags_example = data.frame(unlist(jags_CE$mcmc[[1]]))

sd=3
v = seq(-sd*qnorm(0.975), sd*qnorm(0.975), by=0.01)
z = seq(-2, 2, by=0.2)
theta = sapply(z, function(i){ exp(v)*i } )

png(file.path(file_save, '1_funnel_CE_simple.png'), units='cm', width=35, height=15, res=100)

par(mfrow=c(1, 2))

# stan exploration
plot(NULL, xlim=range(theta), ylim=range(v), xlab='theta', ylab='v', main='(A)')
for(i in 1:ncol(theta)){
  lines(x=theta[,i], y=v, col=col.alpha('black', 0.3), pch=19)
}
points(stan_example[ stan_example$divergent==0, 1:2], col=col.alpha('blue', 0.3), pch=16) 
points(stan_example[ stan_example$divergent==1, 1:2], col=col.alpha('red', 0.3), pch=16) 
mtext( paste0('# divergent = ', sum(stan_example$divergent)), side=3, adj=0, cex=1)

# jags exploration
plot(NULL, xlim=range(theta), ylim=range(v), xlab='theta', ylab='v', main='(B)')
for(i in 1:ncol(theta)){
  lines(x=theta[,i], y=v, col=col.alpha('black', 0.3), pch=19)
}
points(jags_example[,2:1], col=col.alpha('blue', 0.1), pch=16)

par(mfrow=c(1, 1))

dev.off()



  
## 3. Devil's funnel vs priors ####

### 3.1 model ####

#### stan ####
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


#### jags ####
mcmc_code <- "
model{
    v ~ dnorm(0,1)
    theta ~ dnorm(0, exp(v))
}
"
save_code = '2_jags_CE_priors.txt'
writeLines(mcmc_code, con=file.path(file_save, save_code) )

nchains <- 3
nadapt <- 5000
nburn <- 5000
niter <- 5000
parameters <- c('v','theta')

jags_CE = run.jags( model = file.path(file_save, save_code) ,
                    monitor = parameters ,
                    data = list(y=1) ,
                    n.chains = nchains ,
                    adapt = nadapt ,
                    burnin = nburn ,
                    sample = niter)
summary(jags_CE)[,c('Mean','SD','Lower95','Upper95','SSeff','psrf')]



### 3.2 figures ####

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

png(file.path(file_save, '2_jags_CE_priors.png'), units='cm', width=25, height=20, res=100)
plot(jags_CE$mcmc)
dev.off()



stan_samples_CE = stan_CE@sim$samples[[1]]
stan_example = data.frame(stan_samples_CE[,1:2])
stan_example$divergent = attr(stan_samples_CE, 'sampler_params')$divergent__

jags_example = data.frame(unlist(jags_CE$mcmc[[1]]))

sd=1
v = seq(-sd*qnorm(0.975), sd*qnorm(0.975), by=0.01)
z = seq(-2, 2, by=0.2)
theta = sapply(z, function(i){ exp(v)*i } )

png(file.path(file_save, '2_funnel_CE_priors.png'), units='cm', width=35, height=15, res=100)

par(mfrow=c(1, 2))

# stan exploration
plot(NULL, xlim=range(theta), ylim=range(v), xlab='theta', ylab='v', main='(A)')
for(i in 1:ncol(theta)){
  lines(x=theta[,i], y=v, col=col.alpha('black', 0.3), pch=19)
}
points(stan_example[ stan_example$divergent==0, 1:2], col=col.alpha('blue', 0.3), pch=16) 
points(stan_example[ stan_example$divergent==1, 1:2], col=col.alpha('red', 0.3), pch=16) 
mtext( paste0('# divergent = ', sum(stan_example$divergent)), side=3, adj=0, cex=1)

# jags exploration
plot(NULL, xlim=range(theta), ylim=range(v), xlab='theta', ylab='v', main='(B)')
for(i in 1:ncol(theta)){
  lines(x=theta[,i], y=v, col=col.alpha('black', 0.3), pch=19)
}
points(jags_example[,2:1], col=col.alpha('blue', 0.1), pch=16)

par(mfrow=c(1, 1))

dev.off()





## 4. Devil's funnel vs NC ####

### 4.1 model ####

#### stan ####

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


#### jags ####
mcmc_code <- "
model{
    v ~ dnorm(0,1)
    ztheta ~ dnorm(0,1)
    theta = v * ztheta
}
"
save_code = '2_jags_NC.txt'
writeLines(mcmc_code, con=file.path(file_save, save_code) )

nchains <- 3
nadapt <- 5000
nburn <- 5000
niter <- 5000
parameters <- c('v','theta','ztheta')

jags_NC = run.jags( model = file.path(file_save, save_code) ,
                    monitor = parameters ,
                    data = list(y=1) ,
                    n.chains = nchains ,
                    adapt = nadapt ,
                    burnin = nburn ,
                    sample = niter)
summary(jags_CE)[,c('Mean','SD','Lower95','Upper95','SSeff','psrf')]




### 4.2 figures ####

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

png(file.path(file_save, '3_jags_NC.png'), units='cm', width=35, height=30, res=100)
plot(jags_NC$mcmc)
dev.off()


stan_samples_NC = stan_NC@sim$samples[[1]]
stan_example = data.frame(stan_samples_NC[,1:3])
stan_example$divergent = attr(stan_samples_NC, 'sampler_params')$divergent__

jags_example = data.frame(unlist(jags_NC$mcmc[[1]]))


sd=3
v = seq(-sd*qnorm(0.975), sd*qnorm(0.975), by=0.01)
z = seq(-2, 2, by=0.2)
theta = sapply(z, function(i){ exp(v)*i } )

png(file.path(file_save, '3_funnel_NC.png'), units='cm', width=35, height=30, res=100)

par(mfrow=c(2, 2))

# stan
# actual sample space
plot(NULL, xlim=range(z), ylim=range(v), xlab='z', ylab='v', main='(A1)')
for ( l in c(0.1,0.3,0.5,0.8,0.99) ){
  lines( ellipse( diag( c(1, 3) ), centre=c(0,0), level=l), 
         col=col.alpha('black', 0.3))
}
points(stan_example[ stan_example$divergent==0, 1:2], col=col.alpha('blue', 0.3), pch=16) 
points(stan_example[ stan_example$divergent==1, 1:2], col=col.alpha('red', 0.3), pch=16)
mtext( paste0('# divergent = ', sum(stan_example$divergent)), side=3, adj=0, cex=1)

# final sample space
plot(NULL, xlim=range(theta), ylim=range(v), xlab='theta', ylab='v', main='(A2)')
for(i in 1:ncol(theta)){
  lines(x=theta[,i], y=v, col=col.alpha('black', 0.3))
}
points(stan_example[ stan_example$divergent==0, 3:2], col=col.alpha('blue', 0.3), pch=16) 
points(stan_example[ stan_example$divergent==1, 3:2], col=col.alpha('red', 0.3), pch=16)
mtext( paste0('# divergent = ', sum(stan_example$divergent)), side=3, adj=0, cex=1)


# jags
# actual sample space
plot(NULL, xlim=range(z), ylim=range(v), xlab='z', ylab='v', main='(B1)')
for ( l in c(0.1,0.3,0.5,0.8,0.99) ){
  lines( ellipse( diag( c(1, 3) ), centre=c(0,0), level=l), 
         col=col.alpha('black', 0.3))
}
points(jags_example[,c(3,1)], col=col.alpha('blue', 0.1), pch=16)

# final sample space
plot(NULL, xlim=range(theta), ylim=range(v), xlab='theta', ylab='v', main='(B2)')
for(i in 1:ncol(theta)){
  lines(x=theta[,i], y=v, col=col.alpha('black', 0.3))
}
points(jags_example[,c(2,1)], col=col.alpha('blue', 0.1), pch=16)

par(mfrow=c(1, 1))

dev.off()





# Chapter 4 ####

## 1.1 prior elicitation ####


### 1.1.1 FOLV ####

data_load = file.path(getwd(), 'data')
file_load = file.path(getwd(), 'chains_prior')
file_save = file.path(getwd(), 'figures2')

fit_files = 'FOLV_CE_J100_l0.95_Ndata1-1.csv'
stan_model = rstan::read_stan_csv( file.path(file_load, fit_files) )

set.seed(45589)
prior_sim = extract.samples( stan_model )
save(prior_sim, file=file.path(file_save, 'FOLV_CE_J100_l0.95_Ndata1_sim.RData') )
# str(prior_sim)



#### a. ICC's and IIC's #####

png(file.path(file_save, 'FOLV_ICC_prior.png'), 
    units='cm', width=30, height=12, res=200)

par(mfrow=c(1,2))

# ICC
ds_ICC = with(prior_sim, sim_ICC(sim_b=b_k, item_id=1) )
ds_ICC_ave = ave_ICC(sim_b=prior_sim$b_k, item_id=1)
ds_ICC_mar = mar_ICC(sim_ICC_object=ds_ICC)

plot(NULL, xlim=c(-3,3), ylim=range( ds_ICC[,2:ncol(ds_ICC)] ),
     xlab='ability', ylab='probability of endorsing item', main='(A)')
for(s in 1:1000){ 
  lines(ds_ICC$theta, ds_ICC[,s+1], col=col.alpha('black',0.05)) 
}
lines( ds_ICC_ave, col='blue', lwd=2 )
lines( ds_ICC_mar, col='red', lwd=2, lty=2 ) 
legend('topleft', c("Simulated", "Average", 'Marginal'),
       col = c('black','blue','red'), lty=c(1,1,2), bty="n")


# IIF
ds_IIF = with(prior_sim, sim_IIF(sim_b=b_k, item_id=1) )
ds_IIF_ave = ave_IIF(sim_b=prior_sim$b_k, item_id=1)
ds_IIF_mar = mar_IIF(sim_IIF_object=ds_IIF)

plot(NULL, xlim=c(-3,3), ylim=range(ds_IIF[,2:ncol(ds_IIF)]),
     xlab='ability', ylab='information', main='(B)')
for(s in 1:1000){ 
  lines(ds_IIF$theta, ds_IIF[,s+1], col=col.alpha('black',0.05)) 
}
lines( ds_IIF_ave, col='blue', lwd=2 )
lines( ds_IIF_mar, col='red', lwd=2, lty=2 ) 

par(mfrow=c(1,1))

dev.off()




#### b. hit rates #####

data_load = file.path(getwd(), 'data')
file_load = file.path(getwd(), 'chains_prior')
file_save = file.path(getwd(), 'figures2')

load(file.path(data_load, 'ListFormat_J100_l0.95_Ndata5.RData'))

ds_prior = sim_hit_rate(d_true=data_post, d_sim=prior_sim)
# names(ds_prior)



# prior hit rates
png(file.path(file_save, 'FOLV_HitRate1.png'), 
    units='cm', width=30, height=20, res=200)

par(mfrow=c(2,2))

# a. individuals
plot(NULL, xlim=c(0, max(ds_prior$IDind$IDj)+1), ylim=c(0,1), 
     xaxt='n', yaxt='n', main='(A)',
     xlab='Individual ID', ylab='Proportion of correct items')
axis(side=1, at=seq(0, max(ds_prior$IDind$IDj), by=20) ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(s in (1:1000)+1){ 
  points(ds_prior$IDind[,c(1,s)], col=col.alpha('black', 0.005), pch=19) 
}


# b. items
plot(NULL, xlim=c(0, max(ds_prior$IDitem$IDk)+1), ylim=c(0,1), 
     xaxt='n', yaxt='n', main='(B)',
     xlab='Item ID', ylab='Proportion of correct individuals')
axis(side=1, at=seq(1, max(ds_prior$IDitem$IDk), by=1) ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
abline(v=c(1:4)*5+0.5, col=col.alpha('black', 0.5), lty=2)
for(s in (1:1000)+1){ 
  points(ds_prior$IDitem[,c(1,s)], col=col.alpha('black', 0.01), pch=19) 
}


# c. text
plot(NULL, xlim=c(0, max(ds_prior$IDtext$IDl)+1), ylim=c(0,1), 
     xaxt='n', yaxt='n', main='(C)',
     xlab='Text ID', ylab='Proportion of correct individuals')
axis(side=1, at=seq(1, max(ds_prior$IDtext$IDl), by=1) ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(s in (1:1000)+1){ 
  points(ds_prior$IDtext[,c(1,s)], col=col.alpha('black', 0.01), pch=19) 
}


# d. dimensions
plot(NULL, xlim=c(0, max(ds_prior$IDdim$IDd)+1), ylim=c(0,1), 
     xaxt='n', yaxt='n', main='(D)',
     xlab='Dimension ID', ylab='Proportion of correct individuals')
axis(side=1, at=seq(1, max(ds_prior$IDdim$IDd), by=1) ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(s in (1:1000)+1){ 
  points(ds_prior$IDdim[,c(1,s)], col=col.alpha('black', 0.01), pch=19) 
}

par(mfrow=c(1,1))

dev.off()




# prior hit rates per variable
ds_prior = merge(ds_prior$IDind, data.frame( data_post[c('IDind','G','A','E','X')] ), 
                 by.x='IDj', by.y='IDind', all.x=T, all.y=F)

png(file.path(file_save, 'FOLV_HitRate2.png'), 
    units='cm', width=30, height=20, res=200)

par(mfrow=c(2,2))

# gender
plot(NULL, xlim=with(ds_prior, c(min(G)-1, max(G)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.05), pch=19,
     xlab='', ylab='% correct items', main='(A)')
axis(side=1, at=with(ds_prior, min(G):max(G)), labels=c('male','female') ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(s in (1:100)+1 ){
  points(ds_prior[,c('G', paste0('X',s) )], col=col.alpha('black', 0.005), pch=19 )
}

# age
plot(NULL, xlim=with(ds_prior, c(min(A)-1, max(A)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.05), pch=19,
     xlab='', ylab='% correct items', main='(B)')
axis(side=1, at=with(ds_prior, min(A):max(A)) ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(s in (1:500)+1 ){
  points(ds_prior[,c('A', paste0('X',s) )], col=col.alpha('black', 0.005), pch=19 )
}

# edu
plot(NULL, xlim=with(ds_prior, c(min(E)-1, max(E)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.05), pch=19,
     xlab='', ylab='% correct items', main='(C)')
axis(side=1, at=with(ds_prior, min(E):max(E)), labels=c('inst','uni','both') ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(s in (1:100)+1 ){
  points(ds_prior[,c('E', paste0('X',s) )], col=col.alpha('black', 0.005), pch=19 )
}

# exp
plot(NULL, xlim=with(ds_prior, c(min(X)-1, max(X)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.05), pch=19,
     xlab='', ylab='% correct items', main='(D)')
axis(side=1, at=with(ds_prior, min(X):max(X)), 
     labels=c('0y','1-5y','6-10y','+10y') ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(s in (1:100)+1 ){
  points(ds_prior[,c('X', paste0('X',s) )], col=col.alpha('black', 0.005), pch=19 )
}

par(mfrow=c(1,1))

dev.off()




### 1.1.2 SOLV ####

data_load = file.path(getwd(), 'data')
file_load = file.path(getwd(), 'chains_prior')
file_save = file.path(getwd(), 'figures2')

fit_files = 'SOLV_CE_J100_l0.95_Ndata1-1.csv'
stan_model = rstan::read_stan_csv( file.path(file_load, fit_files) )

set.seed(45589)
prior_sim = extract.samples( stan_model )
save(prior_sim, file=file.path(file_save, 'SOLV_CE_J100_l0.95_Ndata1_sim.RData') )
# str(prior_sim)



#### a. ICC's and IIC's #####

png(file.path(file_save, 'SOLV_ICC_prior.png'), 
    units='cm', width=30, height=12, res=200)

par(mfrow=c(1,2))

# ICC
ds_ICC = with(prior_sim, sim_ICC(sim_b=b_k, item_id=1) )
ds_ICC_ave = ave_ICC(sim_b=prior_sim$b_k, item_id=1)
ds_ICC_mar = mar_ICC(sim_ICC_object=ds_ICC)

plot(NULL, xlim=c(-3,3), ylim=range( ds_ICC[,2:ncol(ds_ICC)] ),
     xlab='ability', ylab='probability of endorsing item', main='(A)')
for(s in 1:1000){ 
  lines(ds_ICC$theta, ds_ICC[,s+1], col=col.alpha('black',0.05)) 
}
lines( ds_ICC_ave, col='blue', lwd=2 )
lines( ds_ICC_mar, col='red', lwd=2, lty=2 ) 
legend('topleft', c("Simulated", "Average", 'Marginal'),
       col = c('black','blue','red'), lty=c(1,1,2), bty="n")


# IIF
ds_IIF = with(prior_sim, sim_IIF(sim_b=b_k, item_id=1) )
ds_IIF_ave = ave_IIF(sim_b=prior_sim$b_k, item_id=1)
ds_IIF_mar = mar_IIF(sim_IIF_object=ds_IIF)

plot(NULL, xlim=c(-3,3), ylim=range(ds_IIF[,2:ncol(ds_IIF)]),
     xlab='ability', ylab='information', main='(B)')
for(s in 1:1000){ 
  lines(ds_IIF$theta, ds_IIF[,s+1], col=col.alpha('black',0.05)) 
}
lines( ds_IIF_ave, col='blue', lwd=2 )
lines( ds_IIF_mar, col='red', lwd=2, lty=2 ) 

par(mfrow=c(1,1))

dev.off()




#### b. hit rates #####

data_load = file.path(getwd(), 'data')
file_load = file.path(getwd(), 'chains_prior')
file_save = file.path(getwd(), 'figures2')

load(file.path(data_load, 'ListFormat_J100_l0.95_Ndata5.RData'))

ds_prior = sim_hit_rate(d_true=data_post, d_sim=prior_sim)
# names(ds_prior)



# prior hit rates
png(file.path(file_save, 'SOLV_HitRate1.png'), 
    units='cm', width=30, height=20, res=200)

par(mfrow=c(2,2))

# a. individuals
plot(NULL, xlim=c(0, max(ds_prior$IDind$IDj)+1), ylim=c(0,1), 
     xaxt='n', yaxt='n', main='(A)',
     xlab='Individual ID', ylab='Proportion of correct items')
axis(side=1, at=seq(0, max(ds_prior$IDind$IDj), by=20) ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(s in (1:1000)+1){ 
  points(ds_prior$IDind[,c(1,s)], col=col.alpha('black', 0.005), pch=19) 
}


# b. items
plot(NULL, xlim=c(0, max(ds_prior$IDitem$IDk)+1), ylim=c(0,1), 
     xaxt='n', yaxt='n', main='(B)',
     xlab='Item ID', ylab='Proportion of correct individuals')
axis(side=1, at=seq(1, max(ds_prior$IDitem$IDk), by=1) ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
abline(v=c(1:4)*5+0.5, col=col.alpha('black', 0.5), lty=2)
for(s in (1:1000)+1){ 
  points(ds_prior$IDitem[,c(1,s)], col=col.alpha('black', 0.01), pch=19) 
}


# c. text
plot(NULL, xlim=c(0, max(ds_prior$IDtext$IDl)+1), ylim=c(0,1), 
     xaxt='n', yaxt='n', main='(C)',
     xlab='Text ID', ylab='Proportion of correct individuals')
axis(side=1, at=seq(1, max(ds_prior$IDtext$IDl), by=1) ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(s in (1:1000)+1){ 
  points(ds_prior$IDtext[,c(1,s)], col=col.alpha('black', 0.01), pch=19) 
}


# d. dimensions
plot(NULL, xlim=c(0, max(ds_prior$IDdim$IDd)+1), ylim=c(0,1), 
     xaxt='n', yaxt='n', main='(D)',
     xlab='Dimension ID', ylab='Proportion of correct individuals')
axis(side=1, at=seq(1, max(ds_prior$IDdim$IDd), by=1) ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(s in (1:1000)+1){ 
  points(ds_prior$IDdim[,c(1,s)], col=col.alpha('black', 0.01), pch=19) 
}

par(mfrow=c(1,1))

dev.off()




# prior hit rates per variable
ds_prior = merge(ds_prior$IDind, data.frame( data_post[c('IDind','G','A','E','X')] ), 
                 by.x='IDj', by.y='IDind', all.x=T, all.y=F)

png(file.path(file_save, 'SOLV_HitRate2.png'), 
    units='cm', width=30, height=20, res=200)

par(mfrow=c(2,2))

# gender
plot(NULL, xlim=with(ds_prior, c(min(G)-1, max(G)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.05), pch=19,
     xlab='', ylab='% correct items', main='(A)')
axis(side=1, at=with(ds_prior, min(G):max(G)), labels=c('male','female') ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(s in (1:100)+1 ){
  points(ds_prior[,c('G', paste0('X',s) )], col=col.alpha('black', 0.005), pch=19 )
}

# age
plot(NULL, xlim=with(ds_prior, c(min(A)-1, max(A)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.05), pch=19,
     xlab='', ylab='% correct items', main='(B)')
axis(side=1, at=with(ds_prior, min(A):max(A)) ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(s in (1:500)+1 ){
  points(ds_prior[,c('A', paste0('X',s) )], col=col.alpha('black', 0.005), pch=19 )
}

# edu
plot(NULL, xlim=with(ds_prior, c(min(E)-1, max(E)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.05), pch=19,
     xlab='', ylab='% correct items', main='(C)')
axis(side=1, at=with(ds_prior, min(E):max(E)), labels=c('inst','uni','both') ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(s in (1:100)+1 ){
  points(ds_prior[,c('E', paste0('X',s) )], col=col.alpha('black', 0.005), pch=19 )
}

# exp
plot(NULL, xlim=with(ds_prior, c(min(X)-1, max(X)+1)), ylim=c(0,1), 
     xaxt='n', yaxt='n', col=col.alpha('black', 0.05), pch=19,
     xlab='', ylab='% correct items', main='(D)')
axis(side=1, at=with(ds_prior, min(X):max(X)), 
     labels=c('0y','1-5y','6-10y','+10y') ) 
axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
for(s in (1:100)+1 ){
  points(ds_prior[,c('X', paste0('X',s) )], col=col.alpha('black', 0.005), pch=19 )
}

par(mfrow=c(1,1))

dev.off()





## 1.2 Performance ####

# previous models
chains_path = file.path(getwd(), 'chains_post')
models_int = c('FOLV_CE','FOLV_NC', 'SOLV_CE', 'SOLV_NC')
chains_list = file_id(chains_path, models_int)
idx = str_detect(chains_list$model, '_mod')
chains_list = chains_list[!idx, ]

### trace, trank, acf plots ####

# plots
figures_plot(c_list=chains_list,
             chains_path=chains_path,
             file_save=file.path(getwd(), 'figures3'))

stat_chain(c_list = chains_list,
           chains_path = chains_path,
           file_save = file.path(getwd(), 'figures4'),
           file_name = 'stan_stats_no_mod',
           contr_pars = c('b_G','b_E','b_X'))



# # extract statistics
# chains_path = file.path(getwd(), 'chains_post')
# models_int = c('FOLV_CE_mod','FOLV_CE','FOLV_NC_mod','FOLV_NC', 'SOLV_CE', 'SOLV_NC')
# chains_list = file_id(chains_path, models_int)
# stat_chain(c_list = chains_list,
#            chains_path = chains_path,
#            file_save = file.path(getwd(), 'figures4'),
#            file_name = 'stan_stats_with_mod',
#            contr_pars = c('b_G','b_E','b_X'))
# load( file.path(getwd(), 'figures4', 'stan_stats.RData') )
# idx = with(stan_int, model_type=='FOLV_NC_mod' & 
#              parameter %in% c('Rho_theta_sub[1,2]', 'Rho_theta_sub[1,3]', 'Rho_theta_sub[2,3]'))
# names(stan_int)[5] = 'mean1'
# stan_int[idx,] %>%
#   summarise(mean=mean(mean1), sd=sd(mean1))
# # not enough to justify changes to the document




### neff ####

load( file.path(getwd(), 'figures4', 'stan_stats.RData') )
file_save = file.path(getwd(), 'figures4')


# FOLV
for(ss in c(100, 250, 500)){
  
  png(file.path(file_save, paste0('FOLV_',ss,'_neff1.png')), 
      units='cm', width=30, height=12, res=200)
  par(mfrow=c(1, 2))
  plot_stat(stat_object = stan_int, info = 'n_eff', 
            par_int = c( paste0('m_b[',1:5,']'), paste0('s_b[',1:5,']')),
            model = 'FOLV', ssize = ss, title='(A)')
  
  plot_stat(stat_object = stan_int, info = 'n_eff', 
            par_int = paste0('b_k[',1:25,']'),
            model = 'FOLV', ssize = ss, title='(B)')
  par(mfrow=c(1, 1))
  dev.off()
  
  
  png(file.path(file_save, paste0('FOLV_',ss,'_neff2.png')), 
      units='cm', width=30, height=12, res=200)
  par(mfrow=c(1, 2))
  plot_stat(stat_object = stan_int, info = 'n_eff', 
            par_int = c( paste0('b_G[',1:2,']'),'b_A', paste0('b_E[',1:3,']'),
                         paste0('b_X[',1:4,']')),
            model = 'FOLV', ssize = ss, title='(A)')
  
  plot_stat(stat_object = stan_int, info = 'n_eff', 
            par_int = c('Rho_theta_sub[1,2]','Rho_theta_sub[1,3]','Rho_theta_sub[2,3]'),
            model = 'FOLV', ssize = ss, title='(B)')
  par(mfrow=c(1, 1))
  dev.off()
  
  
  png(file.path(file_save, paste0('FOLV_',ss,'_neff3.png')), 
      units='cm', width=30, height=12, res=200)
  par(mfrow=c(1, 2))
  plot_stat(stat_object = stan_int, info = 'n_eff', 
            par_int = paste0('theta_sub[', 1:ss, ',1]'),
            model = 'FOLV', ssize = ss, title='(A)')
  
  plot_stat(stat_object = stan_int, info = 'n_eff', 
            par_int = paste0('theta_sub[', 1:ss, ',2]'),
            model = 'FOLV', ssize = ss, title='(B)')
  par(mfrow=c(1, 1))
  dev.off()
  
  
  png(file.path(file_save, paste0('FOLV_',ss,'_neff4.png')), 
      units='cm', width=15, height=12, res=200)
  par(mfrow=c(1, 1))
  plot_stat(stat_object = stan_int, info = 'n_eff', 
            par_int = paste0('theta_sub[', 1:ss, ',3]'),
            model = 'FOLV', ssize = ss, title='(A)')
  par(mfrow=c(1, 1))
  dev.off()
  
}



# SOLV
for(ss in c(100, 250, 500)){
  
  png(file.path(file_save, paste0('SOLV_',ss,'_neff1.png')), 
      units='cm', width=30, height=12, res=200)
  par(mfrow=c(1, 2))
  plot_stat(stat_object = stan_int, info = 'n_eff', 
            par_int = c( paste0('m_b[',1:5,']'), paste0('s_b[',1:5,']')),
            model = 'SOLV', ssize = ss, title='(A)')
  plot_stat(stat_object = stan_int, info = 'n_eff', 
            par_int = paste0('b_k[',1:25,']'),
            model = 'SOLV', ssize = ss, title='(B)')
  par(mfrow=c(1, 1))
  dev.off()
  
  
  png(file.path(file_save, paste0('SOLV_',ss,'_neff2.png')), 
      units='cm', width=15, height=12, res=200)
  par(mfrow=c(1, 1))
  plot_stat(stat_object = stan_int, info = 'n_eff', 
            par_int = c( paste0('b_G[',1:2,']'),'b_A', paste0('b_E[',1:3,']'),
                         paste0('b_X[',1:4,']')),
            model = 'SOLV', ssize = ss, title='(A)')
  par(mfrow=c(1, 1))
  dev.off()
  
  png(file.path(file_save, paste0('SOLV_',ss,'_neff3.png')), 
      units='cm', width=30, height=12, res=200)
  par(mfrow=c(1, 2))
  plot_stat(stat_object = stan_int, info = 'n_eff', 
            par_int = c('Rho_theta_sub[1,2]','Rho_theta_sub[1,3]','Rho_theta_sub[2,3]'),
            model = 'SOLV', ssize = ss, title='(A)')
  
  plot_stat(stat_object = stan_int, info = 'n_eff', 
            par_int = paste0('loads[',1:3,']'),
            model = 'SOLV', ssize = ss, title='(B)')
  par(mfrow=c(1, 1))
  dev.off()
  
  png(file.path(file_save, paste0('SOLV_',ss,'_neff4.png')), 
      units='cm', width=30, height=12, res=200)
  par(mfrow=c(1, 2))
  plot_stat(stat_object = stan_int, info = 'n_eff', 
            par_int = paste0('theta_sub[', 1:ss, ',1]'),
            model = 'SOLV', ssize = ss, title='(A)')
  
  plot_stat(stat_object = stan_int, info = 'n_eff', 
            par_int = paste0('theta_sub[', 1:ss, ',2]'),
            model = 'SOLV', ssize = ss, title='(B)')
  par(mfrow=c(1, 1))
  dev.off()
  
  
  png(file.path(file_save, paste0('SOLV_',ss,'_neff5.png')), 
      units='cm', width=15, height=12, res=200)
  par(mfrow=c(1, 1))
  plot_stat(stat_object = stan_int, info = 'n_eff', 
            par_int = paste0('theta_sub[', 1:ss, ',3]'),
            model = 'SOLV', ssize = ss, title='(A)')
  par(mfrow=c(1, 1))
  dev.off()
  
  
  png(file.path(file_save, paste0('SOLV_',ss,'_neff6.png')), 
      units='cm', width=15, height=12, res=200)
  par(mfrow=c(1, 1))
  plot_stat(stat_object = stan_int, info = 'n_eff', 
            par_int = paste0('theta[', 1:ss, ']'),
            model = 'SOLV', ssize = ss, title='(A)')
  par(mfrow=c(1, 1))
  dev.off()
  
}



### Rhat #####

load( file.path(getwd(), 'figures4', 'stan_stats.RData') )
file_save = file.path(getwd(), 'figures4')

# FOLV
for(ss in c(100, 250, 500)){
  
  png(file.path(file_save, paste0('FOLV_',ss,'_Rhat1.png')), 
      units='cm', width=30, height=12, res=200)
  par(mfrow=c(1, 2))
  plot_stat(stat_object = stan_int, info = 'Rhat4', 
            par_int = c( paste0('m_b[',1:5,']'), paste0('s_b[',1:5,']')),
            model = 'FOLV', ssize = ss, title='(C)')
  plot_stat(stat_object = stan_int, info = 'Rhat4', 
            par_int = paste0('b_k[',1:25,']'),
            model = 'FOLV', ssize = ss, title='(D)')
  par(mfrow=c(1, 1))
  dev.off()
  
  
  png(file.path(file_save, paste0('FOLV_',ss,'_Rhat2.png')), 
      units='cm', width=30, height=12, res=200)
  par(mfrow=c(1, 2))
  plot_stat(stat_object = stan_int, info = 'Rhat4', 
            par_int = c( paste0('b_G[',1:2,']'),'b_A', paste0('b_E[',1:3,']'),
                         paste0('b_X[',1:4,']')),
            model = 'FOLV', ssize = ss, title='(C)')
  
  plot_stat(stat_object = stan_int, info = 'Rhat4', 
            par_int = c('Rho_theta_sub[1,2]','Rho_theta_sub[1,3]','Rho_theta_sub[2,3]'),
            model = 'FOLV', ssize = ss, title='(D)')
  par(mfrow=c(1, 1))
  dev.off()
  
  
  png(file.path(file_save, paste0('FOLV_',ss,'_Rhat3.png')), 
      units='cm', width=30, height=12, res=200)
  par(mfrow=c(1, 2))
  plot_stat(stat_object = stan_int, info = 'Rhat4', 
            par_int = paste0('theta_sub[', 1:ss, ',1]'),
            model = 'FOLV', ssize = ss, title='(C)')
  
  plot_stat(stat_object = stan_int, info = 'Rhat4', 
            par_int = paste0('theta_sub[', 1:ss, ',2]'),
            model = 'FOLV', ssize = ss, title='(D)')
  par(mfrow=c(1, 1))
  dev.off()
  
  
  png(file.path(file_save, paste0('FOLV_',ss,'_Rhat4.png')), 
      units='cm', width=15, height=12, res=200)
  par(mfrow=c(1, 1))
  plot_stat(stat_object = stan_int, info = 'Rhat4', 
            par_int = paste0('theta_sub[', 1:ss, ',3]'),
            model = 'FOLV', ssize = ss, title='(B)')
  par(mfrow=c(1, 1))
  dev.off()
  
}



# SOLV
for(ss in c(100, 250, 500)){
  
  png(file.path(file_save, paste0('SOLV_',ss,'_Rhat1.png')), 
      units='cm', width=30, height=12, res=200)
  par(mfrow=c(1, 2))
  plot_stat(stat_object = stan_int, info = 'Rhat4', 
            par_int = c( paste0('m_b[',1:5,']'), paste0('s_b[',1:5,']')),
            model = 'SOLV', ssize = ss, title='(C)')
  
  plot_stat(stat_object = stan_int, info = 'Rhat4', 
            par_int = paste0('b_k[',1:25,']'),
            model = 'SOLV', ssize = ss, title='(D)')
  par(mfrow=c(1, 1))
  dev.off()
  
  
  png(file.path(file_save, paste0('SOLV_',ss,'_Rhat2.png')), 
      units='cm', width=15, height=12, res=200)
  par(mfrow=c(1, 1))
  plot_stat(stat_object = stan_int, info = 'Rhat4', 
            par_int = c( paste0('b_G[',1:2,']'),'b_A', paste0('b_E[',1:3,']'),
                         paste0('b_X[',1:4,']')),
            model = 'SOLV', ssize = ss, title='(B)')
  par(mfrow=c(1, 1))
  dev.off()
  
  
  png(file.path(file_save, paste0('SOLV_',ss,'_Rhat3.png')), 
      units='cm', width=30, height=12, res=200)
  par(mfrow=c(1, 2))
  plot_stat(stat_object = stan_int, info = 'Rhat4', 
            par_int = c('Rho_theta_sub[1,2]','Rho_theta_sub[1,3]','Rho_theta_sub[2,3]'),
            model = 'SOLV', ssize = ss, title='(C)')
  
  plot_stat(stat_object = stan_int, info = 'Rhat4', 
            par_int = paste0('loads[',1:3,']'),
            model = 'SOLV', ssize = ss, title='(D)')
  par(mfrow=c(1, 1))
  dev.off()
  
  
  png(file.path(file_save, paste0('SOLV_',ss,'_Rhat4.png')), 
      units='cm', width=30, height=12, res=200)
  par(mfrow=c(1, 2))
  plot_stat(stat_object = stan_int, info = 'Rhat4', 
            par_int = paste0('theta_sub[', 1:ss, ',1]'),
            model = 'SOLV', ssize = ss, title='(C)')
  
  plot_stat(stat_object = stan_int, info = 'Rhat4', 
            par_int = paste0('theta_sub[', 1:ss, ',2]'),
            model = 'SOLV', ssize = ss, title='(D)')
  par(mfrow=c(1, 1))
  dev.off()
  
  
  png(file.path(file_save, paste0('SOLV_',ss,'_Rhat5.png')), 
      units='cm', width=15, height=12, res=200)
  par(mfrow=c(1, 1))
  plot_stat(stat_object = stan_int, info = 'Rhat4', 
            par_int = paste0('theta_sub[', 1:ss, ',3]'),
            model = 'SOLV', ssize = ss, title='(B)')
  par(mfrow=c(1, 1))
  dev.off()
  
  
  png(file.path(file_save, paste0('SOLV_',ss,'_Rhat6.png')), 
      units='cm', width=15, height=12, res=200)
  par(mfrow=c(1, 1))
  plot_stat(stat_object = stan_int, info = 'Rhat4', 
            par_int = paste0('theta[', 1:ss, ']'),
            model = 'SOLV', ssize = ss, title='(B)')
  par(mfrow=c(1, 1))
  dev.off()
  
}





## 1.3 Recovery ####

### 1.3.1 tables ####

#### calculations ####

# load statistics 
load( file.path(getwd(), 'figures4', 'stan_stats.RData') )
file_save = file.path(getwd(), 'figures5')

# # extract and load true parameters
# extract_true(file_load = file.path(getwd(), 'data'), 
#              file_save = file_save)
load( file.path(getwd(), 'figures5', 'true_pars.RData') )
# true_pars[true_pars$sample_size==250,]



# calculate and load RMSE
# rmse_pars( stats_object = stan_int,
#            true_object = true_pars,
#            file_save = file_save)
load( file.path(getwd(), 'figures5', 'rmse_pars.RData') )



#### regression and contrast ####
par_int = c('a','b_G','b_A','b_E','b_X') 
for(i in 1:length(par_int)){
  idx_mom = str_detect(stan_rmse$parameter, paste0('^', par_int[i]) )
  if(i==1){
    idx = idx_mom
  } else {
    idx = idx | idx_mom
  }
}
idx1 = idx & str_detect(stan_rmse$model_type, '^FOLV')
idx2 = idx & str_detect(stan_rmse$model_type, '^SOLV')

stan_rmse[idx1,] %>%
  spread(key=sample_size, value=RMSE) %>%
  xtable( caption='caption', label = 'tab:', digits=3)

stan_rmse[idx2,] %>%
  spread(key=sample_size, value=RMSE) %>%
  xtable( caption='caption', label = 'tab:')



#### correlations ####
idx = stan_rmse$parameter %in% c('Rho_theta_sub[1,2]', 'Rho_theta_sub[1,3]','Rho_theta_sub[2,3]' )
idx1 = idx & str_detect(stan_rmse$model_type, '^FOLV')
idx2 = idx & str_detect(stan_rmse$model_type, '^SOLV')

stan_rmse[idx1,] %>%
  spread(key=sample_size, value=RMSE) %>%
  xtable( caption='caption', label = 'tab:', digits=3)

stan_rmse[idx2,] %>%
  spread(key=sample_size, value=RMSE) %>%
  xtable( caption='caption', label = 'tab:')


#### loadings ####
idx = stan_rmse$parameter %in% paste0('loads[',1:3,']')
idx2 = idx & str_detect(stan_rmse$model_type, '^SOLV')

stan_rmse[idx2,] %>%
  spread(key=sample_size, value=RMSE) %>%
  xtable( caption='caption', label = 'tab:')


#### texts ####
idx = str_detect(stan_rmse$parameter, 'm_b')
idx1 = idx & str_detect(stan_rmse$model_type, '^FOLV')
idx2 = idx & str_detect(stan_rmse$model_type, '^SOLV')

stan_rmse[idx1,] %>%
  spread(key=sample_size, value=RMSE) %>%
  xtable( caption='caption', label = 'tab:', digits=3)

stan_rmse[idx2,] %>%
  spread(key=sample_size, value=RMSE) %>%
  xtable( caption='caption', label = 'tab:')



idx = str_detect(stan_rmse$parameter, 's_b')
idx1 = idx & str_detect(stan_rmse$model_type, '^FOLV')
idx2 = idx & str_detect(stan_rmse$model_type, '^SOLV')

stan_rmse[idx1,] %>%
  spread(key=sample_size, value=RMSE) %>%
  xtable( caption='caption', label = 'tab:', digits=3)

stan_rmse[idx2,] %>%
  spread(key=sample_size, value=RMSE) %>%
  xtable( caption='caption', label = 'tab:')




#### items ####
idx = str_detect(stan_rmse$parameter, 'b_k')
idx1 = idx & str_detect(stan_rmse$model_type, '^FOLV')
idx2 = idx & str_detect(stan_rmse$model_type, '^SOLV')

stan_rmse[idx1,] %>%
  group_by(model_type, sample_size) %>%
  summarise(n=n(), mean=mean(RMSE), sd=sd(RMSE), min=min(RMSE), max=max(RMSE)) %>%
  xtable( caption='caption', label = 'tab:', digits=3)

stan_rmse[idx2,] %>%
  group_by(model_type, sample_size) %>%
  summarise(n=n(), mean=mean(RMSE), sd=sd(RMSE), min=min(RMSE), max=max(RMSE)) %>%
  xtable( caption='caption', label = 'tab:', digits=3)



#### abilities ####
idx = stan_rmse$parameter %in% paste0('theta_sub[',1:500,',1]')
idx1 = idx & str_detect(stan_rmse$model_type, '^FOLV')
idx2 = idx & str_detect(stan_rmse$model_type, '^SOLV')

stan_rmse[idx1,] %>%
  group_by(model_type, sample_size) %>%
  summarise(mean=mean(RMSE), sd=sd(RMSE), min=min(RMSE), max=max(RMSE)) %>%
  xtable( caption='caption', label = 'tab:', digits=3)

stan_rmse[idx2,] %>%
  group_by(model_type, sample_size) %>%
  summarise(mean=mean(RMSE), sd=sd(RMSE), min=min(RMSE), max=max(RMSE)) %>%
  xtable( caption='caption', label = 'tab:', digits=3)



idx = stan_rmse$parameter %in% paste0('theta_sub[',1:500,',2]')
idx1 = idx & str_detect(stan_rmse$model_type, '^FOLV')
idx2 = idx & str_detect(stan_rmse$model_type, '^SOLV')

stan_rmse[idx1,] %>%
  group_by(model_type, sample_size) %>%
  summarise(mean=mean(RMSE), sd=sd(RMSE), min=min(RMSE), max=max(RMSE)) %>%
  xtable( caption='caption', label = 'tab:', digits=3)

stan_rmse[idx2,] %>%
  group_by(model_type, sample_size) %>%
  summarise(mean=mean(RMSE), sd=sd(RMSE), min=min(RMSE), max=max(RMSE)) %>%
  xtable( caption='caption', label = 'tab:', digits=3)



idx = stan_rmse$parameter %in% paste0('theta_sub[',1:500,',3]')
idx1 = idx & str_detect(stan_rmse$model_type, '^FOLV')
idx2 = idx & str_detect(stan_rmse$model_type, '^SOLV')

stan_rmse[idx1,] %>%
  group_by(model_type, sample_size) %>%
  summarise(mean=mean(RMSE), sd=sd(RMSE), min=min(RMSE), max=max(RMSE)) %>%
  xtable( caption='caption', label = 'tab:', digits=3)

stan_rmse[idx2,] %>%
  group_by(model_type, sample_size) %>%
  summarise(mean=mean(RMSE), sd=sd(RMSE), min=min(RMSE), max=max(RMSE)) %>%
  xtable( caption='caption', label = 'tab:', digits=3)



idx = stan_rmse$parameter %in% paste0('theta[',1:500,']')
idx2 = idx & str_detect(stan_rmse$model_type, '^SOLV')

stan_rmse[idx2,] %>%
  group_by(model_type, sample_size) %>%
  summarise(mean=mean(RMSE), sd=sd(RMSE), min=min(RMSE), max=max(RMSE)) %>%
  xtable( caption='caption', label = 'tab:', digits=3)




### 1.3.2 recovery plots ####
file_save = file.path(getwd(), 'figures5')
load( file.path(file_save, 'results_pars.RData') )

recovery_plots(result_object = res_stan,
               figure_save = file_save)




# correlation plots ####

chains_path = file.path(getwd(), 'chains_post')
models_int = c('SOLV_CE', 'SOLV_NC')
chains_list = file_id(chains_path, models_int)
idx = str_detect(chains_list$model, '_mod')
chains_list = chains_list[!idx, ]

post_corr(c_list = chains_list,
          file_save = file.path(getwd(), 'figures5'))
  




## 1.4 Retrodictive accuracy ####



### 1.4.1 tables ####

chains_path = file.path(getwd(), 'chains_post')
models_int = c('FOLV_CE', 'FOLV_NC', 'SOLV_CE', 'SOLV_NC')
chains_list = file_id(chains_path, models_int)

# # generate aggregation
# # (it takes a lot of time!!)
# retrodictive_agg(c_list = chains_list,
#                  file_save = file.path(getwd(), 'figures6'),
#                  true_load = file.path(getwd(), 'data'))

# # generate within between rmse
# retrodictive_rmse(c_list = chains_list,
#                   file_save = file.path(getwd(), 'figures6'),
#                   true_load = file.path(getwd(), 'figures6'))


load( file.path(getwd(), 'figures6', 'rmse_WB_IDind.RData') )
rmse_final[rmse_final$type=='out',] %>%
  group_by(model_type,sample_size) %>%
  summarise(mean_W=mean(mean_RMSE_within), 
            min_W=min(mean_RMSE_within), 
            max_W=max(mean_RMSE_within),
            mean_B=mean(RMSE_between), 
            min_B=min(RMSE_between), 
            max_B=max(RMSE_between)) %>%
  xtable( caption='caption', label = 'tab:', digits = 3)  


load( file.path(getwd(), 'figures6', 'rmse_WB_IDitems.RData') )
rmse_final[rmse_final$type=='out',] %>%
  group_by(model_type,sample_size) %>%
  summarise(mean_W=mean(mean_RMSE_within), 
            min_W=min(mean_RMSE_within), 
            max_W=max(mean_RMSE_within),
            mean_B=mean(RMSE_between), 
            min_B=min(RMSE_between), 
            max_B=max(RMSE_between)) %>%
  xtable( caption='caption', label = 'tab:', digits=3)


load( file.path(getwd(), 'figures6', 'rmse_WB_IDtext.RData') )
rmse_final[rmse_final$type=='out',]

load( file.path(getwd(), 'figures6', 'rmse_WB_IDdim.RData') )
rmse_final[rmse_final$type=='out',] 




### 1.4.2 plots ####

chains_path = file.path(getwd(), 'chains_post')
models_int = c('FOLV_CE', 'FOLV_NC', 'SOLV_CE', 'SOLV_NC')
chains_list = file_id(chains_path, models_int)

# # ICC and IIF
# # (it takes a lot of time!!)
# item_plots(c_list = chains_list,
#            file_save = file.path(getwd(), 'figures6'),
#            true_load = file.path(getwd(), 'data'))

retrodictive_plots(c_list = chains_list,
                   file_save = file.path(getwd(), 'figures6'),
                   true_load = file.path(getwd(), 'figures6'))








## 1.5 time ####

file_load = file.path(getwd(), 'chains_post')
times = read.csv( file.path(file_load, 'time_elapsed.csv'), stringsAsFactors=F)
times$time_min = times$time/60

times %>%
  group_by(Model, J) %>%
  summarize(mean=mean(time_min), min=min(time_min), max=max(time_min)) %>%
  xtable( caption='caption', label = 'tab:')
