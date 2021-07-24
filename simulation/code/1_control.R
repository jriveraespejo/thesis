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


# loading functions
source('2_functions_sim.R')



# 1. data generating ####
S = 10
condition = expand_grid( J = c(100, 250, 500), load=0.95)
for(i in 1:nrow(condition)){
  for(s in 1:S){
    with(condition[i,],
         data_generation( J=J, loads=rep(load, 3), Ndata=s, seed=4587+s+i,
                          file_dir=file.path(getwd(), 'data') ) )
  }
}



# 2. priors ####
source('4_models_prior.R')
run_prior( model_path = file.path(getwd(), 'models_prior'),
           model_out = file.path(getwd(), 'chains_prior'),
           data_path = file.path(getwd(), 'data') )

# # (to test)
# fitfile = 'FOLV_CE_J100_l0.6_S1-1.csv'
# stan_model = rstan::read_stan_csv( file.path(getwd(), 'priors', fitfile) )
# traceplot_ulam(stan_model)



# 3. posteriors ####
source('5_models_post.R')
run_post( model_path = file.path(getwd(), 'models_post'),
          model_out = file.path(getwd(), 'chains_post'),
          data_path = file.path(getwd(), 'data') )

# # (to test)
# fitfile = 'FOLV_CE_J100_l0.6_S1-1.csv'
# stan_model = rstan::read_stan_csv( file.path(getwd(), 'priors', fitfile) )
# traceplot_ulam(stan_model)