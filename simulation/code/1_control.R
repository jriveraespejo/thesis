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

model_path = file.path(getwd(), 'models_prior')
model_list = list.files( model_path )
model_list = model_list[ str_detect(model_list, '.stan') ]

run_prior( models = model_list,
           model_path = model_path,
           model_out = file.path(getwd(), 'chains_prior'),
           data_path = file.path(getwd(), 'data') )



# 3. posteriors ####
source('5_models_post.R')

model_path = file.path(getwd(), 'models_post')
model_list = list.files( model_path )
model_list = model_list[ str_detect(model_list, '.stan') ]
model_list = model_list[c(1:2,5:6,3:4)]

run_post( models = model_list,
          model_path = model_path,
          model_out = file.path(getwd(), 'chains_post'),
          data_path = file.path(getwd(), 'data') )
