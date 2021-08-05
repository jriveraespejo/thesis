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



# Chapter 4 ####



# function:
#     rmse_parameters
# description:  
#     it detects the names of the parameters in a stanfit object
# arguments:
#     stan_object = object containing a 'precis' object   
#     est_par = character vector with the names of parameters of interest
#
rmse_parameters = functon(model_path, files_int, parameters){
  
  # test
  model_path = file.path(getwd(), 'chains_post')
  parameter_path = file.path(getwd(), 'data')
  files_int = c('FOLV_CE', 'FOLV_NC', 'SOLV_CE', 'SOLV_NC')
  
  model_list = list.files( model_path )
  
  
  # indetify files of interest
  for(i in 1:length(files_int)){
    if(i==1){
      idx = which(str_detect(model_list, files_int[i]))
    } else{
      idx = c(idx, which(str_detect(model_list, files_int[i])) )
    }
  }
  model_list = model_list[idx]
  
  
  # sort files by replicate and chain
  list_model = data.frame(model_string=model_list)
  list_model$model = str_sub(model_list, start=1, end=7)
  
  start = str_locate(model_list, 'J')[,2]
  list_model$sample = as.integer(str_sub(model_list, start=start+1, end=start+3))
  
  start = str_locate(model_list, 'Ndata')[,2]
  end = str_locate(model_list, '-')[,1]
  list_model$data = as.integer(str_sub(model_list, start=start+1, end=end-1))
  
  start = str_locate(model_list, '-')[,1]
  list_model$chain = as.integer(str_sub(model_list, start=start+1, end=start+1))
  
  list_model = list_model[with(list_model, order(model, sample, data, chain)),]
  
  
  
}







fit_files = list_model[1,1]
stan_model = rstan::read_stan_csv( file.path(model_path, fit_files ) )


# true parameters

load('1_04_true_parameters.RData')

# correlation using path analysis (3x3 matrix per row)
extra1 = with(data_true,
              c(1, 
                loads$l_lit * loads$l_inf * var(theta$theta_rc),
                loads$l_lit * loads$l_ref * var(theta$theta_rc),
                loads$l_lit * loads$l_inf * var(theta$theta_rc), 
                1,
                loads$l_inf * loads$l_ref * var(theta$theta_rc),
                loads$l_lit * loads$l_ref * var(theta$theta_rc),
                loads$l_inf * loads$l_ref * var(theta$theta_rc),
                1))

# theta interleaved
theta_mom = data_true$theta[!is.na(data_true$theta$theta_rc_m),]
for(i in 1:nrow(theta_mom) ){
  if(i==1){ 
    extra2 = as.numeric( theta_mom[i, 7:9] ) 
  } else{ 
    extra2 = c(extra2, as.numeric( theta_mom[i, 7:9] ) ) 
  }
}

# vector of all parameters
true_pars = with(data_true,
                 c(texts$m_b, texts$s_b, items$b,
                   0, betas$gender, betas$age, betas$edu, betas$exp,
                   extra1, extra2 ) )

est_pars = c('m_b','s_b','b_k','a','b_gender','b_age','b_edu','b_exp',
             'Rho_theta_sub','theta_sub')


parameter_recovery = function(stan_object, est_par, true_par, prec=3){
  
  # test
  stan_object=diff
  est_par=est_pars
  true_par=true_pars
  prec=3

  # get the point estimates
  res_stan = precis(stan_object, depth=4)
  
  # identify parameters of interest
  for(j in 1:length(est_par)){
    idx_mom = str_detect( row.names(res_stan), paste0('^', est_par[j]) )
    if(j==1){
      idx = which(idx_mom)
    } else{
      idx = c( idx, which(idx_mom) )
    }
  }
  res_stan = round( res_stan[idx,], prec)
  
  # introduce true parameters
  res_stan$true = round(true_par, prec)
  
  # do the parameters have the same sign
  res_stan$same_sign = with(res_stan, as.integer(sign(true) == sign(mean)) )
  
  # identify if parameters are inside compatibility interval
  res_stan$in_CI = with(res_stan, as.integer(true>=`5.5%` & true<=`94.5%`) )
  
  
  
  # calculate rmse for samples vs true parameters
  set.seed(seed)
  post = extract.samples( stan_object )
  idx = names(post) %in% est_par
  post = post[idx]
  # str(post)
  
  rmse_sim = rep(NA, length(true_par))
  
  # j=11
  for(j in 1:length(est_par) ){ #
    
    # extract simulations of parameters
    sim_par = post[[est_par[j]]]
    dimen = dim( sim_par )
    
    if( is.na(dimen[2]) ){
      start = which(is.na(rmse_sim))[1]
      end = start
      rmse_sim[start:end] = sqrt( mean( (sim_par - true_par[start])^2 ) )
      
    } else{
      
      if( !is.na(dimen[3]) ){
        
        for(i in 1:dimen[3]){
          if(i ==1 ){
            sim_par_mom = sim_par[,,i]
          } else{
            sim_par_mom = cbind(sim_par_mom, sim_par[,,i])  
          }
        }
        sim_par = sim_par_mom
        dimen = dim( sim_par )
      }
      
      start = which(is.na(rmse_sim))[1]
      end = (start + dimen[2] - 1)
      true_rep = matrix( rep( true_par[start:end], dimen[1] ),
                         ncol=dimen[2], nrow=dimen[1], byrow=T)
      rmse_sim[start:end] = sqrt( colMeans( (sim_par - true_rep)^2 ) )
    }
  }
  res_stan$RMSE = round(rmse_sim, prec)  
  
  return(res_stan)
}
