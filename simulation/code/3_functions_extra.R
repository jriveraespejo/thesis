
# function:
#     hit_rate
# description:  
#     To calculate hit rate in the outcome, per index. It can be use to calculate 
#     hit rate per individual, or per item. It assumes the data is a long format.
# arguments:
#     idx = variable over which we do the aggregation
#     out = outcome over which we do the aggregation (same dimension as idx)
#
hit_rate <- function(idx, out){
  idx_mom = unique(idx)
  res = sapply(idx_mom, function(i){ mean( out[ idx==i ], na.rm=T ) } )
  
  return( data.frame(id=idx_mom, p=res) )
}


# function:
#     sim_hit_rate
# description:  
#     To calculate hit rate in the outcome, per index. It can be use to calculate 
#     hit rate per individual, or per item.
# arguments:
#     d_true = data list with observed information
#     d_sim = data from prior/posterior simulations
#
sim_hit_rate = function(d_true, d_sim, prob=T, hit_rate_fun=hit_rate){
  # # test
  # d_true = data_post
  # d_sim = prior_sim
  # prob=T
  # hit_rate_fun=hit_rate
  # str(prior_sim)
  
  # calculations
  S = dim(d_sim$m_b)[1]
  
  # evaluation 1
  y = array( dim=c(d_true$N, S) )
  for( i in 1:d_true$N ) {
    y[i,] = with(d_sim, theta_sub[, d_true$IDj[i], d_true$IDd[i]] - b_k[, d_true$IDk[i]] )
    y[i,] = inv_logit( y[i,] )
    if(!prob){
      y[i,] = rbinom(n=length(y[i,]), size=1, prob=y[i,] )
    }
  }
  y = data.frame(y)
  y$IDj = d_true$IDj
  y$IDk = d_true$IDk
  y$IDl = d_true$IDl
  y$IDd = d_true$IDd
  
  # storage
  res = list(
    IDind = y[, c('IDj',paste0('X', 1:S))] %>%
      group_by(IDj) %>%
      summarise_all(mean),
    IDitem = y[, c('IDk', paste0('X', 1:S))] %>%
      group_by(IDk) %>%
      summarise_all(mean),
    IDtext = y[, c('IDl', paste0('X', 1:S))] %>%
      group_by(IDl) %>%
      summarise_all(mean),
    IDdim = y[, c('IDd', paste0('X', 1:S))] %>%
      group_by(IDd) %>%
      summarise_all(mean)
  )
  
  return(res)
}


# function:
#     sim_hit_rate_reg
# description:  
#     To calculate hit rate in the outcome, based on a simple regression model. 
# arguments:
#     d_true = data list with observed information
#     d_sim = data from prior/posterior simulations
#
sim_hit_rate_reg = function(d_true, d_sim, prob=T, hit_rate_fun=hit_rate){
  
  # # test
  # d_true = data_post
  # d_sim = post_sim
  # prob=T
  # hit_rate_fun=hit_rate
  # # str(d_sim)
  
  # calculations
  S = dim(d_sim$b_gender)[1]
  
  # evaluation 1
  y = array( dim=c(d_true$N, S) )
  for( i in 1:d_true$N ) {
    y[i,] = with(d_sim, 
                 ( ind_j[ , d_true$IDj[i] ] + # random effects
                     dim_d[ , d_true$IDd[i] ] -
                     item_k[ , d_true$IDk[i] ] -
                     text_l[ , d_true$IDl[i] ] ) +
                   a + b_gender[ , d_true$gender_d[i] ] + # fixed effects
                   b_age * (d_true$age_d[i] - 29) +
                   b_edu[ , d_true$edu_d[i] ] + 
                   b_exp[ , d_true$exp_d[i] ] )
    y[i,] = inv_logit( y[i,] )
    if(!prob){
      y[i,] = rbinom(n=length(y[i,]), size=1, prob=y[i,] )
    }
  }
  y = data.frame(y)
  y$IDj = d_true$IDj
  y$IDk = d_true$IDk
  y$IDl = d_true$IDl
  y$IDd = d_true$IDd
  
  # storage
  res = list(
    IDind = y[, c('IDj',paste0('X', 1:S))] %>%
      group_by(IDj) %>%
      summarise_all(mean),
    IDitem = y[, c('IDk', paste0('X', 1:S))] %>%
      group_by(IDk) %>%
      summarise_all(mean),
    IDtext = y[, c('IDl', paste0('X', 1:S))] %>%
      group_by(IDl) %>%
      summarise_all(mean),
    IDdim = y[, c('IDd', paste0('X', 1:S))] %>%
      group_by(IDd) %>%
      summarise_all(mean)
  )
  
  return(res)
}


# function:
#     ave_hit_rate
# description:  
#     To calculate marginal hit rate in the outcome, per index. It can be use to 
#     calculate hit rate per individual, or per item.
# arguments:
#     sim_hit_object = data.frame where idx and out/mean is present
#     prob = level of confidence interval (default 0.95)
#
ave_hit_rate = function(sim_hit_object, prob=0.95){
  
  # packages
  require(rethinking)
  
  # calculations
  dimen = dim(sim_hit_object)
  mu = apply(sim_hit_object[, 2:dimen[2]], 1, mean)
  CI = data.frame( t( apply(sim_hit_object[, 2:dimen[2]], 1, PI, prob) ) )
  names(CI) = c('lower', 'upper')
  
  res = data.frame(id=sim_hit_object[,1], p=mu, CI)
  
  return(res)
}



# function:
#     ICC
# description:  
#     It creates the Item Characteristic Curve (ICC) for a dichotomous items
# arguments:
#     theta = denotes the ability of the individuals, set to a default
#     a = item discrimination parameter 
#     b = item difficulty parameter
#
ICC = function(theta=seq(-3,3, by=0.1), b){
  v = theta - b
  p = inv_logit(v)
  
  return( data.frame(theta, p) )
}


# function:
#     sim_ICC
# description:  
#     It creates simulations for ICC's based on discrimination and difficulty 
#     parameters for one item.
# arguments:
#     sim_a = simulations for item discrimination parameters (all items)
#     sim_b = simulations for item difficulty parameters (all items)
#     item_id = IDitem (default = 1)
#
sim_ICC = function(sim_b, item_id=1, ICC_fun=ICC){
  # dimension of simulation
  S = nrow(sim_b)
  
  # calculate for all samples
  for( s in 1:S ){
    res_mom = ICC_fun(b=sim_b[s, item_id])
    if(s == 1){
      res = data.frame(res_mom)
    } else{
      res = cbind(res, p=res_mom$p)
    }
  }
  
  return(res)
}  


# function:
#     ave_ICC
# description:  
#     It creates the Item Characteristic Curve (ICC) for a dichotomous items,
#     considering the average of the parameters.
# arguments:
#     theta = denotes the ability of the individuals, set to a default
#     a = item discrimination parameter 
#     b = item difficulty parameter
#
ave_ICC = function(sim_b, item_id, ICC_fun=ICC){
  mu_b = colMeans(sim_b)
  res = ICC_fun(b=mu_b[item_id])
  
  return( res )
}


# function:
#     mar_ICC
# description:  
#     It creates the Item Characteristic Curve (ICC) for a dichotomous items,
#     considering the average of the parameters.
# arguments:
#     theta = denotes the ability of the individuals, set to a default
#     a = item discrimination parameter 
#     b = item difficulty parameter
#
mar_ICC = function(sim_ICC_object){
  dimen = dim(sim_ICC_object)
  res = data.frame(theta=sim_ICC_object$theta, 
                   p=apply(sim_ICC_object[, 2:dimen[2]], 1, mean))
  
  return( res )
}


# function:
#     IIC
# description:  
#     It creates the Item Information Curve (IIC) for a dichotomous items
# arguments:
#     theta = denotes the ability of the individuals, set to a default
#     a = item discrimination parameter 
#     b = item difficulty parameter
#
IIC = function(theta=seq(-3,3, by=0.1), b, ICC_fun=ICC){
  ICC_res = ICC_fun(b=b)
  ICC_res$inf = with(ICC_res, p*(1-p) )
  
  return( ICC_res[,c('theta','inf')] )
}


# function:
#     sim_IIC
# description:  
#     It creates simulations for IIC's based on discrimination and difficulty 
#     parameters for one item.
# arguments:
#     sim_a = simulations for item discrimination parameters (all items)
#     sim_b = simulations for item difficulty parameters (all items)
#     item_id = IDitem (default = 1)
#
sim_IIC = function(sim_b, item_id=1, IIC_fun=IIC){
  
  # sim_a=prior_sim$a_k
  # sim_b=prior_sim$b_k
  # item_id=1
  # IIC_fun=IIC
  
  # dimension of simulation
  S = nrow(sim_b)
  
  # calculate for all samples
  for( s in 1:S ){
    res_mom = IIC_fun(b=sim_b[s, item_id])
    if(s == 1){
      res = data.frame(res_mom)
    } else{
      res = cbind(res, inf=res_mom$inf)
    }
  }
  
  return(res)
}  


# function:
#     ave_IIC
# description:  
#     It creates the Item Information Curve (IIC) for a dichotomous items,
#     considering the average of the parameters.
# arguments:
#     theta = denotes the ability of the individuals, set to a default
#     a = item discrimination parameter 
#     b = item difficulty parameter
#
ave_IIC = function(sim_b, item_id, IIC_fun=IIC){
  mu_b = colMeans(sim_b)
  res = IIC_fun(b=mu_b[item_id])
  
  return( res )
}


# function:
#     mar_IIC
# description:  
#     It creates the Item Information Curve (IIC) for a dichotomous items,
#     considering the average of the parameters.
# arguments:
#     theta = denotes the ability of the individuals, set to a default
#     a = item discrimination parameter 
#     b = item difficulty parameter
#
mar_IIC = function(sim_IIC_object){
  dimen = dim(sim_IIC_object)
  res = data.frame(theta=sim_IIC_object$theta, 
                   inf=apply(sim_IIC_object[, 2:dimen[2]], 1, mean))
  
  return( res )
}


# function:
#     trace_plot
# description:  
#     It creates a trace plot for a stanfit object 
# arguments:
#     stan_object = object produced after using run.jags
#     pars = character with the name of a parameter
#
trace_plot = function(stan_object, pars) {
  
  # # test
  # stan_object = stan_model
  # pars = idx[4]
  
  # packages 
  require(RColorBrewer)
  require(rethinking)
  
  # posterior
  post = rstan::extract(stan_object, pars=pars, permuted=FALSE)
  # str(post)
  
  # parameters
  n_chains = dim(post)[2]
  chain.cols = rep_len(rethink_palette, n_chains)
  wstart = 1
  wend = dim(post)[1]
  ylim = range(post[wstart:wend, , ])
  ytick = (ylim[2] - ylim[1])/6
  yaxis = round( seq(ylim[1], ylim[2], by=ytick), 2)
  neff = summary(stan_object)$summary[, "n_eff"]
  neff_use <- neff[names(neff) == pars]
  
  # plot
  plot(NULL, type="l", xlim=c(wstart, wend), ylim=ylim,
       xlab="", ylab="", axes=F)
  box(bty="l")
  axis(side=1, at=seq(0, wend, by=100))
  axis(side=2, at=yaxis, las=1 )
  # mtext(paste("n_eff =", round(neff_use, 0)), 3, adj = 1, cex = 1.1)
  mtext(pars, 3, adj = 0, cex=1.1)
  for(c in 1:n_chains){
    lines(1:wend, post[, c, ], col=chain.cols[c], lwd = 0.5)
  }
  
}



# function:
#     trank_plot
# description:  
#     It creates a trank plot for a stanfit object 
# arguments:
#     stan_object = object produced after using run.jags
#     pars = character with the name of a parameter
#     wide = controls the number of iterations (in the chain) considered.
#           (default 50)
#
trank_plot = function(stan_object, pars, wide=50){
  
  # # test
  # stan_object = stan_model
  # pars = idx[1]
  # wide=50
  
  # for colors
  require(RColorBrewer)
  require(rethinking)
  
  # posterior
  post = rstan::extract(stan_object, pars=pars, permuted=FALSE)
  # str(post)
  
  # parameters
  n_chains = dim(post)[2]
  chain.cols = rep_len(rethink_palette, n_chains)
  wstart = 1
  wend = dim(post)[1]
  neff = summary(stan_object)$summary[, "n_eff"]
  neff_use <- neff[names(neff) == pars]
  
  # rank calculation
  ranks = list()
  xrange = rep(1:wide, each=2)
  yrange = vector('list', n_chains)
  for(c in 1:n_chains){
    ranks[[c]] = rank( post[1:(wide+1), c, ] )
    y_ran = c()
    for(i in 2:(wide+1)){
      y_ran = c(y_ran, c( ranks[[c]][i-1], ranks[[c]][i] ) )
    }
    yrange[[c]] = y_ran
  }
  
  
  # plot
  plot(NULL, type='l', xlim=c(0, wide+1), ylim=c(0, wide+1),
       xlab="", ylab="", xaxt ="n", yaxt ="n", axes=F)
  box(bty="l")
  # mtext(paste("n_eff =", round(neff_use, 0)), 3, adj = 1, cex = 0.9)
  # mtext(pars, 3, adj = 0, cex = 1)
  for(c in 1:n_chains){
    lines(xrange, yrange[[c]], col=chain.cols[c], lwd=1.5)  
  }
  
}


# function:
#     acf_plot
# description:  
#     It creates a acf plot for a stanfit object 
# arguments:
#     stan_object = object produced after using run.jags
#     pars = character with the name of a parameter
#
acf_plot = function(stan_object, pars){
  
  # # test
  # stan_object = stan_model
  # pars = idx[1]
  
  # posterior
  post = rstan::extract(stan_object, pars=pars, permuted=FALSE)
  # str(post)
  
  # plot
  acf( post[, 1, ], main='', xlab='', ylab='', mar = c(0, 0, 0, 0) )
  # mtext(paste("n_eff =", round(neff_use, 0)), 3, adj = 1, cex = 0.9)
  # mtext(pars, 3, adj = 0, cex = 1)
  
}


# function:
#     tri_plot
# description:  
#     it plots trace, trank, and ACF plots for a maximul of 5 parameters
# arguments:
#     stan_object = object containing a 'precis' object   
#     pars = character VECTOR with the names of parameters of interest
#             only plots a maximum of five (5) parameters.
#

tri_plot = function(stan_object, pars){

  # # test
  # stan_object = stan_model
  # pars = paste0('m_b[',1:5,']')
  
  # ensure there is only 5 paramaters
  if(length(pars)>5){
    pars = pars[1:5]
  }
  
  # plot
  par(mfrow=c(length(pars), 3), mar=c(3,3.5,1.5,1)+0.1)
  for(i in 1:length(idx)){
    trace_plot(stan_model, pars=idx[i]) 
    trank_plot(stan_model, pars=idx[i])
    acf_plot(stan_model, pars=idx[i])
  }
  par(mfrow=c(1,1), mar=opar$mar)
  
}






# function:
#     detect_parameter
# description:  
#     it detects the names of the parameters in a stanfit object
# arguments:
#     stan_object = object containing a 'precis' object   
#     est_par = character vector with the names of parameters of interest
#
detect_parameter = function(precis_object, est_par){
  
  # identify parameters of interest
  for(j in 1:length(est_par)){
    idx_mom = str_detect( row.names(precis_object), paste0('^', est_par[j]) )
    if(j==1){
      idx = which(idx_mom)
    } else{
      idx = c( idx, which(idx_mom) )
    }
  }
  
  return(idx)
}



# function:
#     file_id
# description:  
#     it detects all stanfit object of interest in a folder
#     (it only applies for a study simulation with multiple conditions)
# arguments:
#     chains_path = location of csv files corresponding to stanfit objects   
#     model_int = character VECTOR with the names of the models of interest
#
file_id = function(chains_path, model_int){
  
  # # test
  # chains_path = file.path(getwd(), 'chains_post')
  # model_int = c('FOLV_CE', 'FOLV_NC', 'SOLV_CE', 'SOLV_NC')
  
  # list all files
  chains_list = list.files( chains_path )
  
  # identify files of interest
  for(i in 1:length(model_int)){
    if(i==1){
      idx = which(str_detect(chains_list, model_int[i]))
    } else{
      idx = c(idx, which(str_detect(chains_list, model_int[i])) )
    }
  }
  chains_list = chains_list[idx]
  
  
  # sort files by replicate and chain
  list_model = data.frame(model_string=chains_list)
  list_model$model = str_sub(chains_list, start=1, end=7)
  
  start = str_locate(chains_list, 'J')[,2]
  list_model$sample = as.integer(str_sub(chains_list, start=start+1, end=start+3))
  
  start = str_locate(chains_list, 'Ndata')[,2]
  end = str_locate(chains_list, '-')[,1]
  list_model$data = as.integer(str_sub(chains_list, start=start+1, end=end-1))
  
  start = str_locate(chains_list, '-')[,1]
  list_model$chain = as.integer(str_sub(chains_list, start=start+1, end=start+1))
  
  list_model = list_model[with(list_model, order(model, sample, data, chain)),]
  # str(list_model)
  
  return(list_model)
}


# function:
#     figures_plot
# description:  
#     it plot all specific parameters for the simulations
# arguments:
#     c_list = data.frame generated with file_id() function
#     chains_path = location of csv files corresponding to stanfit objects
#     file_save = location to save images
#
figures_plot = function(c_list, chains_path, file_save){
  
  # # test
  # c_list=chains_list
  # chains_path=chains_path
  # file_save=file.path(getwd(), 'figures2')
  
  # parameters
  models_int = unique(c_list$model)
  sample_sizes = unique(c_list$sample)
  data_number = unique(c_list$data)
  
  # m=1
  # s=1
  # d=1
  # figures
  for(m in 1:length(models_int)){
    for(s in 1:length(sample_sizes)){
      for(d in 1:length(data_number)){
        
        # identify files
        idx_files = with(chains_list, model==models_int[m] & 
                           sample==sample_sizes[s] &
                           data == data_number[d])
        fit_files = chains_list[idx_files, 1]
        
        # load stan data
        stan_model = rstan::read_stan_csv( file.path(chains_path, fit_files ) )
        
        
        # figures
        file_name = paste0(models_int[m], '_J', sample_sizes[s], '_Ndata', data_number[d],
                           '_', 'mb', '.png')
        file_dir = file.path(file_save, file_name)
        png(file_dir, units='cm', width=30, height=30, res=100)
        idx = paste0('m_b[',1:5,']')
        tri_plot(stan_model, pars=idx)
        dev.off()
        
        
        file_name = paste0(models_int[m], '_J', sample_sizes[s], '_Ndata', data_number[d],
                           '_', 'sb', '.png')
        file_dir = file.path(file_save, file_name)
        png(file_dir, units='cm', width=30, height=30, res=100)
        idx = paste0('s_b[',1:5,']')
        tri_plot(stan_model, pars=idx)
        dev.off()
        
        
        file_name = paste0(models_int[m], '_J', sample_sizes[s], '_Ndata', data_number[d],
                           '_', 'bk1', '.png')
        file_dir = file.path(file_save, file_name)
        png(file_dir, units='cm', width=30, height=30, res=100)
        idx = paste0('b_k[',1:4,']')
        tri_plot(stan_model, pars=idx)
        dev.off()
        
        
        file_name = paste0(models_int[m], '_J', sample_sizes[s], '_Ndata', data_number[d],
                           '_', 'bk2', '.png')
        file_dir = file.path(file_save, file_name)
        png(file_dir, units='cm', width=30, height=30, res=100)
        idx = paste0('b_k[',5:8,']')
        tri_plot(stan_model, pars=idx)
        dev.off()
        
        
        file_name = paste0(models_int[m], '_J', sample_sizes[s], '_Ndata', data_number[d],
                           '_', 'bk3', '.png')
        file_dir = file.path(file_save, file_name)
        png(file_dir, units='cm', width=30, height=30, res=100)
        idx = paste0('b_k[',9:12,']')
        tri_plot(stan_model, pars=idx)
        dev.off()
        
        
        file_name = paste0(models_int[m], '_J', sample_sizes[s], '_Ndata', data_number[d],
                           '_', 'reg1', '.png')
        file_dir = file.path(file_save, file_name)
        png(file_dir, units='cm', width=30, height=30, res=100)
        idx = c('a', paste0('b_G[',1:2,']'), 'b_A')
        tri_plot(stan_model, pars=idx)
        dev.off()
        
        
        file_name = paste0(models_int[m], '_J', sample_sizes[s], '_Ndata', data_number[d],
                           '_', 'reg2', '.png')
        file_dir = file.path(file_save, file_name)
        png(file_dir, units='cm', width=30, height=30, res=100)
        idx = c( paste0('b_E[',1:3,']') )
        tri_plot(stan_model, pars=idx)
        dev.off()
        
        
        file_name = paste0(models_int[m], '_J', sample_sizes[s], '_Ndata', data_number[d],
                           '_', 'reg3', '.png')
        file_dir = file.path(file_save, file_name)
        png(file_dir, units='cm', width=30, height=30, res=100)
        idx = c( paste0('b_X[',1:4,']') )
        tri_plot(stan_model, pars=idx)
        dev.off()
        
        
        file_name = paste0(models_int[m], '_J', sample_sizes[s], '_Ndata', data_number[d],
                           '_', 'Rho', '.png')
        file_dir = file.path(file_save, file_name)
        png(file_dir, units='cm', width=30, height=30, res=100)
        idx = c('Rho_theta_sub[1,2]', 'Rho_theta_sub[1,3]', 'Rho_theta_sub[2,3]')
        tri_plot(stan_model, pars=idx)
        dev.off()
        
        
        # selecting sample of individuals
        idx_num = round( seq(1, sample_sizes[s], length.out=5), 0)
        
        
        file_name = paste0(models_int[m], '_J', sample_sizes[s], '_Ndata', data_number[d],
                           '_', 'theta_sub1', '.png')
        file_dir = file.path(file_save, file_name)
        png(file_dir, units='cm', width=30, height=30, res=100)
        idx = paste0('theta_sub[', idx_num, ',1]')
        tri_plot(stan_model, pars=idx)
        dev.off()
        
        
        file_name = paste0(models_int[m], '_J', sample_sizes[s], '_Ndata', data_number[d],
                           '_', 'theta_sub2', '.png')
        file_dir = file.path(file_save, file_name)
        png(file_dir, units='cm', width=30, height=30, res=100)
        idx = paste0('theta_sub[', idx_num, ',2]')
        tri_plot(stan_model, pars=idx)
        dev.off()
        
        
        file_name = paste0(models_int[m], '_J', sample_sizes[s], '_Ndata', data_number[d],
                           '_', 'theta_sub3', '.png')
        file_dir = file.path(file_save, file_name)
        png(file_dir, units='cm', width=30, height=30, res=100)
        idx = paste0('theta_sub[', idx_num, ',3]')
        tri_plot(stan_model, pars=idx)
        dev.off()
        
        
        if( str_detect(models_int[m], 'SOLV') ){
          
          file_name = paste0(models_int[m], '_J', sample_sizes[s], '_Ndata', data_number[d],
                             '_', 'loads', '.png')
          file_dir = file.path(file_save, file_name)
          png(file_dir, units='cm', width=30, height=30, res=100)
          idx = paste0('loads[', 1:3, ']')
          tri_plot(stan_model, pars=idx)
          dev.off()
          
          
          file_name = paste0(models_int[m], '_J', sample_sizes[s], '_Ndata', data_number[d],
                             '_', 'theta', '.png')
          file_dir = file.path(file_save, file_name)
          png(file_dir, units='cm', width=30, height=30, res=100)
          idx = paste0('theta[', idx_num, ']')
          tri_plot(stan_model, pars=idx)
          dev.off()
          
        }
        
        print( paste0(models_int[m], '_J', sample_sizes[s], '_Ndata', data_number[d]) )
      }
    }
  }
  
}


# function:
#     parameter_recovery
# description:  
#     It display the parameter recovery compared to a set of true parameters.
#     Among the calculations are: mean, sd, compatibility interval, effective
#     sample size Rhat4 for estimated parameters, and comparison against the
#     true parameters like: parameters with same sign (sam_sign), if the 
#     parameter is contained in the compatibility interval (in_CI), and the
#     rmse calculated with a representative sample from the posterior (RMSE).
# arguments:
#     stan_object = object containing a stanfit object (it can be a list also)
#     est_par = character vector with the names of parameters of interest
#     true_par = vector of values for true parameters
#
parameter_recovery = function(stan_object, est_par, true_par, 
                              diff=F, prec=3, seed=1){
  
  # # test
  # stan_object=diff
  # est_par=diff_pars 
  # true_par=diff_true
  # prec=3
  # seed=1
  
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



# function:
#     diff_recovery
# description:  
#     It display the parameter differences recovery compared to a set of true 
#     differences.
# arguments:
#     stan_object = object containing a stanfit object (it can be a list also)
#     est_diff = character vector with the names of parameters of interest
#     true_diff = vector of values for true differences
#
diff_recovery = function(stan_object, est_diff, true_diff, prec=3, seed=1){
  
  # # test
  # stan_object=stan_model_nc
  # est_diff=diff_pars
  # true_diff=diff_true
  # prec=3
  # seed=1
  
  # extract samples
  set.seed(seed)
  post = extract.samples( stan_object )
  idx = names(post) %in% est_diff
  post = post[idx]
  # names(post)
  
  # calculating results
  res_stan = precis(stan_object, depth=4)
  
  # calculations
  # k=1
  for(k in 1:length(est_diff)){
    
    # selecting parameter
    idx = detect_parameter(res_stan, est_diff[k])
    lab_par = rownames(res_stan)[idx]
    npars = length(idx)
    
    # storage
    # i=2
    # j=3
    for(i in 1:npars){
      for(j in 1:npars){
        
        if(j>i){
          if(i == 1 & j==2 & k==1){
            diff = post[[est_diff[k]]][,j] - post[[est_diff[k]]][,i] 
            diff_name = paste( c( lab_par[j], lab_par[i]), collapse=' - ' )
          } else{
            diff = cbind(diff ,
                         post[[est_diff[k]]][,j] - post[[est_diff[k]]][,i] ) 
            diff_name = c(diff_name, 
                          paste( c( lab_par[j], lab_par[i]), collapse=' - ' ) )
          }
        }
        
      }
    }
    
  }
  
  # storage
  diff = as_tibble(diff)
  names(diff) = diff_name
  res_diff = precis( diff, depth=4 )
  res_diff[,1:4] = round( res_diff[,1:4], prec)
  
  # introduce true differences
  res_diff$true = round( true_diff, prec)
  
  # do the parameters have the same sign
  res_diff$same_sign = with(res_diff, as.integer(sign(true) == sign(mean)) )
  
  # identify if parameters are inside compatibility interval
  res_diff$in_CI = with(res_diff, as.integer(true>=`5.5%` & true<=`94.5%`) )
  
  # rmse
  dimen = dim(diff)
  true_rep = matrix( rep( diff_true, dimen[1] ),
                     ncol=dimen[2], nrow=dimen[1], byrow=T)
  rmse_sim = sqrt( colMeans( (diff - true_rep)^2 ) )
  res_diff$RMSE = round(rmse_sim, prec)  
  
  return(res_diff)
}




# function:
#     n_effective
# description:  
#     It joins the number of effective samples from two stan models.
#     normally is used to compared the centered vs non-centered version.
# arguments:
#     stan_object1 = object containing a stanfit object
#     stan_object2 = object containing a stanfit object
#     est_par = character vector with the names of parameters of interest
#
n_effective = function(stan_object1, stan_object2, est_par){
  
  # stan model 1 (centered)
  precis1 = precis(stan_object1, depth=4)
  idx1 = detect_parameter(precis1, est_par)
  
  # stan model 2 (non-centered)
  precis2 = precis(stan_object2, depth=4)
  idx2 = detect_parameter(precis2, est_par)
  
  # plot
  neff_table = data.frame(
    centered=precis1$n_eff[idx1],
    non_centered=precis2$n_eff[idx2]
    )
  
  return(neff_table)
}