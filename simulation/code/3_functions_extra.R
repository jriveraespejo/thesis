
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
  # d_sim = post_sim
  # prob=T
  # hit_rate_fun=hit_rate
  # # str(d_sim)
  
  # calculations
  S = dim(d_sim$m_b)[1]
  N = d_true$N
  
  # evaluation 1
  y = array( dim=c(N, S) )
  for( i in 1:N ) {
    y[i,] = with(d_sim, theta_sub[, d_true$IDj[i], d_true$IDd[i]] - b_k[, d_true$IDk[i]] )
    y[i,] = inv_logit( y[i,] )
    if(!prob){
      y[i,] = rbinom(n=length(y[i,]), size=1, prob=y[i,] )
    }
  }
  var_true = c('IDj','IDk','IDl','IDd','GE','AG','ED','XP')
  y = cbind(data.frame(y), data.frame(d_true[var_true]) )
  
  
  ####
  # calculations
  ####
  
  # individuals
  IDmom = y[, c('IDj',paste0('X', 1:S))] %>%
    group_by(IDj) %>%
    summarise_all(mean)
  IDmom_true = data.frame(d_true[c('IDj','y')]) %>%
    group_by(IDj) %>%
    summarise(true=mean(y))
  IDmom = tibble(merge(IDmom_true, IDmom, by='IDj'))
  IDmom = IDmom[with(IDmom, order(IDj)),]
  mom = IDmom[,3:ncol(IDmom)] - matrix(rep(IDmom$true, S), ncol=S, byrow=F)
  mom = mom^2
  IDmom$RMSE = sqrt(rowMeans(mom))
  IDmom$Sdiff = (rowMeans(IDmom[,3:ncol(IDmom)]) - IDmom$true)^2
  IDmom = IDmom[,c('IDj','true','RMSE','Sdiff',paste0('X',1:S))]
  res = list( IDind = IDmom)
  
  
  # individuals per covariates
  d_mom = data.frame(d_true[c('IDind','G','A','E','X')])
  res$IDind = tibble(merge(res$IDind, d_mom, by.x='IDj', by.y='IDind'))
  res$IDind = res$IDind[, c( 1, (S+5):(S+8), 2:4, 5:(S+4) )]
  
  
  # items
  IDmom = y[, c('IDk',paste0('X', 1:S))] %>%
    group_by(IDk) %>%
    summarise_all(mean)
  IDmom_true = data.frame(d_true[c('IDk','y')]) %>%
    group_by(IDk) %>%
    summarise(true=mean(y))
  IDmom = tibble(merge(IDmom_true, IDmom, by='IDk'))
  IDmom = IDmom[with(IDmom, order(IDk)),]
  mom = IDmom[,3:ncol(IDmom)] - matrix(rep(IDmom$true, S), ncol=S, byrow=F)
  mom = mom^2
  IDmom$RMSE = sqrt(rowMeans(mom))
  IDmom$Sdiff = (rowMeans(IDmom[,3:ncol(IDmom)]) - IDmom$true)^2
  IDmom = IDmom[,c('IDk','true','RMSE','Sdiff',paste0('X',1:S))]
  res$IDitems = IDmom
  
  
  # texts
  IDmom = y[, c('IDl',paste0('X', 1:S))] %>%
    group_by(IDl) %>%
    summarise_all(mean)
  IDmom_true = data.frame(d_true[c('IDl','y')]) %>%
    group_by(IDl) %>%
    summarise(true=mean(y))
  IDmom = tibble(merge(IDmom_true, IDmom, by='IDl'))
  IDmom = IDmom[with(IDmom, order(IDl)),]
  mom = IDmom[,3:ncol(IDmom)] - matrix(rep(IDmom$true, S), ncol=S, byrow=F)
  mom = mom^2
  IDmom$RMSE = sqrt(rowMeans(mom))
  IDmom$Sdiff = (rowMeans(IDmom[,3:ncol(IDmom)]) - IDmom$true)^2
  IDmom = IDmom[,c('IDl','true','RMSE','Sdiff',paste0('X',1:S))]
  res$IDtext = IDmom
  
  
  # dimensions
  IDmom = y[, c('IDd',paste0('X', 1:S))] %>%
    group_by(IDd) %>%
    summarise_all(mean)
  IDmom_true = data.frame(d_true[c('IDd','y')]) %>%
    group_by(IDd) %>%
    summarise(true=mean(y))
  IDmom = tibble(merge(IDmom_true, IDmom, by='IDd'))
  IDmom = IDmom[with(IDmom, order(IDd)),]
  mom = IDmom[,3:ncol(IDmom)] - matrix(rep(IDmom$true, S), ncol=S, byrow=F)
  mom = mom^2
  IDmom$RMSE = sqrt(rowMeans(mom))
  IDmom$Sdiff = (rowMeans(IDmom[,3:ncol(IDmom)]) - IDmom$true)^2
  IDmom = IDmom[,c('IDd','true','RMSE','Sdiff',paste0('X',1:S))]
  res$IDdim = IDmom

  
  # return
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
#     It creates the Item Characteristic Curve (ICC) for dichotomous items
# arguments:
#     theta = denotes the ability of the individuals, set to a default
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
#     It creates simulations for ICC's based on difficulty parameters for 
#     one item.
# arguments:
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
#     It creates the Item Characteristic Curve (ICC) for dichotomous items,
#     considering the average of the parameters.
# arguments:
#     theta = denotes the ability of the individuals, set to a default
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
#     It creates the Item Characteristic Curve (ICC) for dichotomous items,
#     considering the average of the parameters.
# arguments:
#     sim_ICC_object = object produced by sim_ICC function
#
mar_ICC = function(sim_ICC_object){
  dimen = dim(sim_ICC_object)
  res = data.frame(theta=sim_ICC_object$theta, 
                   p=apply(sim_ICC_object[, 2:dimen[2]], 1, mean))
  
  return( res )
}


# function:
#     IIF
# description:  
#     It creates the Item Information Function (IIF) for a dichotomous items
# arguments:
#     theta = denotes the ability of the individuals, set to a default
#     b = item difficulty parameter
#
IIF = function(theta=seq(-3,3, by=0.1), b, ICC_fun=ICC){
  ICC_res = ICC_fun(b=b)
  ICC_res$inf = with(ICC_res, p*(1-p) )
  
  return( ICC_res[,c('theta','inf')] )
}


# function:
#     sim_IIF
# description:  
#     It creates simulations for IIF's based on difficulty parameters for 
#     one item.
# arguments:
#     sim_b = simulations for item difficulty parameters (all items)
#     item_id = IDitem (default = 1)
#
sim_IIF = function(sim_b, item_id=1, IIF_fun=IIF){
  
  # sim_a=prior_sim$a_k
  # sim_b=prior_sim$b_k
  # item_id=1
  # IIF_fun=IIF
  
  # dimension of simulation
  S = nrow(sim_b)
  
  # calculate for all samples
  for( s in 1:S ){
    res_mom = IIF_fun(b=sim_b[s, item_id])
    if(s == 1){
      res = data.frame(res_mom)
    } else{
      res = cbind(res, inf=res_mom$inf)
    }
  }
  
  return(res)
}  


# function:
#     ave_IIF
# description:  
#     It creates the Item Information Function (IIF) for dichotomous items,
#     considering the average of the parameters.
# arguments:
#     theta = denotes the ability of the individuals, set to a default
#     b = item difficulty parameter
#
ave_IIF = function(sim_b, item_id, IIF_fun=IIF){
  mu_b = colMeans(sim_b)
  res = IIF_fun(b=mu_b[item_id])
  
  return( res )
}


# function:
#     mar_IIF
# description:  
#     It creates the Item Information Function (IIF) for a dichotomous items,
#     considering the average of the parameters.
# arguments:
#     sim_IIF_object = object produced by sim_IIF function
#
mar_IIF = function(sim_IIF_object){
  dimen = dim(sim_IIF_object)
  res = data.frame(theta=sim_IIF_object$theta, 
                   inf=apply(sim_IIF_object[, 2:dimen[2]], 1, mean))
  
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
    trace_plot(stan_object, pars=idx[i]) 
    trank_plot(stan_object, pars=idx[i])
    acf_plot(stan_object, pars=idx[i])
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
  # model_int = c('FOLV_CE_mod','FOLV_CE','FOLV_NC_mod','FOLV_NC','SOLV_CE','SOLV_NC')
  # 
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
  
  end = str_locate(chains_list, 'J')[,1]
  list_model$model = str_sub(chains_list, start=1, end=end-2)
  
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
  
  # test
  c_list=chains_list
  chains_path=chains_path
  file_save=file.path(getwd(), 'figures3')
  
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


# function:
#     stat_chain
# description:  
#     it extraxt the information of Rhat and n_eff 
# arguments:
#     c_list = data.frame generated with file_id() function
#     chains_path = location of csv files corresponding to stanfit objects
#     file_save = location to save files
#     file_name = name of the file to save 
#     contr_pars = set of paramerter to calculate contrasts
#
stat_chain = function(c_list, chains_path, file_save, file_name, contr_pars){
  
  # # test
  # c_list=chains_list
  # chains_path=chains_path
  # file_save=file.path(getwd(), 'figures4')
  # file_name = 'stan_stats_no_mod'
  # contr_pars = c('b_G','b_E','b_X')
  
  # parameters
  models_int = unique(c_list$model)
  sample_sizes = unique(c_list$sample)
  data_number = unique(c_list$data)
  
  # m=1
  # s=1
  # d=1
  for(m in 1:length(models_int)){
    for(s in 1:length(sample_sizes)){
      for(d in 1:length(data_number)){
        
        ####
        # simple parameter section
        ####
        
        # identify files of interest
        idx_files = with(chains_list, 
                         model==models_int[m] &
                           sample==sample_sizes[s] &
                           data == data_number[d])
        
        fit_files = chains_list[idx_files, 1]
        stan_model = rstan::read_stan_csv( file.path(chains_path, fit_files) )
        
        # calculate parameters
        stan_result = precis(stan_model, depth=4)
        
        # storage
        if( all(m==1, s==1, d==1) ){
          stan_int = data.frame( model_type = models_int[m],
                                 sample_size = sample_sizes[s],
                                 data_number = data_number[d],
                                 parameter = row.names(stan_result),
                                 stan_result )
        } else {
          stan_int = rbind(stan_int,
                           data.frame( model_type = models_int[m],
                                       sample_size = sample_sizes[s],
                                       data_number = data_number[d],
                                       parameter = row.names(stan_result),
                                       stan_result ) )
        }
        
        ####
        # contrasts section
        ####
        
        # extract samples
        post = extract.samples( stan_model )
        idx = names(post) %in% contr_pars
        
        if( !all(idx==F) ){
          
          post = post[idx]
          # names(post)
          
          # calculations
          # k=1
          for(k in 1:length(contr_pars) ){
            
            # selecting parameter
            idx = detect_parameter(stan_result, contr_pars[k])
            lab_par = rownames(stan_result)[idx]
            npars = length(idx)
            
            # storage
            # i=2
            # j=3
            for(i in 1:npars){
              for(j in 1:npars){
                
                if(j>i){
                  if(i == 1 & j==2 & k==1){
                    diff = post[[contr_pars[k]]][,j] - post[[contr_pars[k]]][,i] 
                    diff_name = paste( c( lab_par[j], lab_par[i]), collapse=' - ' )
                  } else{
                    diff = cbind(diff ,
                                 post[[contr_pars[k]]][,j] - post[[contr_pars[k]]][,i] ) 
                    diff_name = c(diff_name, 
                                  paste( c( lab_par[j], lab_par[i]), collapse=' - ' ) )
                  }
                }
                
              }
            }
            
          }
          
          # calculations
          diff = data.frame(diff)
          names(diff) = diff_name
          res_diff = precis( diff, depth=4 )
          
          # contrast storage
          res_diff = data.frame( model_type = models_int[m],
                                 sample_size = sample_sizes[s],
                                 data_number = data_number[d],
                                 parameter = row.names(res_diff),
                                 res_diff[,-5], 
                                 n_eff = NA,
                                 Rhat4 = NA)
          stan_int = rbind(stan_int, res_diff)
        }
        
        # remove row.names
        row.names(stan_int) = NULL
        
        # remove unnecessary parameters
        uu = c('zb_k','L_Rho_theta_sub','ztheta_sub','ztheta','m_mult','m_theta')
        for(i in 1:length(uu)){
          idx_mom = str_detect(stan_int$parameter, uu[i])
          if(i==1){
            idx = idx_mom
          } else{
            idx = idx | idx_mom
          }
        }
        stan_int = stan_int[!idx,]
        
        # save file
        save(stan_int, file=file.path(file_save, paste0(file_name, '.RData') ) )
        print( paste0(models_int[m], '_J', sample_sizes[s], '_Ndata', data_number[d]) )
      
      }
    }
  }
  
}



# function:
#     plot_stat
# description:  
#     It plots the statistics of interest 
# arguments:
#     stat_object = object produced by the function stat_chain()
#     info = statistics of interest
#     par_int = parameters of interest
#     model = options according to model in stat_chain
#     ssize = sample size of simulation
#     title = main title for the plot 
#
plot_stat = function(stat_object, info, par_int, model, ssize=100, title='(A)'){
  
  # # test
  # stat_object = stan_int
  # info = 'n_eff'
  # par_int = c( paste0('b_G[',1:2,']'),'b_A', paste0('b_E[',1:3,']'),
  #              paste0('b_X[',1:4,']'))
  # model = 'FOLV'
  # ssize = 100
  # title='(A)'

  
  # identify the parameters
  # i=1
  for(i in 1:length(par_int)){
    idx_pars_mom = stat_object$parameter == par_int[i]
    if(i==1){
      idx_pars = idx_pars_mom
    } else{
      idx_pars = idx_pars | idx_pars_mom
    }
  }
  # sum(idx_pars)
  
  # identify models
  model_int = unique(with(stat_object, model_type[str_detect(model_type, model)]))
  for(i in 1:length(model_int)){
    idx = with(stat_object, model_type==model_int[i] & sample_size==ssize & idx_pars)
    stat_mom = stat_object[idx, c('sample_size','data_number','parameter', info)]
    
    if(i==1){
      stat_final = stat_mom
    } else{
      stat_final = merge(stat_final, stat_mom, 
                         by=c('sample_size','data_number','parameter'))
    }
  }
  
  # parameter color
  par_mod = str_locate(stat_final$parameter, '[:digit:]')[,1] - 2
  par_mod = ifelse( is.na(par_mod), 3, par_mod)
  par_mod = str_sub(stat_final$parameter, start = 1, end=par_mod)
  par_mod_un = unique(par_mod)
  col_pars = rep(NA, nrow(stat_final))
  for(i in 1:length(par_mod_un)){
    col_pars = ifelse( par_mod==par_mod_un[i], 
                       col.alpha(i, 0.4), col_pars) 
  }
  
  # plot
  idx_var = str_detect( names(stat_final), info)  
  x_lim = range(stat_final[, idx_var])
  plot(stat_final[,idx_var], xlim=x_lim, ylim=x_lim, main=title,
       col=col_pars, pch=19,
       xlab='Centered parametrization', ylab='Non-centered parametrization')
  abline(a=0, b=1, lty=2)
  
  if(info=='Rhat4'){
    abline(v=1.05, h=1.05, lty=3, col=col.alpha('black', 0.6))
    legend('top', horiz=T, par_mod_un, pch=19, col=unique(col_pars), bty='n')
  } else{
    legend('bottom', horiz=T, par_mod_un, pch=19, col=unique(col_pars), bty='n')
  }
  
}




# function:
#     extract_true
# description:  
#     Extracts the true parameters in simulations 
# arguments:
#     file_load = location of true parameters per replica
#     file_save = location to save file
#
extract_true = function(file_load, file_save){
  
  # # test
  # file_load = file.path(getwd(), 'data')
  # file_save = file.path(getwd(), 'figures5')
  
  # iteration set
  ss = unique(stan_int$sample_size)
  dn = unique(stan_int$data_number)
  
  # s=1
  # d=1
  for(s in 1:length(ss)){
    for(d in 1:length(dn)){
      
      # load parameters
      file_name = paste0('Parameters_J', ss[s], '_l0.95_Ndata', dn[d], '.RData')
      load( file.path(file_load, file_name) )
      
      # extract parameter of interest
      # texts
      par_names = expand_grid( par = names(data_true$texts)[2:ncol(data_true$texts)], 
                               number = data_true$texts$IDtext)
      par_names = paste0(par_names$par,'[', par_names$number, ']')
      true_int = data.frame( parameter = par_names,
                             true = c(data_true$texts$m_b, data_true$texts$s_b) )
      
      # items
      par_names = paste0('b_k[', data_true$items$IDitem, ']')
      true_int = rbind(true_int, data.frame( parameter = par_names,
                                             true = data_true$items$b ) )
      
      # regression
      par_names = c('a', paste0('b_G[',1:2,']'), 'b_A', paste0('b_E[',1:3,']'),
                    paste0('b_X[',1:4,']'), paste0('loads[',1:3,']'),
                    'Rho_theta_sub[1,2]', 'Rho_theta_sub[1,3]', 'Rho_theta_sub[2,3]',
                    'Rho_theta_sub[2,1]', 'Rho_theta_sub[3,1]', 'Rho_theta_sub[3,2]',
                    'Rho_theta_sub[1,1]', 'Rho_theta_sub[2,2]', 'Rho_theta_sub[3,3]')
      true_mom = c(0, unlist(data_true$betas), data_true$betas$exp_corr, c(1,1,1) )
      true_int = rbind(true_int, data.frame( parameter = par_names,
                                             true = true_mom ) )
      # NOTICE: Rho_theta_sub corresponds only to FOLV
      
      # contrasts
      par_names = c("b_G[2] - b_G[1]","b_E[2] - b_E[1]","b_E[3] - b_E[1]",
                    "b_E[3] - b_E[2]","b_X[2] - b_X[1]","b_X[3] - b_X[1]",
                    "b_X[4] - b_X[1]","b_X[3] - b_X[2]","b_X[4] - b_X[2]",
                    "b_X[4] - b_X[3]")
      true_mom = with(data_true$betas, 
                      c( gender[2] - gender[1], 
                         edu[2] - edu[1], edu[3] - edu[1], edu[3] - edu[2],
                         exp[2] - exp[1], exp[3] - exp[1], exp[4] - exp[1],
                         exp[3] - exp[2], exp[4] - exp[2], exp[4] - exp[3]) )
      true_int = rbind(true_int, data.frame( parameter = par_names,
                                             true = true_mom ) )
      
      # abilities
      par_names = c(paste0('theta[',1:data_true$J,']'), 
                    paste0('theta_sub[', 1:data_true$J,',1]'),
                    paste0('theta_sub[', 1:data_true$J,',2]'),
                    paste0('theta_sub[', 1:data_true$J,',3]'))
      true_mom = with(data_true$abilities, c( theta, theta1, theta2, theta3) )
      true_int = rbind(true_int, data.frame( parameter = par_names,
                                             true = true_mom ) )
      
      
      # storage
      rownames(true_int) = NULL
      if(s==1 & d==1){
        true_pars = data.frame(sample_size = ss[s],
                               data_number = dn[d],
                               true_int)
      } else{
        true_pars = rbind(true_pars, 
                          data.frame(sample_size = ss[s],
                                     data_number = dn[d],
                                     true_int) )
      }
      
      # save file
      save(true_pars, file=file.path(file_save, 'true_pars.RData') )
      print( paste0('Parameters_J', ss[s], '_Ndata', dn[d]) )
      
    }
  }
  
}



# function:
#     rmse_pars
# description:  
#     joints statistics and true parameters, and finally calculates the RMSE 
# arguments:
#     stats_object = object generated with stat_chain() function
#     true_object = object generated with extract_true() function
#     file_save = location to save all generated files
#
rmse_pars = function(c_list, chains_path, true_load, file_save){
  
  # # test
  # stats_object = stan_int
  # true_object = true_pars
  # file_save = file.path(getwd(), 'figures5')
  
  # join data
  res_stan = merge(stats_object, true_object, all.x=T,
                  by=c("sample_size","data_number","parameter"))
  
  # change Rho in SOLV
  idx = with(res_stan, 
             str_detect(parameter, '^Rho') & 
               str_detect(model_type, '^SOLV') )
  res_stan$true[idx] = with(res_stan[idx,], ifelse(true!=1, 0, true))
  save(res_stan, file=file.path(file_save, 'results_pars.RData') )
  # simulation contemplated zero conditional correlation under SOLV, however,
  # true parameters generated with extract_true() function only contemplates
  # the implied covariance under FOLV (we need to correct that)
  
  # calculate RMSE
  res_stan$Sdiff = with(res_stan, (mean - true)^2 )
  
  stan_rmse = res_stan %>%
    group_by(model_type, sample_size, parameter) %>%
    summarise( RMSE = sqrt(mean(Sdiff)), n=n() )
  save(stan_rmse, file=file.path(file_save, 'rmse_pars.RData') )
  
  
  # ### NOT IMPLEMENTED
  # 
  # # # test
  # # c_list = chains_list
  # # chains_path = file.path(getwd(), 'chains_post')
  # # true_load = file.path(getwd(), 'figures5')
  # # file_save = file.path(getwd(), 'figures5')
  # 
  # 
  # # parameters
  # mt = unique(c_list$model)
  # ss = unique(c_list$sample)
  # dn = unique(c_list$data)
  # 
  # # plots
  # # m=1
  # # s=1
  # # d=1
  # for(m in 1:length(mt)){
  #   for(s in 1:length(ss)){
  #     for(d in 1:length(dn)){
  #       
  #       ###
  #       # 1. loading required parameters
  #       ###
  #       
  #       # loading true parameters
  #       load( file.path(true_load, 'true_pars.RData') )
  #       
  #       # identify files of stan
  #       idx_files = with(c_list, model==mt[m] & sample==ss[s] & data==dn[d])
  #       fit_files = c_list[idx_files, 1]
  #       stan_model = rstan::read_stan_csv( file.path(chains_path, fit_files) )
  #       
  #       # calculate parameters
  #       stan_mom = precis(stan_model, depth=4)
  #       stan_result = data.frame( model_type=mt[m],
  #                                 sample_size=ss[s],
  #                                 data_number=dn[d],
  #                                 parameter=row.names(stan_mom),
  #                                 stan_mom)
  #       row.names(stan_result) = NULL
  #       
  #       # merge with true
  #       idx = with(true_pars, sample_size==ss[s] & data_number==dn[d] )
  #       stan_result = merge(stan_result, true_pars[idx, ], all.x=T,
  #                           by=c("sample_size","data_number","parameter"))
  #       stan_result$true
  #       
  #       # extracting samples
  #       set.seed(7856+m+s+d)
  #       post_sim = extract.samples( stan_model )
  #       # str(post_sim)
  #     }
  #   }
  # }
  
}




# function:
#     detect_par
# description:  
#     it detects the names of the parameters in an results object generated 
#     with rmse_pars() object
# arguments:
#     result_object = result object generated with rmse_pars() object   
#     est_par = character vector with the names of parameters of interest
#
detect_par = function(result_object, est_par){
  
  # identify parameters of interest
  for(j in 1:length(est_par)){
    idx_mom = str_detect( result_object$parameter, paste0('^', est_par[j]) )
    if(j==1){
      idx = idx_mom
    } else{
      idx = idx | idx_mom
    }
  }
  
  return(idx)
}




# function:
#     recovery_plots
# description:  
#     It plots all relevant recovery plots 
# arguments:
#     result object = object generated with rmse_pars() function
#     file_save = location to save file
#
recovery_plots = function(result_object, figure_save){
  
  # # test
  # result_object = res_stan
  # figure_save = file_save
  
  # figure parameters
  opar = par()
  
  # indexes
  mt = unique(result_object$model_type)
  ss = unique(result_object$sample_size)
  dn = unique(result_object$data_number)
  # names(result_object)
  
  # plots
  # m=1
  # s=1
  # d=1
  for(m in 1:length(mt)){
    for(s in 1:length(ss)){
      for(d in 1:length(dn)){

        ### a. item parameters
        figure_name = paste0(mt[m],'_J',ss[s],'_Ndata',dn[d],'_items.png')
        png( file.path(figure_save, figure_name), 
             units='cm', width=30, height=12, res=200)
        
        par(mfrow=c(1,2))
        
        # identify parameter
        idx = with(res_stan, model_type==mt[m] & sample_size==ss[s] & 
                     data_number==dn[d] & str_detect(parameter, 'b_k'))
        res_mom = res_stan[idx,]
        
        # sorting
        loc_name = str_locate(res_mom$parameter, '[:punct:][:digit:]{1,2}')
        num_name = as.integer(str_sub(res_mom$parameter, start=loc_name[,1]+1, end=loc_name[,2]))
        idx_ord = order(num_name)
        res_mom = res_mom[idx_ord,]
        
        y_lim = range( with(res_mom, c(X5.5., X94.5., true)), na.rm=T )
        
        
        # plot 1: items difficulties
        plot(1:nrow(res_mom), res_mom$mean, ylim=y_lim, xaxt='n', yaxt='n',
             col=col.alpha('blue', 0.3), pch=19, 
             xlab='', ylab='estimates', main='Item difficulty')
        abline(v=c(5,10,15,20)+0.5, col=col.alpha('black', 0.5), lty=2)
        abline(h=0, col=col.alpha('black', 0.3), lty=2)
        axis(side=1, at=1:nrow(res_mom), labels=res_mom$parameter, las=2)
        axis(side=2, at=round( seq(y_lim[1], y_lim[2], by=0.5), 2), las=1)
        for(i in 1:nrow(res_mom)){
          lines(x=rep(i,2), y=with( res_mom[i,], c(X5.5., X94.5.) ) ,
                col=col.alpha('blue', 0.3))
        }
        points(1:nrow(res_mom), res_mom$true, col=col.alpha('black', 0.5))
        legend('topleft', legend=c('estimate','true'), pch=c(19,1), cex=.8, bty='n',
               col=c(col.alpha('blue', 0.3), col.alpha('black', 0.5)) )
        
        
        # identify parameter
        idx_mom = detect_par(res_stan, c('m_b','s_b'))
        idx = with(res_stan, model_type==mt[m] & sample_size==ss[s] & 
                     data_number==dn[d] & idx_mom)
        res_mom = res_stan[idx,]
        
        y_lim = range( with(res_mom, c(X5.5., X94.5., true)), na.rm=T )
        
        
        # plot 2: texts difficulties
        plot(1:nrow(res_mom), res_mom$mean, ylim=y_lim, xaxt='n', yaxt='n',
             col=col.alpha('blue', 0.3), pch=19, 
             xlab='', ylab='estimates', main='Text difficulty and deviation')
        abline(v=5+0.5, col=col.alpha('black', 0.5), lty=2)
        abline(h=0, col=col.alpha('black', 0.3), lty=2)
        axis(side=1, at=1:nrow(res_mom), labels=res_mom$parameter, las=2)
        axis(side=2, at=round( seq(y_lim[1], y_lim[2], by=0.2), 2), las=1)
        for(i in 1:nrow(res_mom)){
          lines(x=rep(i,2), y=with( res_mom[i,], c(X5.5., X94.5.) ) ,
                col=col.alpha('blue', 0.3))
        }
        points(1:nrow(res_mom), res_mom$true, col=col.alpha('black', 0.5))
        
        par(mfrow=c(1,1))
        
        dev.off()
        
        
        
        
        
        ### b. regression and contrasts 
        
        figure_name = paste0(mt[m],'_J',ss[s],'_Ndata',dn[d],'_regression.png')
        png( file.path(figure_save, figure_name), 
             units='cm', width=30, height=12, res=200)
        
        par(mfrow=c(1,2), mar=c(8,4,4,2))
        
        # identify parameter
        idx_mom = detect_par(res_stan, c('a','b_G','b_A','b_E','b_X'))
        idx = with(res_stan, model_type==mt[m] & sample_size==ss[s] & 
                     data_number==dn[d] & idx_mom)
        res_mom = res_stan[idx,]
        
        # sorting
        idx_ord = c(1,9:10,2,3:4,6,12:13,15,18, # it has a specific order
                    11,5,7:8,14,16:17,19:21)
        res_mom = res_mom[idx_ord,]
        
        y_lim = range( with(res_mom, c(X5.5., X94.5., true)), na.rm=T )
        
        # plot 1: regression
        plot(1:11, res_mom$mean[1:11], ylim=y_lim, xaxt='n', yaxt='n',
             col=col.alpha('blue', 0.3), pch=19, 
             xlab='', ylab='estimates', main='Regression parameters')
        # abline(v=c(2,3,6)+0.5, col=col.alpha('black', 0.5), lty=2)
        abline(h=0, col=col.alpha('black', 0.3), lty=2)
        axis(side=1, at=1:11, labels=res_mom$parameter[1:11], las=2)
        axis(side=2, at=round( seq(y_lim[1], y_lim[2], by=1), 2), las=1)
        for(i in 1:11){
          lines(x=rep(i,2), y=with( res_mom[i,], c(X5.5., X94.5.) ) ,
                col=col.alpha('blue', 0.3))
        }
        points(1:11, res_mom$true[1:11], col=col.alpha('black', 0.5))
        legend('topleft', legend=c('estimate','true'), pch=c(19,1), cex=.8, bty='n',
               col=c(col.alpha('blue', 0.3), col.alpha('black', 0.5)) ) 
        
        
        # plot 2: contrast
        plot(12:nrow(res_mom), res_mom$mean[12:nrow(res_mom)], ylim=y_lim, 
             xaxt='n', yaxt='n', col=col.alpha('blue', 0.3), pch=19, 
             xlab='', ylab='estimates', main='Contrasts parameters')
        abline(v=c(12,15)+0.5, col=col.alpha('black', 0.5), lty=2)
        abline(h=0, col=col.alpha('black', 0.3), lty=2)
        axis(side=1, at=12:nrow(res_mom), labels=res_mom$parameter[12:nrow(res_mom)], las=2)
        axis(side=2, at=round( seq(y_lim[1], y_lim[2], by=1), 2), las=1)
        for(i in 12:nrow(res_mom)){
          lines(x=rep(i,2), y=with( res_mom[i,], c(X5.5., X94.5.) ) ,
                col=col.alpha('blue', 0.3))
        }
        points(12:nrow(res_mom), res_mom$true[12:nrow(res_mom)], col=col.alpha('black', 0.5))
        
        par(mfrow=c(1,1), mar=opar$mar)
        
        dev.off()
        
        
        
        ### c. correlations 
        figure_name = paste0(mt[m],'_J',ss[s],'_Ndata',dn[d],'_corr.png')
        png( file.path(figure_save, figure_name), 
             units='cm', width=15, height=12, res=200)
        
        par(mfrow=c(1,1), mar=c(9,4,4,2))
        
        # identify parameter
        idx_mom = res_stan$parameter %in% c('Rho_theta_sub[1,2]','Rho_theta_sub[1,3]','Rho_theta_sub[2,3]')
        idx = with(res_stan, model_type==mt[m] & sample_size==ss[s] & 
                     data_number==dn[d] & idx_mom)
        res_mom = res_stan[idx,]
        
        y_lim = range( with(res_mom, c(X5.5., X94.5., true)), na.rm=T )
        
        # plot 1: correlations
        plot(1:nrow(res_mom), res_mom$mean, xlim=c(0, nrow(res_mom)+1), ylim=y_lim, 
             xaxt='n', yaxt='n', col=col.alpha('blue', 0.3), pch=19, 
             xlab='', ylab='estimates', main='Correlation parameters')
        abline(h=0, col=col.alpha('black', 0.3), lty=2)
        axis(side=1, at=1:nrow(res_mom), labels=res_mom$parameter, las=2)
        axis(side=2, at=round( seq(y_lim[1], y_lim[2], by=0.2), 2), las=1)
        for(i in 1:nrow(res_mom)){
          lines(x=rep(i,2), y=with( res_mom[i,], c(X5.5., X94.5.) ) ,
                col=col.alpha('blue', 0.3))
        }
        points(1:nrow(res_mom), res_mom$true, col=col.alpha('black', 0.5))
        legend('topleft', legend=c('estimate','true'), pch=c(19,1), cex=.8, bty='n',
               col=c(col.alpha('blue', 0.3), col.alpha('black', 0.5)) ) 
      
        par(mfrow=c(1,1), mar=opar$mar)
        
        dev.off()
        
        
        if( str_detect(mt[m], '^SOLV') ){
          
          figure_name = paste0(mt[m],'_J',ss[s],'_Ndata',dn[d],'_loads.png')
          png( file.path(figure_save, figure_name), 
               units='cm', width=15, height=12, res=200)
          
          par(mfrow=c(1,1), mar=c(9,4,4,2))
          
          # identify parameter
          idx_mom = detect_par(res_stan, 'loads')
          idx = with(res_stan, model_type==mt[m] & sample_size==ss[s] & 
                       data_number==dn[d] & idx_mom)
          res_mom = res_stan[idx,]
          
          y_lim = range( with(res_mom, c(X5.5., X94.5., true)), na.rm=T )
          
          # plot 1: correlations
          plot(1:nrow(res_mom), res_mom$mean, xlim=c(0, nrow(res_mom)+1), ylim=c(0, y_lim[2]), 
               xaxt='n', yaxt='n', col=col.alpha('blue', 0.3), pch=19, 
               xlab='', ylab='estimates', main='Loadings')
          abline(h=0, col=col.alpha('black', 0.3), lty=2)
          axis(side=1, at=1:nrow(res_mom), labels=res_mom$parameter, las=2)
          axis(side=2, at=round( seq(0, y_lim[2], by=0.2), 2), las=1)
          for(i in 1:nrow(res_mom)){
            lines(x=rep(i,2), y=with( res_mom[i,], c(X5.5., X94.5.) ) ,
                  col=col.alpha('blue', 0.3))
          }
          points(1:nrow(res_mom), res_mom$true, col=col.alpha('black', 0.5))
          legend('topleft', legend=c('estimate','true'), pch=c(19,1), cex=.8, bty='n',
                 col=c(col.alpha('blue', 0.3), col.alpha('black', 0.5)) ) 
          
          par(mfrow=c(1,1), mar=opar$mar)
          
          dev.off()
          
        }
        
        
        
        ### d. abilities parameters
        
        figure_name = paste0(mt[m],'_J',ss[s],'_Ndata',dn[d],'_abilities.png')
        png( file.path(figure_save, figure_name), 
             units='cm', width=30, height=20, res=200)
        
        par(mfrow=c(2,2), mar=c(9, 4, 4, 2))
    
        # identify parameter
        idx_mom = res_stan$parameter %in% paste0('theta_sub[',1:500,',1]')
        idx = with(res_stan, model_type==mt[m] & sample_size==ss[s] & 
                     data_number==dn[d] & idx_mom)
        res_mom = res_stan[idx,]
        
        # sample
        set.seed(5894) # to make it comparable across all conditions
        idx_s = sample( 1:nrow(res_mom), 25 )
        idx_s = idx_s[order(idx_s)]
        res_mom = res_mom[idx_s,]
        
        # sorting
        loc_name = str_locate(res_mom$parameter, '[:punct:][:digit:]{1,3}')
        num_name = as.integer(str_sub(res_mom$parameter, start=loc_name[,1]+1, end=loc_name[,2]))
        idx_ord = order(num_name)
        res_mom = res_mom[idx_ord,]
        
        y_lim = range( with(res_mom, c(X5.5.,X94.5., true)), na.rm=T )
        
        
        # plot 1: literal
        plot(1:25, res_mom$mean, ylim=y_lim, xaxt='n', yaxt='n',
             col=col.alpha('blue', 0.3), pch=19, 
             xlab='', ylab='estimates', main='Sub-dimension 1')
        abline(h=0, col=col.alpha('black', 0.3), lty=2)
        axis(side=1, at=1:25, labels=res_mom$parameter, las=2)
        axis(side=2, at=round( seq(y_lim[1], y_lim[2], by=1), 2), las=1)
        for(i in 1:25){
          lines(x=rep(i,2), y=with( res_mom[i,], c(X5.5.,X94.5.) ) ,
                col=col.alpha('blue', 0.3))
        }
        points(1:25, res_mom$true, col=col.alpha('black', 0.5))
        
        
        
        
        # identify parameter
        idx_mom = res_stan$parameter %in% paste0('theta_sub[',1:500,',2]')
        idx = with(res_stan, model_type==mt[m] & sample_size==ss[s] & 
                     data_number==dn[d] & idx_mom)
        res_mom = res_stan[idx,]
        res_mom = res_mom[idx_s,] # sample
        
        # sorting
        loc_name = str_locate(res_mom$parameter, '[:punct:][:digit:]{1,3}')
        num_name = as.integer(str_sub(res_mom$parameter, start=loc_name[,1]+1, end=loc_name[,2]))
        idx_ord = order(num_name)
        res_mom = res_mom[idx_ord,]
        
        y_lim = range( with(res_mom, c(X5.5.,X94.5., true)), na.rm=T )
        
        
        # plot 1: inferential
        plot(1:25, res_mom$mean, ylim=y_lim, xaxt='n', yaxt='n',
             col=col.alpha('blue', 0.3), pch=19, 
             xlab='', ylab='estimates', main='Sub-dimension 2')
        abline(h=0, col=col.alpha('black', 0.3), lty=2)
        axis(side=1, at=1:25, labels=res_mom$parameter, las=2)
        axis(side=2, at=round( seq(y_lim[1], y_lim[2], by=1), 2), las=1)
        for(i in 1:25){
          lines(x=rep(i,2), y=with( res_mom[i,], c(X5.5.,X94.5.) ) ,
                col=col.alpha('blue', 0.3))
        }
        points(1:25, res_mom$true, col=col.alpha('black', 0.5))
        
        
        
        
        # identify parameter
        idx_mom = res_stan$parameter %in% paste0('theta_sub[',1:500,',3]')
        idx = with(res_stan, model_type==mt[m] & sample_size==ss[s] & 
                     data_number==dn[d] & idx_mom)
        res_mom = res_stan[idx,]
        res_mom = res_mom[idx_s,] # sample
        
        # sorting
        loc_name = str_locate(res_mom$parameter, '[:punct:][:digit:]{1,3}')
        num_name = as.integer(str_sub(res_mom$parameter, start=loc_name[,1]+1, end=loc_name[,2]))
        idx_ord = order(num_name)
        res_mom = res_mom[idx_ord,]
        
        y_lim = range( with(res_mom, c(X5.5.,X94.5., true)), na.rm=T )
        
        
        # plot 1: reflective
        plot(1:25, res_mom$mean, ylim=y_lim, xaxt='n', yaxt='n',
             col=col.alpha('blue', 0.3), pch=19, 
             xlab='', ylab='estimates', main='Sub-dimension 3')
        abline(h=0, col=col.alpha('black', 0.3), lty=2)
        axis(side=1, at=1:25, labels=res_mom$parameter, las=2)
        axis(side=2, at=round( seq(y_lim[1], y_lim[2], by=1), 2), las=1)
        for(i in 1:25){
          lines(x=rep(i,2), y=with( res_mom[i,], c(X5.5.,X94.5.) ) ,
                col=col.alpha('blue', 0.3))
        }
        points(1:25, res_mom$true, col=col.alpha('black', 0.5))
        
        
        if( str_detect(mt[m], '^SOLV') ){
          
          # identify parameter
          idx_mom = res_stan$parameter %in% paste0('theta[',1:500,']')
          idx = with(res_stan, model_type==mt[m] & sample_size==ss[s] & 
                       data_number==dn[d] & idx_mom)
          res_mom = res_stan[idx,]
          res_mom = res_mom[idx_s,] # sample
          
          # sorting
          loc_name = str_locate(res_mom$parameter, '[:punct:][:digit:]{1,3}')
          num_name = as.integer(str_sub(res_mom$parameter, start=loc_name[,1]+1, end=loc_name[,2]))
          idx_ord = order(num_name)
          res_mom = res_mom[idx_ord,]
          
          y_lim = range( with(res_mom, c(X5.5.,X94.5., true)), na.rm=T )
          
          
          # plot 1: reflective
          plot(1:25, res_mom$mean, ylim=y_lim, xaxt='n', yaxt='n',
               col=col.alpha('blue', 0.3), pch=19, 
               xlab='', ylab='estimates', main='Higher-order dimension')
          abline(h=0, col=col.alpha('black', 0.3), lty=2)
          axis(side=1, at=1:25, labels=res_mom$parameter, las=2)
          axis(side=2, at=round( seq(y_lim[1], y_lim[2], by=1), 2), las=1)
          for(i in 1:25){
            lines(x=rep(i,2), y=with( res_mom[i,], c(X5.5.,X94.5.) ) ,
                  col=col.alpha('blue', 0.3))
          }
          points(1:25, res_mom$true, col=col.alpha('black', 0.5))
          
        }
        
        par(mfrow=c(1,1), mar=opar$mar)
        
        dev.off()
        
        # advance
        print( paste0(mt[m],'_J', ss[s], '_Ndata', dn[d]) )
        
      }
    }
  }
  
}




# function:
#     item_plots
# description:  
#     It generates the ICC and IIF from each condition
# arguments:
#     c_list = data.frame generated with file_id() function
#     file_save = location to save file
#     true_load = location of true data
#
item_plots = function(c_list, file_save, true_load){
  
  # # test
  # c_list = chains_list
  # file_save = file.path(getwd(), 'figures6')
  # true_load = file.path(getwd(), 'data')
  
  # figure parameters
  opar = par()
  
  # parameters
  mt = unique(c_list$model)
  ss = unique(c_list$sample)
  dn = unique(c_list$data)
  
  # plots
  # m=1
  # s=1
  # d=1
  for(m in 1:length(mt)){
    for(s in 1:length(ss)){
      for(d in 1:length(dn)){
      
        ####
        # 1. loading required parameters
        ####
        
        # loading true parameters
        file_name = paste0('Parameters_J',ss[s],'_l0.95_Ndata',dn[d],'.RData')
        load( file.path(true_load, file_name) )
        # data_true
        
        # identify files of stan
        idx_files = with(chains_list, model==mt[m] & sample==ss[s] & data==dn[d])
        fit_files = chains_list[idx_files, 1]
        stan_model = rstan::read_stan_csv( file.path(chains_path, fit_files) )
        
        # extracting samples
        set.seed(45589+m+s+d)
        post_sim = extract.samples( stan_model )
        # str(post_sim)
        
        
        
        ####
        # 2. plots
        ####
        
        ## 6.1 ICC
        figure_name = paste0(mt[m],'_J',ss[s],'_Ndata',dn[d],'_ICC.png')
        png( file.path(file_save, figure_name), 
             units='cm', width=55, height=80, res=100)
        
        par(mfrow=c(7,4))
        for(k in 1:25){
          
          # simulated, marginal, average, and true ICC
          ds_ICC = sim_ICC(sim_b=post_sim$b_k, item_id=k)
          ds_ICC_mar = mar_ICC(sim_ICC_object=ds_ICC)
          ds_ICC_ave = ave_ICC(sim_b=post_sim$b_k, item_id=k)
          ds_ICC_true = ICC(b=data_true$items$b[k])
          
          # plot
          plot(NULL, xlim=c(-3,3), ylim=c(0,1), main=paste0('item ', k),
               xlab='individual ability', ylab='probability of correct item')
          abline(h=0.5, v=0, lty=2, col=col.alpha('black', 0.5))
          for(i in 1:500){ 
            lines(ds_ICC$theta, ds_ICC[,i+1], col=col.alpha('black',0.03)) 
          }
          lines( ds_ICC_mar, col='blue', lwd=3, lty=2 )
          lines( ds_ICC_ave, col='red', lwd=2, lty=3 ) 
          lines( ds_ICC_true, col='black', lwd=3 ) 
          
        }
        par(mfrow=c(1,1))
        
        dev.off()
        
        
        
        ## 6.2 IIF
        figure_name = paste0(mt[m],'_J',ss[s],'_Ndata',dn[d],'_IIF.png')
        png( file.path(file_save, figure_name), 
             units='cm', width=55, height=80, res=100)
        
        par(mfrow=c(7,4))
        for(k in 1:25){
          
          # simulated, marginal, average, and true IIF
          ds_IIF = sim_IIF(sim_b=post_sim$b_k, item_id=k)
          ds_IIF_mar = mar_IIF(sim_IIF_object=ds_IIF)
          ds_IIF_ave = ave_IIF(sim_b=post_sim$b_k, item_id=k)
          ds_IIF_true = IIF(b=data_true$items$b[k])
          
          # plot
          plot(NULL, xlim=c(-3,3), ylim=c(0, 0.3), main=paste0('item ', k),
               xlab='individual ability', ylab='probability of correct item')
          abline(h=0.5, v=0, lty=2, col=col.alpha('black', 0.5))
          for(i in 1:500){ 
            lines(ds_IIF$theta, ds_IIF[,i+1], col=col.alpha('black',0.03)) 
          }
          lines( ds_IIF_mar, col='blue', lwd=2, lty=2 )
          lines( ds_IIF_ave, col='red', lwd=2, lty=2 ) 
          lines( ds_IIF_true, col='black', lwd=3 )
          
        }
        par(mfrow=c(1,1))
        
        dev.off()
        
        # advance
        print( paste0(mt[m],'_J', ss[s], '_Ndata', dn[d]) )
        
      }
    }
  }
}       




# function:
#     retrodictive_agg
# description:  
#     It generates aggregation of data 
# arguments:
#     c_list = data.frame generated with file_id() function
#     file_save = location to save file
#     true_load = location of true data
#
retrodictive_agg = function(c_list, file_save, true_load){
  
  # # test
  # c_list = chains_list
  # file_save = file.path(getwd(), 'figures6')
  # true_load = file.path(getwd(), 'data')
  
  # figure parameters
  opar = par()
  
  # parameters
  mt = unique(c_list$model)
  ss = unique(c_list$sample)
  dn = unique(c_list$data)
  
  # plots
  # m=1
  # s=1
  # d=1
  for(m in 1:length(mt)){
    for(s in 1:length(ss)){
      for(d in 1:length(dn)){
        
        ####
        # 1. loading required parameters
        ####
        
        # loading true parameters
        file_name = paste0('ListFormat_J',ss[s],'_l0.95_Ndata',dn[d],'.RData')
        load( file.path(true_load, file_name) )
        # data_post
        
        # identify files of stan
        idx_files = with(c_list, model==mt[m] & sample==ss[s] & data==dn[d])
        fit_files = c_list[idx_files, 1]
        stan_model = rstan::read_stan_csv( file.path(chains_path, fit_files) )
        
        # extracting samples
        set.seed(7856+m+s+d)
        post_sim = extract.samples( stan_model )
        # str(post_sim)
        
        
        ### 6.3 individuals hit rate
        ds_prob_post = sim_hit_rate(d_true=data_post, d_sim=post_sim, prob=T)
        file_name = paste0(mt[m],'_J',ss[s],'_Ndata',dn[d],'_HitRate_prob.RData')
        save(ds_prob_post, file=file.path(file_save, file_name))
        
        ds_out_post = sim_hit_rate(d_true=data_post, d_sim=post_sim, prob=F)
        file_name = paste0(mt[m],'_J',ss[s],'_Ndata',dn[d],'_HitRate_out.RData')
        save(ds_out_post, file=file.path(file_save, file_name))
        # str(post_sim)
        
        # advance
        print( paste0(mt[m],'_J', ss[s], '_Ndata', dn[d]) )
        
      }
    }
  }
}

 

# function:
#     retrodictive_plots
# description:  
#     It generates aggregation of data 
# arguments:
#     c_list = data.frame generated with file_id() function
#     file_save = location to save file
#     true_load = location of true data
#
retrodictive_plots = function(c_list, file_save, true_load){
  
  # # test
  # c_list = chains_list
  # file_save = file.path(getwd(), 'figures6')
  # true_load = file.path(getwd(), 'figures6')
  
  # figure parameters
  opar = par()
  
  # parameters
  mt = unique(c_list$model)
  ss = unique(c_list$sample)
  dn = unique(c_list$data)
  
  # plots
  # m=1
  # s=1
  # d=1
  for(m in 1:length(mt)){
    for(s in 1:length(ss)){
      for(d in 1:length(dn)){
        
        ####
        # 1. loading required parameters
        ####
        
        # loading stats
        file_name = paste0(mt[m],'_J',ss[s],'_Ndata',dn[d],'_HitRate_prob.RData')
        load( file.path(true_load, file_name) )
        
        file_name = paste0(mt[m],'_J',ss[s],'_Ndata',dn[d],'_HitRate_out.RData')
        load( file.path(true_load, file_name) )
        
        
        ####
        # 2. plots
        ####
        
        # 2.1 individuals
        
        # average and marginal hit rate
        ds_ave = ave_hit_rate(sim_hit_object=ds_prob_post$IDind, prob=0.95)
        ds_mar = ave_hit_rate(sim_hit_object=ds_out_post$IDind, prob=0.95)
        
        # random sample
        n = 50
        set.seed(2546) # to make all comparable
        idx = sample(ds_ave$IDj, size=n)
        idx = idx[order(idx)]
        
        ds_s = ds_prob_post$IDind[ ds_prob_post$IDind$IDj %in% idx,]
        ds_ave_s = ds_ave[ds_ave$IDj %in% idx,]
        ds_mar_s = ds_mar[ds_mar$IDj %in% idx,]
        
        
        # plot
        figure_name = paste0(mt[m],'_J',ss[s],'_Ndata',dn[d],'_HitRate_ind.png')
        png( file.path(file_save, figure_name), 
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
        
        
        
        # observed hit rates per variable
        figure_name = paste0(mt[m],'_J',ss[s],'_Ndata',dn[d],'_HitRate_var.png')
        png( file.path(file_save, figure_name), 
             units='cm', width=30, height=25, res=100)

        par(mfrow=c(2,2))
        
        # gender
        plot(ds_s$G-0.1, ds_s$true, xlim=with(ds_s, c(min(G)-1, max(G)+1)), ylim=c(0,1), 
             xaxt='n', yaxt='n', col=col.alpha('black', 0.5), pch=1,
             xlab='', ylab='% correct items', main='gender')
        axis(side=1, at=with(ds_s, min(G):max(G)), labels=c('male','female') ) 
        axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
        for(j in (1:200)+1 ){
          points(ds_s[,c('G', paste0('X',j) )], col=col.alpha('blue', 0.002), pch=19 )
        }
        legend('topleft', legend=c('observed','predicted'), pch=c(1,19,19), 
               cex=.8, bty='n', col=c('black', col.alpha('blue',0.3), col.alpha('red',0.3) ) )
        
        
        # age
        plot(ds_s$A-0.2, ds_s$true, xlim=with(ds_s, c(min(A)-1, max(A)+1)), ylim=c(0,1), 
             xaxt='n', yaxt='n', col=col.alpha('black', 0.5), pch=1,
             xlab='', ylab='% correct items', main='age')
        axis(side=1, at=with(ds_s, min(A):max(A)) ) 
        axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
        for(j in (1:500)+1 ){
          points(ds_s[,c('A', paste0('X',j) )], col=col.alpha('blue', 0.002), pch=19 )
        }
        
        # edu
        plot(ds_s$E-0.1, ds_s$true, xlim=with(ds_s, c(min(E)-1, max(E)+1)), ylim=c(0,1), 
             xaxt='n', yaxt='n', col=col.alpha('black', 0.5), pch=1,
             xlab='', ylab='% correct items', main='education')
        axis(side=1, at=with(ds_s, min(E):max(E)), labels=c('inst','uni','both') ) 
        axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
        for(j in (1:200)+1 ){
          points(ds_s[,c('E', paste0('X',j) )], col=col.alpha('blue', 0.002), pch=19 )
        }
        
        
        # exp
        plot(ds_s$X-0.1, ds_s$true, xlim=with(ds_s, c(min(X)-1, max(X)+1)), ylim=c(0,1), 
             xaxt='n', yaxt='n', col=col.alpha('black', 0.5), pch=1,
             xlab='', ylab='% correct items', main='experience')
        axis(side=1, at=with(ds_s, min(X):max(X)), 
             labels=c('0y','1-5y','6-10y','+10y') ) 
        axis(side=2, at=seq(0, 1, by=0.1), las=1 ) 
        for(j in (1:200)+1 ){
          points(ds_s[,c('X', paste0('X',j) )], col=col.alpha('blue', 0.002), pch=19 )
        }
        
        par(mfrow=c(1,1))
        
        dev.off()
        
        
        
        
        # 2.2 items
        
        # average and marginal hit rate
        ds_s = ds_prob_post$IDitems
        ds_ave_s = ave_hit_rate(sim_hit_object=ds_prob_post$IDitems, prob=0.95)
        ds_mar_s = ave_hit_rate(sim_hit_object=ds_out_post$IDitems, prob=0.95)
        
        
        # plot
        figure_name = paste0(mt[m],'_J',ss[s],'_Ndata',dn[d],'_HitRate_items.png')
        png( file.path(file_save, figure_name), 
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
        
        
        
        
        # 2.3 texts
        
        # average and marginal hit rate
        ds_s = ds_prob_post$IDtext
        ds_ave_s = ave_hit_rate(sim_hit_object=ds_prob_post$IDtext, prob=0.95)
        ds_mar_s = ave_hit_rate(sim_hit_object=ds_out_post$IDtext, prob=0.95)
        
        
        # plot
        figure_name = paste0(mt[m],'_J',ss[s],'_Ndata',dn[d],'_HitRate_texts.png')
        png( file.path(file_save, figure_name), 
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
        
        
        
        
        # 2.4 dimensions
        
        # average and marginal hit rate
        ds_s = ds_prob_post$IDdim
        ds_ave_s = ave_hit_rate(sim_hit_object=ds_prob_post$IDdim, prob=0.95)
        ds_mar_s = ave_hit_rate(sim_hit_object=ds_out_post$IDdim, prob=0.95)
        
        
        # plot
        figure_name = paste0(mt[m],'_J',ss[s],'_Ndata',dn[d],'_HitRate_dims.png')
        png( file.path(file_save, figure_name), 
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
        
        # advance
        print( paste0(mt[m],'_J', ss[s], '_Ndata', dn[d]) )
        
      }
    }
  }
}



# function:
#     retrodictive_rmse
# description:  
#     It generates rmse aggregated by diverse dimensions 
# arguments:
#     c_list = data.frame generated with file_id() function
#     file_save = location to save file
#     true_load = location of true data
#
retrodictive_rmse = function(c_list, file_save, true_load){
  
  # # test
  # c_list = chains_list
  # file_save = file.path(getwd(), 'figures6')
  # true_load = file.path(getwd(), 'figures6')
  
  # parameters
  mt = unique(c_list$model)
  ss = unique(c_list$sample)
  dn = unique(c_list$data)
  
  # plots
  # m=1
  # s=1
  # d=1
  for(m in 1:length(mt)){
    for(s in 1:length(ss)){
      for(d in 1:length(dn)){
        
        ####
        # 1. loading required parameters
        ####
        
        # loading stats
        file_name = paste0(mt[m],'_J',ss[s],'_Ndata',dn[d],'_HitRate_prob.RData')
        load( file.path(true_load, file_name) )
        
        file_name = paste0(mt[m],'_J',ss[s],'_Ndata',dn[d],'_HitRate_out.RData')
        load( file.path(true_load, file_name) )
        
        
        ####
        # 2. storage
        ####
        
        # 2.1 individuals
        idx_var = c('IDj','G','A','E','X','RMSE','Sdiff')
        rmse_mom = data.frame(model_type = mt[m],
                              sample_size = ss[s],
                              data_number = dn[d],
                              type = 'prob',
                              ds_prob_post$IDind[,idx_var])
        rmse_mom = rbind(rmse_mom,
                         data.frame(model_type = mt[m],
                                    sample_size = ss[s],
                                    data_number = dn[d],
                                    type = 'out',
                                    ds_out_post$IDind[,idx_var]) )
        names(rmse_mom)[10] = 'RMSE_within'
        
        if(all(m==1, s==1, d==1)){
          rmse_IDind = rmse_mom
        } else{
          rmse_IDind = rbind(rmse_IDind, rmse_mom)
        }
        
        
        # 2.2 items
        idx_var = c('IDk','RMSE','Sdiff')
        rmse_mom = data.frame(model_type = mt[m],
                              sample_size = ss[s],
                              data_number = dn[d],
                              type = 'prob',
                              ds_prob_post$IDitems[,idx_var])
        rmse_mom = rbind(rmse_mom,
                         data.frame(model_type = mt[m],
                                    sample_size = ss[s],
                                    data_number = dn[d],
                                    type = 'out',
                                    ds_out_post$IDitems[,idx_var]) )
        names(rmse_mom)[6] = 'RMSE_within'
        
        if(all(m==1, s==1, d==1)){
          rmse_IDitems = rmse_mom
        } else{
          rmse_IDitems = rbind(rmse_IDitems, rmse_mom)
        }
        
        
        # 2.3 texts
        idx_var = c('IDl','RMSE','Sdiff')
        rmse_mom = data.frame(model_type = mt[m],
                              sample_size = ss[s],
                              data_number = dn[d],
                              type = 'prob',
                              ds_prob_post$IDtext[,idx_var])
        rmse_mom = rbind(rmse_mom,
                         data.frame(model_type = mt[m],
                                    sample_size = ss[s],
                                    data_number = dn[d],
                                    type = 'out',
                                    ds_out_post$IDtext[,idx_var]) )
        names(rmse_mom)[6] = 'RMSE_within'
        
        if(all(m==1, s==1, d==1)){
          rmse_IDtext = rmse_mom
        } else{
          rmse_IDtext = rbind(rmse_IDtext, rmse_mom)
        }
        
        
        # 2.4 dimensions
        idx_var = c('IDd','RMSE','Sdiff')
        rmse_mom = data.frame(model_type = mt[m],
                              sample_size = ss[s],
                              data_number = dn[d],
                              type = 'prob',
                              ds_prob_post$IDdim[,idx_var])
        rmse_mom = rbind(rmse_mom,
                         data.frame(model_type = mt[m],
                                    sample_size = ss[s],
                                    data_number = dn[d],
                                    type = 'out',
                                    ds_out_post$IDdim[,idx_var]) )
        names(rmse_mom)[6] = 'RMSE_within'
        
        if(all(m==1, s==1, d==1)){
          rmse_IDdim = rmse_mom
        } else{
          rmse_IDdim = rbind(rmse_IDdim, rmse_mom)
        }
        
      }
    }
  }
  
  
  # saving previous storage
  save(rmse_IDind, file=file.path(file_save, 'rmse_IDind.RData') )
  save(rmse_IDitems, file=file.path(file_save, 'rmse_IDitems.RData') )
  save(rmse_IDtext, file=file.path(file_save, 'rmse_IDtext.RData') )
  save(rmse_IDdim, file=file.path(file_save, 'rmse_IDdim.RData') )
  
  
  # calculate mean within and between 
  rmse_final = rmse_IDind %>%
    group_by(model_type, sample_size, type, IDj) %>%
    summarize(n=n(), 
              mean_RMSE_within = mean(RMSE_within),
              RMSE_between = sqrt(mean(Sdiff)) )
  save(rmse_final, file=file.path(file_save, 'rmse_WB_IDind.RData') )
  
  rmse_final = rmse_IDitems %>%
    group_by(model_type, sample_size, type, IDk) %>%
    summarize(n=n(), 
              mean_RMSE_within = mean(RMSE_within),
              RMSE_between = sqrt(mean(Sdiff)) )
  save(rmse_final, file=file.path(file_save, 'rmse_WB_IDitems.RData') )
  
  rmse_final = rmse_IDtext %>%
    group_by(model_type, sample_size, type, IDl) %>%
    summarize(n=n(), 
              mean_RMSE_within = mean(RMSE_within),
              RMSE_between = sqrt(mean(Sdiff)) )
  save(rmse_final, file=file.path(file_save, 'rmse_WB_IDtext.RData') )
  
  rmse_final = rmse_IDdim %>%
    group_by(model_type, sample_size, type, IDd) %>%
    summarize(n=n(), 
              mean_RMSE_within = mean(RMSE_within),
              RMSE_between = sqrt(mean(Sdiff)) )
  save(rmse_final, file=file.path(file_save, 'rmse_WB_IDdim.RData') )
  
}