
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
#     trank_plot
# description:  
#     It creates a trank plot for a jags object 
# arguments:
#     jags_object = object produced after using run.jags
#     pars = character vector with the names of parameters to choose
#     wide = controls the number of iterations (in the chain) considered.
#           (default 50)
#
trank_plot = function(jags_object, pars, wide=50){
  
  # # test
  # jags_object = jags_model
  # pars = idx
  # wide = 50
  
  # for colors
  require(RColorBrewer)
  require(rethinking)
  
  # previous details
  jags_mcmc = jags_object$mcmc[,pars]
  n_chains = length(jags_mcmc)
  n_pars = length(pars)
  
  jags_summ = summary(jags_object)
  jags_SSeff = jags_summ[ rownames(jags_summ) %in% pars, 'SSeff'] 
  jags_SSeff = round( jags_SSeff / n_chains, 0 )
  
  ranks = list()
  xrange = rep(1:wide, each=2)
  yrange = vector('list', 4)
  
  # plot
  for(c in 1:n_chains){
    ranks[[c]] = apply(jags_mcmc[[c]][1:(wide+1),], 2, rank)
    
    for(p in 1:n_pars){
      yran = c()
      for(i in 2:(wide+1)){
        yran = c(yran, c(ranks[[c]][i-1, p], ranks[[c]][i, p] ) )
      }
      
      if(p==1){
        y_ran = yran 
      } else{
        y_ran = cbind(y_ran, yran)
      }
      
    }
    yrange[[c]] = y_ran
  }
  
  # plot
  trank_col = rethink_palette[1:n_chains]
  
  for(p in 1:n_pars){
    plot(NULL, type='l', xlim=c(0, wide+1), ylim=c(0, wide+1),
         xlab="", ylab="", xaxt ="n", yaxt ="n", axes=F)
    box(bty="l")
    for(c in 1:n_chains){
      lines(xrange, yrange[[c]][,p], col=trank_col[c], lwd=1.5)  
    }
    mtext( names(jags_SSeff)[p], side=3, adj=0, cex=1)
    mtext( paste("n_eff =", jags_SSeff[p]), side=3, adj=1, cex=1)
  }
}



# function:
#     acf_plot
# description:  
#     It creates a acf plot for a mcmc object 
# arguments:
#     mcmc_object = object containing either stanfit or mcmc (jags) object   
#     pars = character vector with the names of parameters to choose
#
acf_plot = function(mcmc_object, pars){
  
  # mcmc_object = stan_model
  n_pars = length(pars)
  
  if(class(mcmc_object) == 'stanfit'){
    sim_object = mcmc_object@sim[[1]][[1]]
    
    for(p in 1:n_pars){
      acf(sim_object[pars[p]][[1]], main='', 
          xlab='', ylab='', mar = c(0, 0, 0, 0) )
      mtext( pars[p], side=3, adj=0, cex=1)
      # mtext( paste("n_eff =", sim_object[p]), side=3, adj=1, cex=1)
    }
    
  } else{
    sim_object = mcmc_object$mcmc[[1]]
    
    for(p in 1:n_pars){
      acf(sim_object[colnames(sim_object) == pars[p]], 
          main='', xlab='', ylab='', mar = c(0, 0, 0, 0) )
      mtext( pars[p], side=3, adj=0, cex=1)
      # mtext( paste("n_eff =", sim_object[p]), side=3, adj=1, cex=1)
    }
    
  }
  
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