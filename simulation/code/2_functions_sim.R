
# function:
#     data_generation
# description:  
#     To generate data based on different parameter settings.
#     Only two parameters are effectively controlled in the experimentation: 
#     sample size (J), and the loading from the SOLV to the FOLV (loads).
# characteristics of the sample design:
#   - one instrument
#   - hierarchical measurement scales with FOLV and SOLV
#   - one evaluation time
#   - one sample of individuals
#   - testlet items, multiple items come from one text
#   - with covariates
#   - NO missingness
# arguments:
#     J = individual sample sizes
#     loads = loadings from SOLV to FOLV (it control correlation)
#     Ndata = defines the number of data generated
#     file_dir = path to save generated data
#     s_theta = sd to simulate SOLV and FOLV
#     s_text = sd in generating items from a specific text
#     D = number of dimensions (default 3)
#     K = number of items (default 25)
#     L = number of texts (default 5), it has to be a multiple of K
#     seed = seed used to generate the simulation (default 1)
#     prec = rounding in abilities and item generation parameters (default 3)

data_generation =function( J=100, loads=rep(0.1, 3), Ndata=1, file_dir,
                           s_theta=0.5, s_text=0.5,
                           D=3, K=25, L=5, seed=1, prec=3){

  # # test
  # J=100
  # D=3
  # K=25
  # L=5
  # loads=c(0.1, 0.2, 0.3)
  # Ndata=1
  # s_theta=0.5
  # s_text=0.2
  # seed=1
  # prec=3

  
  #_______________
  # 1. generation 
  #_______________
  set.seed(seed)
  
  
  ## 1.1. regression parameters
  mom = expand_grid(a=loads[-3], b=loads[-1])
  mom = mom[-3,]
  
  betas = list( gender = c(0, 0.5),
                age = -0.02,
                edu = c(-0.5, 0.5, 0),
                exp = c(-0.5, 0, 0.35, 0.5),
                loads = loads, # loadings
                exp_corr = with(mom, a*b) ) # expected correlation
  
  
  ## 1.2 covariates
  abilities = data.frame( IDind=1:J,
                          gender = sample(c(1,2), size=J, replace=T),
                          age = sample(30:65, size=J, replace=T),
                          edu = sample(c(1,2,3), size=J, replace=T),
                          exp = sample(c(1,2,3,4), size=J, replace=T),
                          theta=rep(NA,J), theta1=rep(NA,J), 
                          theta2=rep(NA,J), theta3=rep(NA,J) )
  
  # variable indices
  SOLV_loc = which(names(abilities)=='theta')
  FOLV_loc = which(names(abilities)=='theta1'):which(names(abilities)==paste0('theta', D))
  
  
  ## 1.3. abilities
  # SOLV
  m_theta = rep(NA, J)
  for(j in 1:J){
    ageC = with(abilities, age[j] - min(age) + 1 )
    m_theta[j] = with(abilities,
                      betas$gender[ gender[j] ] + 
                        betas$age * ageC +
                        betas$edu[ edu[j] ] + 
                        betas$exp[ exp[j] ] )  
  }
  abilities[, SOLV_loc] = round( rnorm(J, m_theta, s_theta), prec)
  
  # FOLV
  s_mult = diag( rep(s_theta, D) ) # independence after considering SOLV
  m_mult = data.frame( theta1=rep(NA, J), theta2=rep(NA, J), theta3=rep(NA, J) )
  for(j in 1:J){
    m_mult[j,] = betas$loads * abilities$theta[j]
    abilities[j, FOLV_loc] = round( mvrnorm(n=1, mu=unlist(m_mult[j,]), Sigma=s_mult), prec)
  }

  
  ## 1.4 texts
  texts = data.frame(IDtext=1:L, m_b=rep(NA, L), s_b=rep(NA, L))
  texts$m_b = seq(-1.5, 1.5, length.out=L)
  texts$s_b = rep(s_text, L)
  
  
  ## 1.5 items
  items = data.frame(IDitem=1:K, IDtext=rep(NA,K), IDdim=rep(NA,K), b=rep(NA, K))
  kl = K/L # items per text
  for(k in 1:nrow(texts)){
    items$b[(1:kl) + (k-1)*kl] = with(texts, round( rnorm(kl, m_b[k], s_b[k]), prec ) )
    items$IDtext[(1:kl) + (k-1)*kl] = k
  }
  
  
  ## 1.6 dimensions
  items$IDdim = sample(1:3, K, replace=T)
  
  
  
  #_______________
  # 2. storage
  #_______________
  
  ## 2.1 parameters
  data_true = list(
    
    seed=seed,
    
    # indices
    J = J,
    D = D, 
    K = K, 
    L = L,
    
    # abilities
    betas = betas,
    abilities = abilities,
    
    # texts and items
    texts = texts,
    items = items
    
  )
  file_name = paste0('Parameters_J',J,'_l',loads[1],'_Ndata',Ndata,'.RData')
  save(data_true, file=file.path(file_dir, file_name) )
  
  

  ## 2.2 long format
  data_eval = data.frame(
    
    # individual data
    IDind = rep(abilities$IDind, K),
    gender = rep(abilities$gender, K),
    age = rep(abilities$age, K),
    edu = rep(abilities$edu, K),
    exp = rep(abilities$exp, K),
    theta = rep(abilities$theta, K),
    theta1 = rep(abilities$theta1, K),
    theta2 = rep(abilities$theta2, K),
    theta3 = rep(abilities$theta3, K),
    
    # items
    IDitem = rep(items$IDitem, each=J),
    IDtext = rep(items$IDtext, each=J),
    IDdim = rep(items$IDdim, each=J),
    b = rep(items$b, each=J) 
    
  )
  
  # appropriate ability
  data_eval$theta_jkld = NA
  for(i in 1:nrow(data_eval) ){
    data_eval$theta_jkld[i] = data_eval[i, FOLV_loc[ data_eval$IDdim[i] ] ]
  }
  
  # linear predictor, probability, and outcome
  data_eval$v_jkld = with(data_eval, theta_jkld - b)
  data_eval$p_jkld = inv_logit( data_eval$v )
  data_eval$y_jkld = rbinom( nrow(data_eval), 1, prob=data_eval$p)

  file_name = paste0('LongFormat_J',J,'_l',loads[1],'_Ndata',Ndata,'.RData')
  save(data_eval, file=file.path(file_dir, file_name) )
  
  
  
  ## 2.3 estimation data
  data_post = list(
    
    # evaluation data
    N = nrow(data_eval),
    J = J,
    K = K,
    L = L,
    D = D,
    IDj = data_eval$IDind,
    IDk = data_eval$IDitem,
    IDl = data_eval$IDtext,
    IDd = data_eval$IDdim,
    GE = data_eval$gender,
    AG = data_eval$age,
    ED = data_eval$edu,
    XP = data_eval$exp,
    y = data_eval$y_jkld,
    
    # individual data
    IDind = abilities$IDind,
    G = abilities$gender,
    A = abilities$age,
    E = abilities$edu,
    X = abilities$exp,
    
    # items
    IDitem = items$IDitem,
    IDtext = items$IDtext,
    IDdim = items$IDdim
    
  )
  
  file_name = paste0('ListFormat_J',J,'_l',loads[1],'_Ndata',Ndata,'.RData')
  save(data_post, file=file.path(file_dir, file_name) )
  
}




# function:
#     run_prior
# description:  
#     To run each prior simulation for each model
# arguments:
#     model_path = path where all models are located
#     model_out = path where chains have to be saved
#     data_path = path where data per condition is located

run_prior = function(model_path, model_out, data_path){
  
  # # tests
  # model_path = file.path(getwd(), 'models_prior')
  # model_out = file.path(getwd(), 'priors')
  # data_path = file.path(getwd(), 'data')
  
  # model list
  model_list = list.files( model_path )
  model_list = model_list[ str_detect(model_list, '.stan') ]
  # length(model_list)
  
  # data list
  data_list = list.files( data_path )
  data_list = data_list[str_detect(data_list, 'ListFormat')]
  # length(data_list)
  
  # j=1
  for(j in 1:length(model_list) ){
    
    # compile model
    set_cmdstan_path('~/cmdstan')
    mod = cmdstan_model( file.path(model_path, model_list[j] ) )
    
    # generate base name
    base_nam = str_replace( data_list[1], 'ListFormat_', 
                            str_replace( model_list[j], '.stan', '_') )
    base_nam = str_replace( base_nam, '.RData', '' )
    
    # load data
    load( file.path(data_path, data_list[1] ) )
    
    # run model
    mod$sample(data = data_post, 
               output_dir = model_out, 
               output_basename = base_nam,
               chains=1, parallel_chains=1, init=0, adapt_delta=0.99)
  }
  
}




# function:
#     run_post
# description:  
#     To run each model for each data in all conditions
# arguments:
#     model_path = path where all models are located
#     model_out = path where chains have to be saved
#     data_path = path where data per condition is located

run_post = function(model_path, model_out, data_path){
  
  # # test
  # model_path = file.path(getwd(), 'models_post')
  # model_out = file.path(getwd(), 'chains_post')
  # data_path = file.path(getwd(), 'data')
  
  # model list
  model_list = list.files( model_path )
  model_list = model_list[ str_detect(model_list, '.stan') ]
  # length(model_list)
  
  # data list
  data_list = list.files( data_path )
  data_list = data_list[str_detect(data_list, 'ListFormat')]
  # length(data_list)
  
  # j=1
  # i=2
  for(j in 1:length(model_list) ){
    
    # compile model
    set_cmdstan_path('~/cmdstan')
    mod = cmdstan_model( file.path(model_path, model_list[j] ) )
    
    for(i in 1:length(data_list) ){
      
      # generate base name
      base_nam = str_replace( data_list[i], 'ListFormat_', 
                              str_replace(model_list[j], '.stan', '_') )
      base_nam = str_replace( base_nam, '.RData', '' )
      
      # load data
      load( file.path(data_path, data_list[i] ) )
      
      # run model
      t0 = proc.time()
      mod$sample(data = data_post, 
                 output_dir = model_out, 
                 output_basename = base_nam,
                 chains=3, parallel_chains=3, init=0, adapt_delta=0.99)
      t1 = proc.time()
      
      # saving time
      start = str_locate(base_nam, 'J')[2] 
      J = str_sub( base_nam, start=start+1, end=(start+3) )
      
      start = str_locate(base_nam, 'l')[2] 
      l = str_sub( base_nam, start=start+1, end=(start+3) )
      
      start = str_locate(base_nam, 'Ndata')[2]
      S = str_sub( base_nam, start=start+1, end=(start+2) )
      
      elapsed = (t1-t0)
      
      if(j==1 & i==1){
        time_elap = data.frame( J=J, load=l, data=S, time=unlist(elapsed['elapsed']) )
      } else{
        time_elap = rbind(time_elap, c(J, l, S, unlist(elapsed['elapsed']) ) )
      }
      
      # save time elapsed (saved at each iteration)
      write.csv( time_elap, row.names=F, 
                 file=file.path(model_out, 'time_elapsed.csv' ) )
      
    }
  }

}
