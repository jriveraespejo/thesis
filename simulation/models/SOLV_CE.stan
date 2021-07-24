
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
    int GE[N];
    real AG[N];
    int ED[N];
    int XP[N];
    int y[N];
    
    // individuals
    int IDind[J];
    int G[J];
    real A[J];
    int E[J];
    int X[J];
 
    // items
    int IDitem[K];
    int IDtext[K];
    int IDdim[K];
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
    real b_X[4];
    
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
    b_A ~ normal(0, 0.5);
    b_E ~ normal(0, 1);
    b_X ~ normal(0, 0.5);
    for(j in 1:J){
      m_theta[j] = a + b_G[ G[j] ] + 
                    b_A * ( A[j] - min(A) ) +
                    b_E[ E[j] ] + 
                    b_X[ X[j] ];
    }
    theta ~ normal(m_theta, 1); 
    // setting the scale at this level
    
    // FOLV
    loads ~ lognormal(0, 0.5);
    for(j in 1:J){
      m_mult[j,] = [ loads[1]*theta[j], 
                      loads[2]*theta[j], 
                      loads[3]*theta[j] ]';
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
//# generated quantities{
//#     vector[N] log_lik;
//#     real v;
//#     real p;
//#     
//#     // likelihood
//#     for( i in 1:N ) {
//#       v = theta_sub[ IDj[i], IDd[i] ] - b_k[ IDk[i] ];
//#       p = inv_logit(v);
//#       log_lik[i] = bernoulli_lpmf( y[i] | p);
//#     }
//# }

