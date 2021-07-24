
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
    real zb_k[K];

    // betas
    real a;
    real b_G[2];
    real b_A;
    real b_E[3];
    real b_X[4];
    
    // abilities
    vector[J] ztheta;             // reading comprehension
    real<lower=0> loads[D];       // loadings
    cholesky_factor_corr[D] L_Rho_theta_sub; // sub-dimensions: lit, inf, ref
    matrix[D, J] ztheta_sub; 
}
transformed parameters{
    real b_k[K];
    vector[J] m_theta;
    vector[J] theta;             // reading comprehension
    matrix[D, D] Rho_theta_sub;
    matrix[J, D] m_mult;
    matrix[J, D] theta_sub;
    
    // items
    for(k in 1:K){
      b_k[k] = m_b[ IDtext[k] ] + s_b[ IDtext[k] ] * zb_k[k]; 
    }
    
    // SOLV
    for(j in 1:J){
      m_theta[j] = a + b_G[ G[j] ] + 
                    b_A * ( A[j] - min(A) ) +
                    b_E[ E[j] ] + 
                    b_X[ X[j] ];
    }
    theta = m_theta + 1 * ztheta;
    // setting the scale at this level
    
    // FOLV
    Rho_theta_sub = multiply_lower_tri_self_transpose(L_Rho_theta_sub);
    for(j in 1:J){
      m_mult[j,] = to_row_vector( [ loads[1]*theta[j], 
                                    loads[2]*theta[j], 
                                    loads[3]*theta[j] ] );
    }
    theta_sub = m_mult + (L_Rho_theta_sub * ztheta_sub)';
    // also setting scale here
}
model{
    // declare
    real v;
    real p;
              
    // items
    m_b ~ normal(0, 1);
    s_b ~ exponential(2);
    zb_k ~ normal(0, 1);
    
    // abilities
    a ~ normal(0, 0.5);
    b_G ~ normal(0, 0.5);
    b_A ~ normal(0, 0.5);
    b_E ~ normal(0, 1);
    b_X ~ normal(0, 0.5);
    ztheta ~ normal(0,1);
    L_Rho_theta_sub ~ lkj_corr_cholesky(2);
    loads ~ lognormal(0, 0.5);
    to_vector(ztheta_sub) ~ normal(0,1);
    
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

