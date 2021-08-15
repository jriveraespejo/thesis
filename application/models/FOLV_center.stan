
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
    int G[N];
    real A[N];
    int E[N];
    int S[N];
    int Xpu[N];
    int Xpr[N];
    int Di[N];
    int y[N];
    
    // individuals
    int IDind[J];
    int GE[J];
    real AG[J];
    int ED[J];
    int SP[J];
    int XPpu[J];
    int XPpr[J];
    int DI[J];
 
    // items
    int IDitem[K];
    int IDtext[K];
    int IDdim[K];
} 
transformed data {
  real min_AG;
  min_AG = min(AG);
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
    real b_S[3];
    real b_Xpu[4];
    real b_Xpr[4];
    real b_Di[4];
    
    // abilities
    corr_matrix[D] Rho_theta_sub; // sub-dimensions: lit, inf, ref
    vector[D] theta_sub[J]; 
}
model{
    // declare
    vector[D] m_mult[J];
    real v;
    real p;
              
    // items
    m_b ~ normal(0, 1);
    s_b ~ exponential(2);
    for(k in 1:K){            // priors
      b_k[k] ~ normal( m_b[ IDtext[k] ], s_b[ IDtext[k] ]); 
    }
    
    // abilities
    a ~ normal(0, 0.5);
    b_G ~ normal(0, 0.5);
    b_A~ normal(0, 0.5);
    b_E ~ normal(0, 1);
    b_S ~ normal(0, 0.5);
    b_Xpu ~ normal(0, 0.5);
    b_Xpr ~ normal(0, 0.5);
    b_Di ~ normal(0, 0.5);
    Rho_theta_sub ~ lkj_corr(2);
    for(j in 1:J){
      m_mult[j,] = rep_vector( a + b_G[ GE[j] ] + b_A*(AG[j] - min_AG) +
                                b_E[ ED[j] ] + b_S[ SP[j] ] + b_Di[ DI[j] ] +
                                b_Xpu[ XPpu[j] ] + b_Xpr[ XPpr[j] ] , D);
      // using the same predictor for the three latents
    }
    theta_sub ~ multi_normal( m_mult, Rho_theta_sub );
    
    // model
    for( i in 1:N ) {
      v = theta_sub[ IDj[i], IDd[i] ] - b_k[ IDk[i] ];
      p = inv_logit(v);
      y[i] ~ bernoulli(p);
    }
}
generated quantities{
    vector[N] log_lik;
    real v;
    real p;

    // likelihood
    for( i in 1:N ) {
      v = theta_sub[ IDj[i], IDd[i] ] - b_k[ IDk[i] ];
      p = inv_logit(v);
      log_lik[i] = bernoulli_lpmf( y[i] | p);
    }
}

