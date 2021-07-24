
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
    // individuals
    real<lower=0> s_ind;
    vector[J] zind_j;
    real<lower=0> s_dim;
    vector[D] zdim_d;
    
    // items and text
    real<lower=0> s_item;   
    vector[K] zitem_k;
    real<lower=0> s_text;
    vector[L] ztext_l;

    // betas
    real a;
    real b_G[2];
    real b_A;
    real b_E[3];
    real b_X[4];
}
transformed parameters{
    vector[J] m_ind;
    vector[J] ind_j;
    vector[D] dim_d;
    vector[K] item_k;
    vector[L] text_l;
    
    // individuals
    for(j in 1:J){
    m_ind[j] = a + b_G[ G[j] ] +
                b_A * ( A[j] - min(A) ) +
                b_E[ E[j] ] +
                b_X[ X[j] ];
    }  
    ind_j = m_ind + s_ind * zind_j;
    dim_d = s_dim * zdim_d;
    
    // items and texts
    item_k = s_item * zitem_k;
    text_l = s_text * ztext_l;
}
model{
    // declare
    real v;
    real p; 
    
    // abilities
    s_ind ~ exponential(2);
    zind_j ~ normal(0, 1);
    s_dim ~ exponential(2);
    zdim_d ~ normal(0, 1);
    
    a ~ normal(0, 0.5);
    b_G ~ normal(0, 0.5);
    b_A ~ normal(0, 0.5);
    b_E ~ normal(0, 1);
    b_X ~ normal(0, 0.5);
    
    // items and texts
    s_item ~ exponential(2);
    zitem_k ~ normal(0, 1);
    s_text ~ exponential(2);
    ztext_l ~ normal(0, 1);

    // model
    for( i in 1:N ) {
      v = ( ind_j[ IDj[i] ] + dim_d[ IDd[i] ] ) -
          ( item_k[ IDk[i] ] + text_l[ IDl[i] ] );
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
//#       v = ( ind_j[ IDj[i] ] + dim_d[ IDd[i] ] ) -
//#           ( item_k[ IDk[i] ] + text_l[ IDl[i] ] );
//#       p = inv_logit(v);
//#       log_lik[i] = bernoulli_lpmf( y[i] | p);
//#     }
//# }

