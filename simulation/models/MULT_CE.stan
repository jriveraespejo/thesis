
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
    vector[J] ind_j;
    real<lower=0> s_dim;
    vector[D] dim_d;
    
    // items and text
    real<lower=0> s_item;   
    vector[K] item_k;
    real<lower=0> s_text;
    vector[L] text_l;

    // betas
    real a;
    real b_G[2];
    real b_A;
    real b_E[3];
    real b_X[4];
}
model{
    // declare
    real v;
    real p; 
    
    // abilities
    s_ind ~ exponential(2);
    ind_j ~ normal(0, s_ind);
    s_dim ~ exponential(2);
    dim_d ~ normal(0, s_dim);
    
    a ~ normal(0, 0.5);
    b_G ~ normal(0, 0.5);
    b_A ~ normal(0, 0.5);
    b_E ~ normal(0, 1);
    b_X ~ normal(0, 0.5);
    
    // items and texts
    s_item ~ exponential(2);
    item_k ~ normal(0, s_item);
    s_text ~ exponential(2);
    text_l ~ normal(0, s_text);

    // model
    for( i in 1:N ) {
      v = ( ind_j[ IDj[i] ] + dim_d[ IDd[i] ] ) -
          ( item_k[ IDk[i] ] + text_l[ IDl[i] ] ) +
          a + b_G[ GE[i] ] +
                b_A * ( AG[i] - min(A) ) +
                b_E[ ED[i] ] +
                b_X[ XP[i] ];
      p = inv_logit(v);
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
//#           ( item_k[ IDk[i] ] + text_l[ IDl[i] ] ) +
//#           a + b_G[ GE[i] ] +
//#                 b_A * ( AG[i] - min(A) ) +
//#                 b_E[ ED[i] ] +
//#                 b_X[ XP[i] ];
//#       p = inv_logit(v);
//#       log_lik[i] = bernoulli_lpmf( y[i] | p);
//#     }
//# }

