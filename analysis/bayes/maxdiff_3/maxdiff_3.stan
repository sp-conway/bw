data {
  int<lower=0> n_subs; // number of subjects
  int N; // total data points
  int sub_ns[N]; // vector of subject indices NOTE THAT THESE WERE CHANGED TO BE SEQUENTIAL, SEE KEY
  vector[N] t_w;
  vector[N] c_w;
  vector[N] d_w;
  vector[N] distance5;
  vector[N] distance9;
  vector[N] distance14;
  array[N,6] int counts; // counts of all rankings, IN ORDER OF TCD, TDC, CTD, CDT, DTC, DCT
}

parameters {
  real b_w; // effect of wide rectangle
  real b_distance5; // effect of distance = 5
  real b_distance9; // effect of distance = 9
  real b_distance14; // effect of distance = 14
  real b_competitor; // effect of competitor  
  real b_target;
  real b_decoy; // effect of decoy  
  real<lower=0> b_target_s_sigma;
  real<lower=0> b_competitor_s_sigma;
  real<lower=0> b_decoy_s_sigma;
  vector[n_subs] b_target_s; // 
  vector[n_subs] b_competitor_s; // 
  vector[n_subs] b_decoy_s; // 
}

transformed parameters{
  matrix [N,3] U;  // all utilities
  array[N] matrix[3, 3] Ud; // diffs
  array[N] matrix[3, 3] p; // probability of all rankings
  matrix[N,6] p_rank; // probability of all rankings removing 0s and put in order for model
  for(i in 1:N){
    
    // Utility of target option
    U[i,1] = b_target+
            b_target_s[sub_ns[i]]+
            b_w*t_w[i] + 
              b_distance5*distance5[i] +
              b_distance9*distance9[i] + 
              b_distance14*distance14[i];
              
    // Utility of competitor option          
    U[i,2] = b_competitor+
              b_competitor_s[sub_ns[i]]+
              b_w*c_w[i] + 
              b_distance5*distance5[i] +
              b_distance9*distance9[i] + 
              b_distance14*distance14[i];
    
    // Utility of decoy option
    U[i,3] = b_decoy+
              b_decoy_s[sub_ns[i]]+
              b_w*d_w[i];
              
    // compute utility differences
    for(j in 1:3){
      for(k in 1:3){
        if(j==k){
          Ud[i][j,k] = 0; // can't pick same option for best and worst
        }else{
          Ud[i][j,k] = U[i,j]-U[i,k];
        }
      }
    }
    
    // need to sum up denominator to ensure we don't include cases when j=k
    // likely slows down sampler a lot but oh well
    real sum_exp_Ud = 0;
    for(j in 1:3){
      for(k in 1:3){
        if(j!=k){
          sum_exp_Ud += exp(Ud[i][j,k]);
        }
      }
    }
    
    // compute probabilities for each ranking
    for(j in 1:3){
      for(k in 1:3){
        if(j==k){
          p[i][j,k]=0;
        }else{
          p[i][j,k]=exp(Ud[i][j,k])/sum_exp_Ud;
        }
      }
    }
    
    
    // flatten probabilities and assign CORRECTLY to rankings
    p_rank[i, 1] = p[i][1, 2]; 
    p_rank[i, 2] = p[i][1, 3];
    p_rank[i, 3] = p[i][2, 1]; 
    p_rank[i, 4] = p[i][2, 3]; 
    p_rank[i, 5] = p[i][3, 1];
    p_rank[i, 6] = p[i][3, 2];
  }
}

model {
  b_w ~ normal(0,5);
  b_distance5 ~ normal(0,5); // effect of distance = 5
  b_distance9 ~ normal(0,5); // effect of distance = 9
  b_distance14 ~ normal(0,5); // effect of distance = 14
  b_target ~ normal(0,5); //  
  b_competitor ~ normal(0,5); //  
  b_decoy ~ normal(0,5); //  
  b_target_s_sigma ~ cauchy(0,2.5);
  b_target_s ~ normal(0, b_target_s_sigma); // 
  
  b_competitor_s_sigma ~ cauchy(0,2.5);
  b_competitor_s ~ normal(0, b_competitor_s_sigma); // 
  
  b_decoy_s_sigma ~ cauchy(0,2.5);
  b_decoy_s ~ normal(0, b_decoy_s_sigma); // 
  for(m in 1:N){
   counts[m,]~ multinomial(to_vector(p_rank[m, ]));
  }
}

// get best and worst probs
// also do ppc
generated quantities {
  matrix[N,3] p_best;
  matrix[N,3] p_worst;
  array[N,6] int counts_rank_rep;
  for(m in 1:N){
    p_best[m,1]=p_rank[m,1]+p_rank[m,2];
    p_best[m,2]=p_rank[m,3]+p_rank[m,4];
    p_best[m,3]=p_rank[m,5]+p_rank[m,6];
    
    p_worst[m,1]=p_rank[m,3]+p_rank[m,5]; 
    p_worst[m,2]=p_rank[m,1]+p_rank[m,6];
    p_worst[m,3]=p_rank[m,2]+p_rank[m,4];
    
    counts_rank_rep[m,]=multinomial_rng(to_vector(p_rank[m, ]), sum(counts[m,]));
  }
}
