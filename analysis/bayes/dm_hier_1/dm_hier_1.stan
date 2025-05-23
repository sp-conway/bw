data {
  int S; // Number of subjects
  int D; // Number of distances
  int O; // target decoy orientaiton
  int K; // number of options (6)
  array[S,D,O,K] int counts; // counts of all rankings, IN ORDER OF TCD, TDC, CTD, CDT, DTC, DCT
}

parameters {
  array[S,D,O] simplex[K] theta;
  array[D,O] simplex[K] alpha_prob;
  array[D,O] real concentration; 
}

transformed parameters{
  array[D,O] vector<lower=0> [K] alpha;
  for(d in 1:D){
    for(o in 1:O){
      alpha[d,o]=alpha_prob[d,o]*concentration[d,o];
    }
  }
}

model{
  for(d in 1:D){
    for(o in 1:O){
      concentration[d,o] ~ uniform(0,500);
      alpha_prob[d,o] ~ dirichlet([1,1,1,1,1,1]);
    }
  }
  for(s in 1:S){
    for(d in 1:D){
      for(o in 1:O){
        theta[s,d,o] ~ dirichlet(alpha[d,o]);
      }
    }
  }
}

generated quantities {
  array[S,D,O] vector<lower=0,upper=1> [3] p_best;
  array[S,D,O] vector<lower=0,upper=1> [3] p_worst;
  array[S,D,O,6] int counts_rank_rep;
  for(s in 1:S){
    for(d in 1:D){
      for(o in 1:O){
        p_best[s,d,o][1]=theta[s,d,o][1]+theta[s,d,o][2];
        p_best[s,d,o][2]=theta[s,d,o][3]+theta[s,d,o][4];
        p_best[s,d,o][3]=theta[s,d,o][5]+theta[s,d,o][6];
        
        p_worst[s,d,o][1]=theta[s,d,o][1]+theta[s,d,o][2];
        p_worst[s,d,o][2]=theta[s,d,o][3]+theta[s,d,o][4];
        p_worst[s,d,o][3]=theta[s,d,o][5]+theta[s,d,o][6];
        
        counts_rank_rep[s,d,o,] = multinomial_rng(to_vector(theta[s,d,o, ]), sum(counts[s,d,o,]));
      }
    }
  }
}
