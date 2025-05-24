# running Bayesian dirichlet-multinomial model of best-worst choice

# pre-setup ====================================================================================
rm(list=ls())
library(here)
library(tidyverse)
library(glue)
library(fs)
library(rstan)
library(bayesplot)


# general setup ====================================================================================
# model settings
which_model <- "dm_hier_1"
debug_model <- F

# paths
model_dir <- here("analysis","bayes",which_model)
stan_file <- path(model_dir,glue("{which_model}.stan"))

# sampler settings
n_chains <- ifelse(debug_model, 1, 4)
n_iter <- ifelse(debug_model, 100, 3000)

# DATA SETUP =================================================================================
d <- here("data","cleaned","bw_all.csv") %>%
  read_csv() %>%
  filter(str_detect(effect,"attraction")) %>%
  arrange(sub_n, block_n, trial_n) %>%
  rowwise() %>%
  mutate(min=which.min(c(a1,a2,a3))) %>%
  ungroup() %>%
  mutate(best=case_when(
    choice_best==min~"d",
    choice_best==1 & h1>w1 ~ "h",
    choice_best==1 & w1>h1 ~ "w",
    choice_best==2 & h2>w2 ~ "h",
    choice_best==2 & w2>h2 ~ "w",
    choice_best==3 & h3>w3 ~ "h",
    choice_best==3 & w3>h3 ~ "w",
  ),
  worst=case_when(
    choice_worst==min~"d",
    choice_worst==1 & h1>w1 ~ "h",
    choice_worst==1 & w1>h1 ~ "w",
    choice_worst==2 & h2>w2 ~ "h",
    choice_worst==2 & w2>h2 ~ "w",
    choice_worst==3 & h3>w3 ~ "h",
    choice_worst==3 & w3>h3 ~ "w",
  ),
  best_att=case_when(
    best=="h" & set=="h"~"t",
    best=="h" & set=="w"~"c",
    best=="w" & set=="w"~"t",
    best=="w" & set=="h"~"c",
    best=="d"~"d"
  ),
  worst_att=case_when(
    worst=="h" & set=="h"~"t",
    worst=="h" & set=="w"~"c",
    worst=="w" & set=="w"~"t",
    worst=="w" & set=="h"~"c",
    worst=="d"~"d"
  )) %>%
  rowwise() %>%
  mutate(middle_att=setdiff(c("t","c","d"),c(best_att,worst_att))) %>%
  ungroup() %>%
  mutate(order=str_c(best_att,middle_att,worst_att)) %>%
  select(sub_n, set, distance,order) 

# tcd tdc ctd cdt dtc dct
d_counts <- d %>%
  group_by(sub_n,set,distance,order) %>%
  summarise(N=n()) %>%
  ungroup() %>%
  pivot_wider(names_from = order, values_from = N, values_fill = 0) %>%
  relocate(c(tcd, tdc, ctd, cdt, dtc, dct),.after=distance)

n_subs <- length(unique(d_counts$sub_n))

# re-numbering subjects to be sequential
subs_key <- tibble(
  sub_n = sort(unique(d_counts$sub_n)),
  sub_n_new = seq(1,n_subs,1)
)
d_counts_clean <- d_counts %>%
  left_join(subs_key) %>%
  relocate(sub_n_new,.after=sub_n)
sub_ns <- unique(d_counts_clean$sub_n_new) # NEW SUBJECT NUMBERS
counts <- array(NA_integer_,dim=c(n_subs,4,2,6))
distance <- unique(d_counts_clean$distance) # 2 5 9 14
set <- unique(d_counts$set) # h w
for(s in sub_ns){
  print(s)
  for(d in 1:length(distance)){
    for(or in 1:length(set)){
      tmp <- d_counts_clean %>%
        filter(sub_n_new==s & distance==distance[d] & set==set[or]) %>%
        select(-c(sub_n,sub_n_new,set,distance)) %>%
        as.matrix()
      counts[s,d,or,] <- tmp
    }
  }
}
stan_data <- list(S=n_subs,
                  D=length(distance),
                  O=length(set),
                  K=6,
                  counts=counts)

# compile model and sample from posterior =================================================================================
m <- stan_model(stan_file)
fit_file <- path(model_dir,glue("{which_model}_fit.RData"))
if(!file_exists(fit_file) | debug_model){
  fit <- sampling(m, data=stan_data,
                  chains=n_chains,
                  iter=n_iter,
                  cores=n_chains)
  if(!debug_model){
    save(fit, file=fit_file)
    fit_summary <- summary(fit, probs=c(.025, .975))
    save(stan_data, d_counts_clean, subs_key, file=path(model_dir,glue("{which_model}_data_for_model.RData")))
    save(fit_summary, file=path(model_dir,glue("{which_model}_fit_summary.RData")))
    diagnostics <- get_sampler_params(fit)
    save(diagnostics, file=path(model_dir,glue("{which_model}_diagnostics.RData")))
  }
  
  # get model predictions ================================================================================================================================================
  p_best <- extract(fit, pars="p_best")$p_best
  p_worst <- extract(fit, pars="p_worst")$p_worst
  counts_rank_rep <- extract(fit, pars="counts_rank_rep")$counts_rank_rep # ppc 
  
  save(p_best,file=path(model_dir,glue("{which_model}_p_best.RData")))
  save(p_worst,file=path(model_dir,glue("{which_model}_p_worst.RData")))
  save(counts_rank_rep,file=path(model_dir,glue("{which_model}_counts_rank_rep.RData")))
}



