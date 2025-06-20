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
which_model <- "maxdiff_1"
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

d_counts <- d %>%
  group_by(sub_n,set,distance,order) %>%
  summarise(N=n()) %>%
  ungroup() %>%
  pivot_wider(names_from = order, values_from = N, values_fill = 0) 
n_subs <- length(unique(d_counts$sub_n))

# re-numbering subjects to be sequential
subs_key <- tibble(
  sub_n = sort(unique(d_counts$sub_n)),
  sub_n_new = seq(1,n_subs,1)
)
d_counts_clean <- d_counts %>%
  left_join(subs_key) %>%
  relocate(sub_n_new,.after=sub_n)

N <- nrow(d_counts_clean)
sub_ns <- d_counts_clean$sub_n_new
t_w <- d_w <- if_else(d_counts_clean$set=="w", 1, 0)
c_w <- if_else(d_counts_clean$set=="w", 0, 1)
distance2 <- if_else(d_counts_clean$distance==2,1,0)
distance5 <- if_else(d_counts_clean$distance==5,1,0)
distance9 <- if_else(d_counts_clean$distance==9,1,0)
distance14 <- if_else(d_counts_clean$distance==14,1,0)
# IN ORDER OF TCD, TDC, CTD, CDT, DTC, DCT
counts <- cbind(d_counts_clean$tcd,
                d_counts_clean$tdc,
                d_counts_clean$ctd,
                d_counts_clean$cdt,
                d_counts_clean$dtc,
                d_counts_clean$dct)
stan_data <- list(N=N,
                  n_subs=n_subs,
                  sub_ns=sub_ns,
                  t_w=t_w,
                  d_w=d_w,
                  c_w=c_w,
                  distance2=distance2,
                  distance5=distance5,
                  distance9=distance9,
                  distance14=distance14,
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
  
  params <- c("lp__",
              "b_0",
              "b_competitor",
              "b_decoy",
              "b_competitor_best",
              "b_decoy_best",
              "b_competitor_worst",
              "b_decoy_worst")
  for(par in params){
    try({
      p <- mcmc_trace(fit,par)
      p
      ggsave(p, filename=path(model_dir,glue("{which_model}_{par}_trace.jpeg")),width=5,height=4)
      rm(p)
    })
    try({
      p <- mcmc_hist(fit,par)
      p
      ggsave(p, filename=path(model_dir,glue("{which_model}_{par}_hist.jpeg")),width=5,height=4)
      rm(p)
    })
    try({
      p <- mcmc_intervals(fit,par,prob=.95,prob_outer=1,point_est="mean")
      p
      ggsave(p, filename=path(model_dir,glue("{which_model}_{par}_intervals.jpeg")),width=5,height=4)
      rm(p)
    })
    try({
      p <- mcmc_dens(fit,par)
      p
      ggsave(p, filename=path(model_dir,glue("{which_model}_{par}_dens.jpeg")),width=5,height=4)
      rm(p)
    })
  }
  
  try({
    p <- mcmc_trace(fit,regex_pars = "sigma")
    p
    ggsave(p, filename=path(model_dir,glue("{which_model}_sigma_trace.jpeg")),width=5,height=4)
    rm(p)
  })
  
  try({
    p <- mcmc_hist(fit,regex_pars = "sigma")
    p
    ggsave(p, filename=path(model_dir,glue("{which_model}_sigma_hist.jpeg")),width=5,height=4)
    rm(p)
  })
  
  try({
    p <- mcmc_trace(fit,c("b_distance_2",
                          "b_distance_9",
                          "b_distance_14"))
    p
    ggsave(p, filename=path(model_dir,glue("{which_model}_distance_trace.jpeg")),width=5,height=4)
    rm(p)
  })
  
  try({
    p <- mcmc_hist(fit,c("b_distance_2",
                         "b_distance_9",
                         "b_distance_14"))
    p
    ggsave(p, filename=path(model_dir,glue("{which_model}_distance_hist.jpeg")),width=5,height=4)
    rm(p)
  })
  
  # get model predictions ================================================================================================================================================
  p_best <- extract(fit, pars="p_best")$p_best
  p_worst <- extract(fit, pars="p_worst")$p_worst
  counts_rank_rep <- extract(fit, pars="counts_rank_rep")$counts_rank_rep
  
  save(p_best,file=path(model_dir,glue("{which_model}_p_best.RData")))
  save(p_worst,file=path(model_dir,glue("{which_model}_p_worst.RData")))
  save(counts_rank_rep,file=path(model_dir,glue("{which_model}_counts_rank_rep.RData")))
}else{
  load(fit_file)
  U <- extract(fit,pars="U")$U
  distances <- case_when(d_counts_clean$distance==2~1,
                         d_counts_clean$distance==5~2,
                        d_counts_clean$distance==9~3,
                        d_counts_clean$distance==14~4)
  U_summary <- tibble()
  for(i in 1:4){
    print(i)
    for(j in 1:3){
      tmp <- rowMeans(U[,distances==distances[i],j])
      U_summary <- bind_rows(U_summary,
                             tibble(
                               option=c("t","c","d")[j],
                               distance=c(2,5,9,14)[i],
                               m=mean(tmp),
                               l=HDInterval::hdi(tmp)[1],
                               u=HDInterval::hdi(tmp)[2]
                             ))
    }
  }
  U_summary %>%
    ggplot(aes(distance,m,fill=option))+
    geom_col(position="dodge",width=1.2)+
    geom_errorbar(aes(ymin=l,ymax=u),position=position_dodge(width=1.2),width=.05)+
    geom_hline(yintercept=0,alpha=.5)+
    scale_x_continuous(breaks=c(2,5,9,14),limits=c(0,16),labels=c("2%","5%","9%","14%"))+
    # scale_y_continuous(limits=c(-.1,.05))+
    # scale_shape_manual(values=c(1,4),name="")+
    ggsci::scale_fill_startrek()+
    labs(x="tdd",y="mean U")+
    ggthemes::theme_few()
  ggsave(filename=path(model_dir,glue("{which_model}_U_mean_posterior.jpeg")),width=5,height=4)
  

}



