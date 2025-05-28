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
which_model <- "dm_hier_2"
debug_model <- F

# paths
model_dir <- here("analysis","bayes",which_model)
stan_file <- path(model_dir,glue("{which_model}.stan"))

# sampler settings
n_chains <- ifelse(debug_model, 1, 4)
n_iter <- ifelse(debug_model, 100, 3000)

# DATA SETUP =================================================================================
dat <- here("data","cleaned","bw_all.csv") %>%
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
  ))
dat1 <- dat %>%
  rowwise() %>%
  mutate(middle_att=setdiff(c("t","c","d"),c(best_att,worst_att))) %>%
  ungroup() %>%
  mutate(order=str_c(best_att,middle_att,worst_att)) %>%
  select(sub_n, set, distance,order) 

# tcd tdc ctd cdt dtc dct
d_counts <- dat1 %>%
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
                  K=6, # n options
                  counts=counts)

# compile model and sample from posterior =================================================================================
fit_file <- path(model_dir,glue("{which_model}_fit.RData"))
if(!file_exists(fit_file) | debug_model){
  m <- stan_model(stan_file)
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
}else{
  load(fit_file)
}

fit_summary <- summary(fit, probs=c(.025, .975))$summary
diagnostics <- get_sampler_params(fit)
p_best <- extract(fit, pars="p_best")$p_best
p_worst <- extract(fit, pars="p_worst")$p_worst
counts_rank_rep <- extract(fit, pars="counts_rank_rep")$counts_rank_rep # ppc 

color_scheme_set("mix-blue-pink")
mcmc_trace(fit, "lp__")
ggsave(filename=path(model_dir,"lp_trace.jpeg"),width=5,height=4)

mcmc_trace(fit,regex_pars = "alpha\\[1,")
ggsave(filename=path(model_dir,"alpha1_trace.jpeg"),width=8,height=8)

mcmc_trace(fit,regex_pars = "alpha\\[2,")
ggsave(filename=path(model_dir,"alpha2_trace.jpeg"),width=8,height=8)

mcmc_trace(fit,regex_pars = "alpha\\[3,")
ggsave(filename=path(model_dir,"alpha3_trace.jpeg"),width=8,height=8)

mcmc_trace(fit,regex_pars = "alpha\\[4,")
ggsave(filename=path(model_dir,"alpha4_trace.jpeg"),width=8,height=8)


mcmc_trace(fit,regex_pars = "alpha_prob\\[1,")
ggsave(filename=path(model_dir,"alphaprob1_trace.jpeg"),width=8,height=8)

mcmc_trace(fit,regex_pars = "alpha_prob\\[2,")
ggsave(filename=path(model_dir,"alphaprob2_trace.jpeg"),width=8,height=8)

mcmc_trace(fit,regex_pars = "alpha_prob\\[3,")
ggsave(filename=path(model_dir,"alphaprob3_trace.jpeg"),width=8,height=8)

mcmc_trace(fit,regex_pars = "alpha_prob\\[4,")
ggsave(filename=path(model_dir,"alphaprob4_trace.jpeg"),width=8,height=8)

mcmc_trace(fit,regex_pars = "concentration")
ggsave(filename=path(model_dir,"conc_trace.jpeg"),width=5,height=4)

summary(fit_summary[,"Rhat"])

# Getting group means and subject-level means ===============================================================
# mean p_best and p_worst over subjects for each iteration
# Resulting dimensions: [iterations, D, K]
p_best_iter_means <- apply(p_best, c(1, 3, 5), mean)
p_worst_iter_means <- apply(p_worst, c(1, 3, 5), mean)


mean_preds <- tibble()
for(D in 1:4){
    for(K in 1:3){
      mean_preds <- bind_rows(mean_preds,
                tibble(distance=distance[D],
                       option=c("t","c","d")[K],
                       p_best_m=mean(p_best_iter_means[,D,K]),
                       p_best_m_l=HDInterval::hdi(p_best_iter_means[,D,K])[1],
                       p_best_m_u=HDInterval::hdi(p_best_iter_means[,D,K])[2],
                       p_worst_m=mean(p_worst_iter_means[,D,K]),
                       p_worst_m_l=HDInterval::hdi(p_worst_iter_means[,D,K])[1],
                       p_worst_m_u=HDInterval::hdi(p_worst_iter_means[,D,K])[2])
      )
    }
}
save(mean_preds,file=path(model_dir,"model_mean_preds.RData"))

p_best_iter_sub_means <- apply(p_best, c(1, 2, 3, 5), mean)
p_worst_iter_sub_means <- apply(p_worst, c(1,2, 3, 5), mean)
sub_preds <- tibble()
for(S in 1:n_subs){
  print(S)
  for(D in 1:4){
    for(K in 1:3){
      sub_preds <- bind_rows(sub_preds,
                              tibble(sub_n=sub_ns[S],
                                     distance=distance[D],
                                     option=c("t","c","d")[K],
                                     p_best_m=mean(p_best_iter_sub_means[,S,D,K]),
                                     p_best_m_l=HDInterval::hdi(p_best_iter_sub_means[,S,D,K])[1],
                                     p_best_m_u=HDInterval::hdi(p_best_iter_sub_means[,S,D,K])[2],
                                     p_worst_m=mean(p_worst_iter_sub_means[,S,D,K]),
                                     p_worst_m_l=HDInterval::hdi(p_worst_iter_sub_means[,S,D,K])[1],
                                     p_worst_m_u=HDInterval::hdi(p_worst_iter_sub_means[,S,D,K])[2])
      )
    }
  }
}
save(sub_preds,file=path(model_dir,"model_sub_mean_preds.RData"))

# Compare model predictions to data ======================================================================

sub_data <- dat %>%
  select(-c(h1,w1,h2,w2,h3,w3,a1,a2,a3,rt_best,rt_worst,choice_best,choice_worst,min,best,worst)) %>%
  pivot_longer(contains("att"),names_to = "type",values_to = "option") %>%
  mutate(type=str_remove(type,"_att")) %>%
  group_by(distance,type,option,sub_n) %>%
  summarise(n=n()) %>%
  group_by(distance,type,sub_n) %>%
  mutate(p=n/sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  left_join(subs_key) %>%
  select(-sub_n) %>%
  rename(sub_n=sub_n_new) 

sub_data_to_model <- sub_data %>%
  pivot_wider(names_from = type, 
              values_from = p,
              names_glue = "p_{type}_m") %>%
  mutate(source="data") %>%
  bind_rows(mutate(sub_preds,source="model"))
sub_data_to_model %>%
  pivot_wider(names_from = source,
              values_from = contains("p_")) %>%
  ggplot(aes(p_best_m_data,p_best_m_model,col=option))+
  geom_point(aes(alpha=.25),shape=".")+
  labs(x="p(best) data",y="p(best) model")+
  # geom_errorbar(aes(ymin=p_best_m_l_model,ymax=p_best_m_u_model))+
  coord_fixed(xlim=c(0,1),ylim=c(0,1))+
  facet_grid(distance~.)+
  ggsci::scale_color_startrek()+
  ggthemes::theme_few()
ggsave(filename=path(model_dir,"pbest_sub_model_data_preds.jpeg"),width=5,height = 6)

sub_data_to_model %>%
  pivot_wider(names_from = source,
              values_from = contains("p_")) %>%
  ggplot(aes(p_worst_m_data,p_worst_m_model,col=option))+
  geom_point(aes(alpha=.25),shape=".")+
  labs(x="p(worst) data",y="p(worst) model")+
  # geom_errorbar(aes(ymin=p_worst_m_l_model,ymax=p_worst_m_u_model))+
  coord_fixed(xlim=c(0,1),ylim=c(0,1))+
  facet_grid(distance~.)+
  ggsci::scale_color_startrek()+
  ggthemes::theme_few()
ggsave(filename=path(model_dir,"pworst_sub_model_data_preds.jpeg"),width=5,height = 6)

sub_data %>%
  group_by(distance,option,type) %>%
  summarise(m=mean(p)) %>%
  ungroup() %>%
  mutate(source="data") %>%
  bind_rows(
    mutate(mean_preds,source="model") %>%
      pivot_longer(contains("p_")) %>%
      mutate(name=str_remove(name,"p_")) %>%
      mutate(name=str_replace_all(name,c("m_l"="l","m_u"="u"))) %>%
      separate(name,into=c("type","stat")) %>%
      pivot_wider(names_from = stat, values_from=value)
  ) %>%
  pivot_wider(names_from = type,
              values_from = c(m,l,u)) %>%
  ggplot(aes(m_worst,m_best,col=option,shape=source))+
  geom_point(alpha=.8)+
  coord_fixed(xlim=c(0,1),ylim=c(0,1))+
  labs(x="mean p(worst)",y="mean p(best)")+
  ggsci::scale_color_startrek()+
  scale_shape_manual(values=c(1,4))+
  facet_grid(distance~.)+
  ggthemes::theme_few()
ggsave(filename=path(model_dir,"mean_probs_model_data_preds.jpeg"),width=5,height = 6)

theta <- extract(fit,"theta")$theta      
mean(theta[,,,1,1])

d_props_clean <- d_counts_clean %>%
  rowwise() %>%
  mutate(N=sum(c(tcd,tdc,ctd,cdt,dtc,dct))) %>%
  ungroup() %>%
  mutate(tcd=tcd/N,
         tdc=tdc/N,
         ctd=ctd/N,
         cdt=cdt/N,
         dtc=dtc/N,
         dct=dct/N) %>%
  select(-N) %>%
  pivot_longer(-c(sub_n,sub_n_new,set,distance),names_to = "order",values_to = "prop")


mrank_model <- c()
for(i in 1:length(distance)){
  for(j in 1:6){
    mrank_model <- c(mrank_model, mean(theta[,,i,,j]))
  }
}

d_props_clean %>%
  group_by(distance,order) %>%
  summarise(m=mean(prop)) %>%
  ungroup() %>%
  mutate(m_model=mrank_model) %>%
  ggplot(aes(m,m_model,col=order))+
  geom_point()+
  labs(x="m data",y="m model")+
  ggsci::scale_color_startrek()+
  scale_shape_manual(values=c(1,4))+
  coord_fixed(xlim=c(0,1),ylim=c(0,1))+
  facet_grid(distance~.)+
  ggthemes::theme_few()
ggsave(filename=path(model_dir,"mean_rank_probs_model_data_preds.jpeg"),width=5,height = 6)
