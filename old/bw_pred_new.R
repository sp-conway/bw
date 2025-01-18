# Setup =================================================================================
rm(list=ls())
library(tidyverse)
library(here)
library(glue)
library(posterior)
library(bayesplot)
library(fs)
library(qs)
library(mvtnorm)

# results directory
results_dir <- here("new_bayes_modeling_results")
preds_dir <- here("bw_sim_from_new_bayes")
# CONTROL PARAMS
N <- 5000

# CONSTANTS
n_chain <- 5 # CANNOT CHANGE THESE - NUMBER OF CHAINS USED WHEN FITTING MODEL
n_iter <- 4000 # CANNOT CHANGE THESE - NUMBER OF CHAINS USED WHEN FITTING MODEL

# number of total draws
n_draws <- n_chain*n_iter

# load in modeling results ====================================================================
dat <- path(results_dir,"dat_for_model.qs") %>%
  qread() # # data 
d_clean <- dat[[1]]
load(path(results_dir,"mu_cv_collapsed.RData"))

# function for simulation ===================================================================
sim_choice <- function(N,mu,cv){
  samp <- mvtnorm::rmvnorm(N,mu,cv)
  mx <- apply(samp, 1, which.max)
  mn <- apply(samp, 1, which.min)
  pbest <- c(sum(mx==1),
             sum(mx==2),
             sum(mx==3))/N
  pworst <- c(sum(mn==1),
              sum(mn==2),
              sum(mn==3))/N
  return(cbind(pbest,pworst))
}

# get unique trial type indices per subject =======================================================
# grouping data by all trial types
# then getting an index for each
# removing repeats from mu later for more efficient simulations
# don't need to simualte the same trial type for each subject more than once
d_clean_w_ind <- d_clean %>%
  group_by(sub_n_new,set,distance,diag) %>% 
  mutate(iter=1:n()) %>%
  ungroup() %>%
  mutate(is_unique=iter==1)

# INDICES IN DATA
dat_ind <- which(d_clean_w_ind$is_unique)

# CLEAN DF FOR LATER
d_clean_unique <- filter(d_clean_w_ind, is_unique)

# NUMBER OF DATA POINTS FOR MU
n_dat <- length(dat_ind)


# do simulation !!!! ====================================================================================================
if(!file_exists(path(preds_dir,"bw_preds.RData"))){
  p <- array(NA, dim = c(n_dat,3,2))
  for(i in 1:n_dat){
    cat(i,"/",n_dat,"\n")
    p[i,,] <- sim_choice(N, mu_collapsed[i,], cv_avgd)
  }
  
  tryCatch(qsave(list(p,d_clean_unique),file=path(preds_dir,"bw_preds.qs")),error=function(e) print(e))
  tryCatch(save(p,d_clean_unique,file=path(preds_dir,"bw_preds.RData")),
           error=function(e) print(e))
}else{
  load(path(preds_dir,"bw_preds.RData"))
  if(!file_exists(path(preds_dir,"bw_preds_clean.RData"))){
    # (try) and clean results ====================================================================================================
    try({
      res_df <- bind_cols(
        disp_cond=d_clean_unique$disp_cond,
        sub_n=d_clean_unique$sub_n,
        set=d_clean_unique$set,
        distance=d_clean_unique$distance,
        diag=d_clean_unique$diag,
        target_best=p[,1,1],
        competitor_best=p[,2,1],
        decoy_best=p[,3,1],
        target_worst=p[,1,2],
        competitor_worst=p[,2,2],
        decoy_worst=p[,3,2]) %>%
        pivot_longer(c(contains("best"),contains("worst")),values_to = "prop") %>%
        separate(name,into = c("choice","type")) %>%
        pivot_wider(names_from = type,
                    values_from = prop)
      tryCatch(qsave(res_df,file=path(preds_dir,"bw_preds_clean.qs")),
               error=function(e) print(e))
      tryCatch(save(res_df,file=path(preds_dir,"bw_preds_clean.RData")),
               error=function(e) print(e))
      tryCatch(write_csv(res_df,path(preds_dir,"bw_preds_clean.csv")),
               error=function(e) print(e))
    })
  }else{
    load(path(preds_dir,"bw_preds_clean.RData"))
  }
}

analyze <- function(d,groups){
  p <- d %>%
    pivot_longer(c(best,worst),names_to = "type",values_to = "prop") %>%
    group_by(sub_n,type,choice,!!!groups) %>%
    summarise(prop=mean(prop)) %>%
    group_by(type,choice,!!!groups) %>%
    summarise(m=mean(prop)) %>%
    ungroup() %>%
    pivot_wider(names_from = type,
                values_from = m) %>% 
    mutate(choice=str_sub(choice,1,1)) 
  return(p)
}
res_df %>%
  analyze(groups=vars(set,distance)) %>%
  ggplot(aes(worst,best))+
  geom_text(aes(label=choice), col="darkblue",alpha=.5, size=2.8)+
  coord_fixed(xlim=c(0,.6),ylim=c(0,.6))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +  # Round x-axis labels to 1 decimal
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+
  labs(x="p(worst)",y="p(best)",title = "model predictions")+
  ggthemes::theme_few()+
  facet_grid(as.factor(as.numeric(distance))~set)+
  theme(text = element_text(size=14),plot.title = element_text(hjust=0.5))
ggsave(filename = path(preds_dir,"bw_preds_by_set_dist.jpeg"),width=4,height=5)
res_df %>%
  analyze(groups=vars(distance)) %>%
  ggplot(aes(worst,best))+
  geom_text(aes(label=choice), col="darkblue",alpha=.5, size=2.8)+
  coord_fixed(xlim=c(0,.6),ylim=c(0,.6))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +  # Round x-axis labels to 1 decimal
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))  +
  labs(x="p(worst)",y="p(best)",title = "model predictions")+
  ggthemes::theme_few()+
  facet_grid(as.factor(as.numeric(distance))~.)+
  theme(text = element_text(size=14),plot.title = element_text(hjust=0.5))
ggsave(filename = path(preds_dir,"bw_preds_by_dist.jpeg"),width=4,height=5)
