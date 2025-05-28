# pre-setup ====================================================================================
rm(list=ls())
library(here)
library(tidyverse)
library(glue)
library(HDInterval)
library(fs)
library(rstan)
library(bayesplot)


# general setup ====================================================================================
# model settings
which_model <- "maxdiff_1"

# paths
model_dir <- here("analysis","bayes",which_model)
stan_file <- path(model_dir,glue("{which_model}.stan"))

load(path(model_dir,glue("{which_model}_p_best.RData")))
load(path(model_dir,glue("{which_model}_p_worst.RData")))
load(path(model_dir,glue("{which_model}_data_for_model.RData")))

get_model_mean_preds <- function(preds, data, type){
  # browser()
  sub_ns_new <- unique(d_counts_clean$sub_n_new)
  distances <- unique(d_counts_clean$distance)
  n_dists <- length(distances)
  mean_preds <- tibble()
  i <- 1
  for(d in distances){
    print(d)
    t_tmp <- apply(preds[,d_counts_clean$distance==d,1],1,mean)
    c_tmp <- apply(preds[,d_counts_clean$distance==d,2],1,mean)
    d_tmp <- apply(preds[,d_counts_clean$distance==d,3],1,mean)
    mean_preds <- bind_rows(
      mean_preds,
      tibble(
        type=type,
        distance=d,
        choice=c("t","c","d"),
        m=c(mean(t_tmp),
            mean(c_tmp),
            mean(d_tmp)),
        lower=c(hdi(t_tmp)[1],
                    hdi(c_tmp)[1],
                    hdi(d_tmp)[1]),
        upper=c(hdi(t_tmp)[2],
                    hdi(c_tmp)[2],
                    hdi(d_tmp)[2])
      )
    )
    i <- i+1
  }
  return(mean_preds)
}

preds_all <- bind_rows(get_model_mean_preds(p_best,d_counts_clean,"best"),
                       get_model_mean_preds(p_worst,d_counts_clean,"worst")) %>%
  pivot_wider(names_from = type,
              values_from = c(m,lower,upper)) %>%
  mutate(source="model")

analyze_data <- function(){
  # read in data, only include critical trials, recode bw cond to be more understandable
  d <- here("data","cleaned","bw_all.csv") %>%
    read_csv() %>%
    mutate(bw_cond=recode(bw_cond,"bw"="best-worst","wb"="worst-best")) %>%
    filter(str_detect(effect,"attraction")) %>%
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
  d %>%
    pivot_longer(contains("att"),names_to = "type",values_to = "option") %>%
    mutate(type=str_remove(type,"_att")) %>%
    group_by(sub_n,distance,type,option) %>%
    summarise(n=n()) %>%
    group_by(sub_n,distance,type) %>%
    mutate(prop=n/sum(n)) %>%
    ungroup() %>%
    select(-n) %>%
    left_join(distinct(d,sub_n,bw_cond)) %>%
    group_by(distance,type,option) %>%
    summarise(m_prop=mean(prop),
              n=n(),
              se=sd(prop)/sqrt(n),
              lower=m_prop-qt(.975,n-1)*se,
              upper=m_prop-qt(.975,n-1)*se) %>%
    ungroup() %>%
    select(-c(n,se)) %>%
    pivot_wider(names_from = type,
                values_from = c(m_prop,lower,upper)) %>%
    rename(choice=option,
           m_best=m_prop_best,
           m_worst=m_prop_worst) %>%
    mutate(source="data")
}

data_preds_all <- analyze_data() %>%
  bind_rows(preds_all)

# 
p <- data_preds_all %>%
  ggplot(aes(m_worst,m_best,col=choice,shape=source))+
  geom_point(alpha=.5)+
  coord_fixed(xlim=c(0,.75),ylim=c(0,.75))+
  scale_shape_manual(name="",
                     values = c(1,4))+
  facet_grid(distance~.)+
  ggsci::scale_color_startrek(name="")+
  labs(x="mean p(worst)",y="mean p(best)")+
  ggthemes::theme_few()+
  theme(text = element_text(size=19))
p
ggsave(p, filename = path(model_dir,glue("{which_model}_means_model_v_data.jpeg")),width=6,height=8)  
