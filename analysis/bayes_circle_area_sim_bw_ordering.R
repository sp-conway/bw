# simulating choice predictions for bw choice, conditioned on best fitting mean params from circle exp
# ONLY simulating triangle condition
# NOW ONLY INCLUDING sigma_constant_comp_effect
# removed other models from thsi project

rm(list=ls())
library(tidyverse)
library(here)
library(glue)
library(mvtnorm)
library(fs)

N <- 1e7
outl <- "no_outliers"
which_model <- "sigma_constant_comp_effect"
const_tc <- T
eq_cors <- F
vary_decoy_means <- T

load_and_run_model <- function(N,cond,const_tc=F,eq_cors=F,outl="no_outliers",vary_decoy_means=F){
  print(cond)
  dir <- here("analysis","sim_from_bayes_circle_area","bayes_circle_area",which_model,cond,outl)
  load(path(dir,"mu_avgs.RData"))
  load(path(dir,"fit_summary.RData"))
  s <- fit_summary %>%
    filter(str_detect(variable,"s\\[")) %>%
    pull(mean)
  cors <- fit_summary %>%
    filter(str_detect(variable,"cor\\[")) %>%
    mutate(
      pair=case_when(
        str_detect(variable,"cor\\[1,2\\]")~"tc",
        str_detect(variable,"cor\\[1,3\\]")~"td",
        str_detect(variable,"cor\\[2,3\\]")~"cd"
      )
    ) %>%
    filter(!is.na(pair)) %>%
    select(pair,mean)
  if(eq_cors){
    cors$mean <- mean(cors$mean) # if correlations are supposed to be equal constrain them this way
  }
  if(vary_decoy_means){
    mu_avgs_fully_collapsed_w_hdis <- mu_avgs_fully_collapsed_w_hdis %>%
      mutate(m=case_when( # from conway and cohen paper
        stim=="d" & distance==2 ~ -0.025,
        stim=="d" & distance==5 ~ -0.095,
        stim=="d" & distance==9 ~ -0.20,
        stim=="d" & distance==14 ~ -0.26,
        T~m
      ))
  }
  a <- ( s %*% t(s) )
  cv <- matrix(c(1, cors$mean[1], cors$mean[2],
                 cors$mean[1], 1, cors$mean[3],
                 cors$mean[2], cors$mean[3], 1), nrow=3, ncol=3,byrow=T)*a
  dists <- unique(mu_avgs_fully_collapsed_w_hdis$distance)
  sim <- tibble()
  do_order <- function(x){
    y <- str_replace_all(as.character(order(x,decreasing = T)),
                         c("1"="t",
                           "2"="c",
                           "3"="d"))
    z <- str_flatten(y)
    return(z)
  }
  for(d in dists){
    print(d)
    mu_tmp <- mu_avgs_fully_collapsed_w_hdis %>%
      filter(distance==d) 
    mu_tmp1 <- c(mu_tmp %>%
                   filter(stim=="t") %>%
                   pull(m),
                 mu_tmp %>%
                   filter(stim=="c") %>%
                   pull(m),
                 mu_tmp %>%
                   filter(stim=="d") %>%
                   pull(m))
    if(const_tc){
      mu_tmp1[1:2] <- mean(mu_tmp1[1:2])
    }
    x <- rmvnorm(N,mu_tmp1,cv)
    
    
    y <- apply(x,1,do_order)
    sim <- bind_rows(sim,
                     tibble(
                        order=y
                      ) %>%
                        group_by(order) %>%
                        summarise(n=n()) %>%
                        ungroup() %>%
                        mutate(prop=n/sum(n),
                               distance=d,
                               disp_cond=cond)
    )
    
  }
  return(sim)
}

f <- here("analysis","sim_from_bayes_circle_area","bayes_circle_area",which_model,glue("bw_preds_ordering_{which_model}_{outl}{tc}{cor}{dec}.csv",
                                                                                       tc=ifelse(const_tc,"_const_tc",""),
                                                                                       cor=ifelse(eq_cors,"_eq_cors",""),
                                                                                       dec=ifelse(vary_decoy_means,"_vary_decoy_means","")))
if(file_exists(f)){
  model_sims <- read_csv(f)
}else{
  model_sims <- load_and_run_model(N,cond="triangle",const_tc=const_tc,eq_cors=eq_cors,outl = outl,vary_decoy_means=vary_decoy_means) # only simming triangle cond
  write_csv(model_sims, file=f)
}
ggplot(model_sims,aes(reorder(order,-prop),prop))+
  geom_col(position="dodge",fill="lightblue")+
  facet_grid(distance~.)+
  labs(x="order")+
  ggthemes::theme_few()

ggsave(filename=here("analysis","sim_from_bayes_circle_area","bayes_circle_area",which_model,glue("bw_preds_ordering_{which_model}_{outl}{tc}{cor}{dec}.jpeg",
                                                                                                  tc=ifelse(const_tc,"_const_tc",""),
                                                                                                  cor=ifelse(eq_cors,"_eq_cors",""),
                                                                                                  dec=ifelse(vary_decoy_means,"_vary_decoy_means",""))),width=6,height=8)

model_sims %>% 
  arrange(distance,prop) %>%
  print(n=24)
