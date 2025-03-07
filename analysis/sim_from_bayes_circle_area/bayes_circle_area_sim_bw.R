rm(list=ls())
library(tidyverse)
library(here)
library(glue)
library(mvtnorm)
library(fs)

N <- 20000
outl <- "no_outliers"

for(which_model in c("sigma_constant","sigma_constant_target_effect","sigma_constant_comp_effect")){
  load_and_run_model <- function(N,cond,outl="no_outliers"){
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
    a <- ( s %*% t(s) )
    cv <- matrix(c(1, cors$mean[1], cors$mean[2],
                   cors$mean[1], 1, cors$mean[3],
                   cors$mean[2], cors$mean[3], 1), nrow=3, ncol=3,byrow=T)*a
    dists <- unique(mu_avgs_fully_collapsed_w_hdis$distance)
    sim <- tibble()
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
      x <- rmvnorm(N,mu_tmp1,cv)
      mx <- apply(x, 1, which.max)
      mn <- apply(x, 1, which.min)
      pmax <- c(sum(mx==1)/N,
                sum(mx==2)/N,
                sum(mx==3)/N)
      pmin <- c(sum(mn==1)/N,
                sum(mn==2)/N,
                sum(mn==3)/N)
      sim <- bind_rows(sim,
                       tibble(option=c("t","c","d"),
                              disp_cond=cond,
                              distance=d,
                              p=pmax,
                              type="best"),
                       tibble(option=c("t","c","d"),
                              disp_cond=cond,
                              distance=d,
                              p=pmin,
                              type="worst"))
    }
    return(sim)
  }
  
  model_sims <- map(c("horizontal","triangle"),~load_and_run_model(N,.x,outl = outl)) %>%
    list_rbind()
  model_sims_wide <- model_sims %>%
    mutate(disp_cond=factor(disp_cond,levels=c("triangle","horizontal"))) %>%
    pivot_wider(names_from = type, values_from = p, names_prefix = "m_prop_") 
  
  model_sims_wide %>%
    ggplot(aes(m_prop_worst,m_prop_best,col=option,label=str_sub(option,1,1)))+
    geom_text(size=3,alpha=.8)+
    coord_fixed(xlim=c(0,1),ylim=c(0,1))+
    ggsci::scale_color_startrek(name="stimulus")+
    labs(x="p(worst)",y="p(best)")+
    facet_grid(distance~disp_cond)+
    ggthemes::theme_few()+
    theme(text=element_text(size=28),
          legend.position = "none")
  ggsave(filename=here("analysis","sim_from_bayes_circle_area","bayes_circle_area",which_model,glue("bw_preds_{which_model}_{outl}.jpeg")),width=12,height=6)
  write_csv(model_sims_wide, file=here("analysis","sim_from_bayes_circle_area","bayes_circle_area",which_model,glue("bw_preds_{which_model}_{outl}.csv")))
  
}
