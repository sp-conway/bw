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

N <- 1e6
outl <- "no_outliers"
which_model <- "sigma_constant_comp_effect"
vary_decoy_means <- F
load_and_run_model <- function(N,cond,outl="no_outliers",vary_decoy_means=F){
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
  if(vary_decoy_means){
    mu_avgs_fully_collapsed_w_hdis <- mu_avgs_fully_collapsed_w_hdis %>%
      mutate(m=case_when(
        stim=="d" & distance==2 ~ -0.025,
        stim=="d" & distance==5 ~ -0.095,
        stim=="d" & distance==9 ~ -0.20,
        stim=="d" & distance==14 ~ -0.26,
        T~m
      ))
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

model_sims <- load_and_run_model(N,"triangle",outl = outl,vary_decoy_means=vary_decoy_means) # only simming triangle cond
model_sims_wide <- model_sims %>%
  pivot_wider(names_from = type, values_from = p, names_prefix = "m_prop_")  %>%
  mutate(distance=str_glue("{distance}% TDD"),
         distance=factor(distance,levels=c("2% TDD",
                                           "5% TDD",
                                           "9% TDD",
                                           "14% TDD")),
         option=factor(option,levels=c("t","c","d")))

model_sims_wide %>%
  ggplot(aes(m_prop_worst,m_prop_best,col=option,label=str_sub(option,1,1)))+
  geom_text(size=6,alpha=.8)+
  # scale_x_continuous(breaks=c(25,.5,1))+
  # scale_y_continuous(breaks=c(0,.5,1))+
  coord_fixed(xlim=c(0,.5),ylim=c(0,.5))+
  ggsci::scale_color_startrek(name="stimulus")+
  labs(x="p(worst)",y="p(best)")+
  facet_grid(distance~.)+
  ggthemes::theme_few()+
  theme(text=element_text(size=18),
        legend.position = "none")
ggsave(filename=here("analysis","sim_from_bayes_circle_area","bayes_circle_area",which_model,glue("bw_preds_{which_model}_{outl}{dec}.jpeg",
                                                                                                  dec=ifelse(vary_decoy_means,"_vary_decoy_means",""))),width=6,height=8)
write_csv(model_sims_wide, file=here("analysis","sim_from_bayes_circle_area","bayes_circle_area",which_model,glue("bw_preds_{which_model}_{outl}{dec}.csv",
                                                                                                                  dec=ifelse(vary_decoy_means,"_vary_decoy_means",""))))

model_sims_wide %>%
  filter(option!="d") %>%
  select(-m_prop_best) %>%
  pivot_wider(names_from = option,
              values_from = m_prop_worst) %>%
  mutate(diff=c-t)

