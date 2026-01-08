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
vary_decoy_means <- T
plot_data <- T
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
        stim=="d" & distance==2 ~ -0.025,# from conway and cohen paper
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
              upper=m_prop+qt(.975,n-1)*se) %>%
    ungroup() %>%
    select(-c(n,se)) %>%
    pivot_wider(names_from = type,
                values_from = c(m_prop,lower,upper)) %>%
    mutate(source="data",
           distance=str_glue("{distance}% TDD"),
           distance=factor(distance,levels=c("2% TDD",
                                             "5% TDD",
                                             "9% TDD",
                                             "14% TDD")))
}

model_output_file <- here("analysis","sim_from_bayes_circle_area","bayes_circle_area",which_model,glue("bw_preds_{which_model}_{outl}{dec}.csv",
                                                                                                                        dec=ifelse(vary_decoy_means,"_vary_decoy_means","")))
if(file_exists(model_output_file)){
  model_sims_wide <- read_csv(model_output_file)
}else{
  model_sims <- load_and_run_model(N,"triangle",outl = outl,vary_decoy_means=vary_decoy_means) # only simming triangle cond
  model_sims_wide <- model_sims %>%
    pivot_wider(names_from = type, values_from = p, names_prefix = "m_prop_")  %>%
    mutate(distance=str_glue("{distance}% TDD"),
           distance=factor(distance,levels=c("2% TDD",
                                             "5% TDD",
                                             "9% TDD",
                                             "14% TDD")),
           option=factor(option,levels=c("t","c","d")))
  write_csv(model_sims_wide,file=model_output_file)
}

if(plot_data){
  data <- analyze_data()
  model_sims_wide <- model_sims_wide %>%
    mutate(source="mtcm") %>%
    bind_rows(data)
}

model_sims_wide %>%
  mutate(option=factor(option,levels=c("t","c","d")),
         distance=factor(distance,levels=c("2% TDD",
                                           "5% TDD",
                                           "9% TDD",
                                           "14% TDD"))) %>%
  ggplot(aes(m_prop_worst,m_prop_best))+
  geom_point(aes(col=option,shape=source),size=3.5)+
  scale_shape_manual(name="",values=c(1,4)) +
  coord_fixed(xlim=c(0,.8),ylim=c(0,.8))+
  # scale_x_continuous(breaks=c(0,.5,1))+
  # scale_y_continuous(breaks=c(0,.5,1))+
  # annotate("point", x = 0.1, y = 0.75, shape = 17, size = 3) +
  # annotate("text", x = 0.12, y = 0.75, label = "Data", hjust = 0, size = 5) +
  # annotate("point", x = 0.1, y = 0.7, shape = 4, size = 3) +
  # annotate("text", x = 0.12, y = 0.7, label = "Model", hjust = 0, size = 5)+
  ggsci::scale_color_startrek(name="")+
  labs(x="p(worst)",y="p(best)")+
  facet_wrap(vars(distance),nrow=2,ncol=2)+
  ggthemes::theme_few()+
  theme(text=element_text(size=22),
        legend.spacing.x = unit(0, 'cm'),
        legend.direction = 'horizontal',
        legend.title = element_blank(),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.width = unit(.01,'cm'),
        legend.key.height = unit(.01,'cm'),
        legend.position="inside",
        legend.position.inside=c(.25,.91))
ggsave(filename=here("analysis","sim_from_bayes_circle_area","bayes_circle_area",which_model,glue("bw_preds_{which_model}_{outl}{dec}.jpeg",
                                                                                                  dec=ifelse(vary_decoy_means,"_vary_decoy_means",""))),width=8,height=8)

