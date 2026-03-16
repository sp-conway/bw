rm(list=ls())
library(tidyverse)
library(fs)
library(mvtnorm)
library(here)
library(matrixcalc)
library(latex2exp)
library(furrr)
mu_T_all <- seq(.5,1,.1)

sample <- function(N, mu, cv){
  X <- rmvnorm(N, mu, cv)
  # find which max and min
  mx <- apply(X, 1, which.max)
  mn <- apply(X, 1, which.min)
  pmax <- c(sum(mx==1)/N,
            sum(mx==2)/N,
            sum(mx==3)/N)
  pmin <- c(sum(mn==1)/N,
            sum(mn==2)/N,
            sum(mn==3)/N)
  tibble(
    option=c("t","c","d"),
    b=pmax,
    w=pmin
  )
}

sample_wrapper <- function(mu_T){
  # fixed params
  s <- c(1,1,1)
  a <- ( s %*% t(s) )
  # varying params
  mu_C_all <- seq(.5,1,.1)
  mu_D_diff_all <- seq(.5,0,-.1)
  rho_TD_all <- seq(.5,1,.1)
  rho_TC_CD_all <- seq(.5,1,.1)
  i <- 1
  sim <- vector("list")
  N <- 500000
  K <- length(mu_T_all)*length(mu_C_all)*length(mu_D_diff_all)*length(rho_TD_all)*length(rho_TC_CD_all)
  for(mu_C in mu_C_all){
    for(mu_D_diff in mu_D_diff_all){
      for(rho_TD in rho_TD_all){
        for(rho_TC in rho_TC_CD_all){
          cat(i,"/",K,"\n==========================\n")
          mu_D <- mu_T-mu_D_diff
          mu <- c(mu_T, mu_C, mu_D)
          rho_CD <- rho_TC
          cv <- matrix(c(1, rho_TC, rho_TD,
                         rho_TC, 1, rho_CD,
                         rho_TD, rho_CD, 1), nrow=3, ncol=3,byrow=T)*a
          if(is.positive.semi.definite(cv)){
            # st <- Sys.time()
            sim[[i]] <- sample(N, mu, cv) %>%
              mutate(mu_T=mu_T,
                     mu_C=mu_C,
                     mu_D=mu_D,
                     rho_TD=rho_TD,
                     rho_TC_CD=rho_CD)
            i <- i+1
            end <- Sys.time()
            # print(end-st)
          }
        }
      }
    }
  }
  return(sim)
}
f <- here("analysis/sim_effects_params/sim_effects.RData")
if(!file_exists(f)){
  plan(multisession)
  sim <- future_map(mu_T_all,~sample_wrapper(.x),.options = furrr_options(seed = NULL))
  save(sim,file=f)
}else{
  load(f)
}


sim_clean <- list_rbind(sim)
sim_effects <- sim_clean %>%
  mutate(across(c(b,w),\(x) round(x, digits=2))) %>%
  pivot_wider(names_from = option,
              values_from = c(b,w)) %>%
  mutate(effect=case_when(
    b_t>b_c & w_c>w_t & b_c>b_d & b_t>b_d & w_d>w_t & w_d>w_c ~ "attraction",
    b_c>b_t & w_c>w_t & b_c>b_d & w_d>w_c & b_t>b_d & w_d>w_t ~ "non-monotonic repulsion",
    b_c>b_t & w_t>w_c & b_c>b_d & w_d>w_c & b_t>b_d & w_d>w_t ~ "repulsion",
    b_t>b_c & w_t>w_c & b_c>b_d & w_d>w_c & b_t>b_d & w_d>w_t ~ "non-monotonic attraction",
    b_c>b_t & w_c>w_t & b_c>b_d & w_c>w_d & abs(b_t-b_d)<=.01 & abs(w_t-w_d)<=.01 ~ "similarity",
    T~"other"
  ),
  mu_td=
  case_when(
    mu_T>mu_D ~ "mu_T>mu_D",
    mu_T==mu_D ~ "mu_T=mu_D"
  ),
  mu_tc=case_when(
    mu_T>mu_C ~ "mu_T>mu_C",
    mu_C>mu_T ~ "mu_C>mu_T",
    mu_T==mu_C~ "mu_T=mu_C"
  ),
  rho_pattern=case_when(
    rho_TD>rho_TC_CD~"rho_TD>rho_TC=rho_CD",
    rho_TD<rho_TC_CD~"rho_TD<rho_TC=rho_CD",
    rho_TD==rho_TC_CD~"rho_TD=rho_TC=rho_CD"
  ),
  mu_pattern=str_glue("{mu_td}, {mu_tc}"))
#
#
# table(sim_effects$effect,useNA = "ifany")
sim_params_given_effects_counts <- sim_effects %>%
  group_by(effect,rho_pattern,mu_pattern) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  filter(n>10) %>%
  pivot_wider(names_from = effect, values_from = n, values_fill = 0) %>%
  pivot_longer(-c(rho_pattern,mu_pattern),
               names_to = "effect",
               values_to = "n") %>%
  group_by(effect) %>%
  mutate(prop=n/sum(n)) %>%
  ungroup() %>%
  mutate(effect=factor(effect,
                       levels=c("attraction","non-monotonic attraction",
                                "repulsion","non-monotonic repulsion",
                                "similarity","other"),
                       labels=c("attraction","non-monotonic\nattraction",
                                "repulsion","non-monotonic\nrepulsion",
                                "similarity","other")))
sim_effects_given_params_counts <- sim_effects %>%
  group_by(effect,rho_pattern,mu_pattern) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  filter(n>10) %>%
  pivot_wider(names_from = effect, values_from = n, values_fill = 0) %>%
  pivot_longer(-c(rho_pattern,mu_pattern),
               names_to = "effect",
               values_to = "n") %>%
  group_by(mu_pattern,rho_pattern) %>%
  mutate(prop=n/sum(n)) %>%
  ungroup() %>%
  mutate(effect=factor(effect,
                       levels=c("attraction","non-monotonic attraction",
                                "repulsion","non-monotonic repulsion",
                                "similarity","other"),
                       labels=c("attraction","non-monotonic\nattraction",
                                "repulsion","non-monotonic\nrepulsion",
                                "similarity","other")))
params_given_effects_plot <- ggplot(sim_params_given_effects_counts,
                                    aes(mu_pattern,rho_pattern,fill=prop))+
  geom_tile()+
  scale_fill_gradient(low="black",high="white")+
  coord_flip()+
  labs(x=TeX("$\\mu$"),y=TeX("$\\rho$"),
       title="Params | Effects")+
  scale_x_discrete(labels=rev(c(
    TeX("$\\mu_{T}>\\mu_{D},\\mu_{T}>\\mu_{C}$"),
    TeX("$\\mu_{T}>\\mu_{D},\\mu_{T}=\\mu_{C}$"),
    TeX("$\\mu_{T}>\\mu_{D},\\mu_{C}>\\mu_{T}$"),
    TeX("$\\mu_{T}=\\mu_{D},\\mu_{T}>\\mu_{C}$"),
    TeX("$\\mu_{T}=\\mu_{D},\\mu_{T}=\\mu_{C}$"),
    TeX("$\\mu_{T}=\\mu_{D},\\mu_{C}>\\mu_{T}$")
  )))+
  scale_y_discrete(labels=c(
    TeX("$\\rho_{TD}<\\rho_{TC}=\\rho_{CD}$"),
    TeX("$\\rho_{TD}=\\rho_{TC}=\\rho_{CD}$"),
    TeX("$\\rho_{TD}>\\rho_{TC}=\\rho_{CD}$")
  ))+
  facet_grid(effect~.)+
  ggthemes::theme_few()+
  theme(text=element_text(size=12))
params_given_effects_plot
ggsave(params_given_effects_plot,filename=here("analysis/sim_effects_params/sim_params_given_effects.jpeg"),
       width=8,height=8)


effects_given_params_plot <- ggplot(sim_effects_given_params_counts,
                                    aes(mu_pattern,rho_pattern,fill=prop))+
  geom_tile()+
  scale_fill_gradient(low="black",high="white")+
  coord_flip()+
  labs(x=TeX("$\\mu$"),y=TeX("$\\rho$"),title="Effects | Params")+
  scale_x_discrete(labels=rev(c(
    TeX("$\\mu_{T}>\\mu_{D},\\mu_{T}>\\mu_{C}$"),
    TeX("$\\mu_{T}>\\mu_{D},\\mu_{T}=\\mu_{C}$"),
    TeX("$\\mu_{T}>\\mu_{D},\\mu_{C}>\\mu_{T}$"),
    TeX("$\\mu_{T}=\\mu_{D},\\mu_{T}>\\mu_{C}$"),
    TeX("$\\mu_{T}=\\mu_{D},\\mu_{T}=\\mu_{C}$"),
    TeX("$\\mu_{T}=\\mu_{D},\\mu_{C}>\\mu_{T}$")
  )))+
  scale_y_discrete(labels=c(
    TeX("$\\rho_{TD}<\\rho_{TC}=\\rho_{CD}$"),
    TeX("$\\rho_{TD}=\\rho_{TC}=\\rho_{CD}$"),
    TeX("$\\rho_{TD}>\\rho_{TC}=\\rho_{CD}$")
  ))+
  facet_grid(effect~.)+
  ggthemes::theme_few()+
  theme(text=element_text(size=12))
ggsave(effects_given_params_plot,filename=here("analysis/sim_effects_params/sim_effects_given_params.jpeg"),
       width=8,height=8)

