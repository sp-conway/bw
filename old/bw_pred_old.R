rm(list=ls())
library(tidyverse)
library(here)
library(fs)
library(mvtnorm)
library(glue)

# options
which_disp <- "triangle"
which_diag <- "both" # "lower", "upper", or "both" diagonals. "both" just includes both but collapses across
do_ranking <- F
n_sim <- 100000 # number of simulations per trial type

# results directory
save_dir <- dir_create(glue("bw_sim_from_old_{which_disp}"))

# file / directory names
model_folder <- here("old_modeling_results",glue("{which_disp}_att_rep"))

# load model params
load(path(model_folder,"sigma_mvn.RData"))
load(path(model_folder,"cor_mat_mu_rep.RData"))
load(path(model_folder,"mu_pa_z_rep_collapsed.RData"))

# correlation matrix for the equal cor case
# arbitrarily picking .6 for all
cor_mat_mu_eq <- matrix(c(1,.6,.6,.6,1,.6,.6,.6,1),nrow=3,byrow = T)

# indices for trials
ind <- switch (which_diag,
  "upper" = str_detect(trial_type,"upper"),
  "lower" = str_detect(trial_type,"lower"),
  "both" = 1:length(trial_type)
)

# grab data we're using
mu_pa_z_rep_collapsed <- mu_pa_z_rep_collapsed[ind,]
trial_type <- trial_type[ind]

# compute variance-covariance matrix
get_cov_mat <- function(Sigma, Omega){
  sigma_mu <- apply(Sigma,2,mean)
  a <- sigma_mu %*% t(sigma_mu)
  cov_mat <- Omega*a
  return(cov_mat)
}
cov_mat_fitted_cor <- get_cov_mat(sigma_mvn, cor_mat_mu)
cov_mat_eq_cor <- get_cov_mat(sigma_mvn, cor_mat_mu_eq)

# sim trials for a single condition
sim_trial <- function(n, mu, cv, trial_type, do_ranking=F){
  print(trial_type)
  s <- rmvnorm(n, mu, cv)
  mn <- apply(s, 1, which.min)
  mx <- apply(s, 1, which.max)
  x <- tibble(
    min=recode(mn,
               `1`="t",
               `2`="c",
               `3`="d"),
    max=recode(mx,
               `1`="t",
               `2`="c",
               `3`="d"),
    trial_type=trial_type
  )
  if(do_ranking){
    get_second <- function(s) which(s==(sort(s)[2]) )
    scnd <- apply(s, 1, get_second)
    x <- x %>% 
      mutate(second=recode(scnd,
                               `1`="t",
                               `2`="c",
                               `3`="d"))
  }
  return(x)
}

# wrapper function to sim trials for all conds.
sim_trial_wrapper <- function(n_sim, trial_type, mu_rep, cov_mat, cor_type, do_ranking=F){
  print(cor_type)
  # browser()
  results <- map_dfr(1:length(trial_type),
                     ~sim_trial(n_sim, 
                                mu_rep[.x,],
                                cov_mat, trial_type[.x], do_ranking=do_ranking))
  if(do_ranking){
    results1 <- results %>%
      unite("rank",c(max, second, min),sep="-") %>%
      separate(trial_type,into=c("set","distance","diag"),sep="_") %>%
      group_by(set, distance, rank) %>%
      summarise(n=n()) %>% 
      group_by(set, distance) %>%
      mutate(prop=n/sum(n)) %>%
      ungroup() %>%
      select(-n) %>%
      pivot_wider(names_from = rank, values_from = prop)  %>%
      mutate(cor_type=cor_type)
  }else{
    # browser()
    results1 <- results %>%
      pivot_longer(-trial_type,values_to = "option",names_to = "category")  %>%
      separate(trial_type,into=c("set","distance","diag"),sep="_") %>%
      group_by(set, distance, category, option) %>%
      summarise(n=n()) %>% 
      group_by(set, distance, category) %>%
      mutate(prop=n/sum(n)) %>%
      ungroup() %>%
      select(-n) %>%
      pivot_wider(names_from = category, values_from = prop)  %>%
      mutate(cor_type=cor_type)
  }
  return(results1)
}

results <- map2_dfr(list(cov_mat_fitted_cor, cov_mat_eq_cor), 
                    list("fitted","equal"), 
                    ~sim_trial_wrapper(n_sim, trial_type, mu_pa_z_rep_collapsed, .x, .y, do_ranking=do_ranking))
if(do_ranking){
  save(results,file=path(save_dir,glue("{which_diag}_rank_sim.RData")))
  do_rank_plot <- function(d, cond, diag){
    r1 <- d %>%
      filter(cor_type==cond) %>%
      pivot_longer(-c(set,distance,cor_type),values_to = "prop",names_to = "rank") %>%
      group_by(set, distance) %>%
      arrange(desc(prop)) %>%
      mutate(ord=1:n()) %>%
      ungroup()
    ggplot(r1,aes(x=ord, y=prop))+
      geom_col(fill="lightblue",col="black")+
      geom_text(aes(label=rank, y=prop+.05))+
      facet_grid(as.numeric(distance)~set,scales="free",drop=T)+
      # scale_x_discrete(labels = r1[, setNames(as.character(id), ord)]) + 
      ggthemes::theme_few()
    ggsave(filename=path(save_dir,glue("rank_props_{cond}{diag}.jpeg")),width=6,height=5)
  }
  map(c("equal","fitted"),~do_rank_plot(results, .x, which_diag))
  
  results %>%
    pivot_longer(-c(cor_type,set,distance),names_to = "rank",values_to = "prop") %>%
    pivot_wider(names_from = cor_type, values_from = prop) %>%
    mutate(diff=fitted-equal)%>%
    group_by(set, distance) %>%
    arrange(desc(diff)) %>%
    mutate(ord=1:n()) %>%
    ungroup() %>%
    ggplot(aes(x=ord, y=diff))+
    geom_col(fill="lightblue",col="black")+
    geom_text(aes(label=rank, y=if_else(diff<0, .01, diff+.01)))+
    facet_grid(as.numeric(distance)~set,scales="free",drop=T)+
    # scale_x_discrete(labels = r1[, setNames(as.character(id), ord)]) + 
    ggthemes::theme_few()
  ggsave(filename=path(save_dir,glue("rank_diff_{which_diag}.jpeg")),width=6,height=5)
  
}else{
  save(results,file=path(save_dir,glue("{which_diag}_bw_sim.RData")))
  do_str_plot <- function(d, cond, diag){
    d %>%
      filter(cor_type==cond) %>%
      ggplot(aes(min,max))+
      geom_text(aes(label=option), col="darkblue",alpha=.75, size=5)+
      # geom_point(alpha=.7,size=3)+
      see::scale_color_okabeito()+
      coord_fixed(xlim=c(.15,.55),ylim=c(.15,.55))+
      labs(x="p(min)",y="p(max)")+
      ggthemes::theme_few()+
      facet_grid(as.factor(as.numeric(distance))~set)+
      theme(text = element_text(size=14))
    ggsave(filename=path(save_dir,glue("bw_props_statetrace_{cond}{diag}.jpeg")),width=6,height=5)
  }
  map(c("equal","fitted"),~do_str_plot(results, .x, which_diag))
  
  results %>%
    pivot_longer(cor_type) %>%
    pivot_wider(names_from = value,
                values_from = c(max,min)) %>%
    select(-name) %>%
    mutate(max_diff=max_fitted-max_equal,
           min_diff=min_fitted-min_equal) %>%
    ggplot(aes(min_diff, max_diff))+
    geom_vline(xintercept=0,linetype="dashed",alpha=.4)+
    geom_hline(yintercept=0,linetype="dashed",alpha=.4)+
    geom_text(aes(label=option),size=5, col="purple",alpha=.65)+
    # geom_point(aes(col=option),alpha=.7,size=3)+
    see::scale_color_okabeito()+
    coord_fixed(xlim=c(-.05, .05),ylim=c(-.05, .05))+
    scale_x_continuous(breaks = c(-.05,0,.05))+
    scale_y_continuous(breaks = c(-.05,0,.05))+
    facet_grid(as.factor(as.numeric(distance))~set)+
    labs(y="p(max|fitted)-p(max|equal)",x="p(min|fitted)-p(min|equal)")+
    ggthemes::theme_few()+
    theme(text = element_text(size=14))
  ggsave(filename=path(save_dir,glue("bw_diffprops_{which_diag}.jpeg")),width=6,height=5)
  
  results %>%
    filter(option!="d") %>%
    pivot_wider(names_from = option,
                values_from = c(max,min)) %>%
    mutate(max_rst=max_t/(max_t+max_c),
           min_rst=min_t/(min_t+min_c)) %>%
    select(-c(max_t,min_t,max_c,min_c)) %>%
    group_by(distance,cor_type) %>%
    summarise(max_rst=mean(max_rst),
              min_rst=mean(min_rst)) %>%
    ggplot(aes(min_rst,max_rst))+
    geom_point(col="purple",alpha=.5)+
    # geom_text(aes(label=distance,col=set),alpha=.75, size=5)+
    # geom_point(alpha=.7,size=3)+
    see::scale_color_okabeito()+
    coord_fixed(xlim=c(.45,.52),ylim=c(.45,.52))+
    labs(x="p(min)",y="p(max)")+
    ggthemes::theme_few()+
    facet_grid(.~cor_type)+
    theme(text = element_text(size=14))
}
