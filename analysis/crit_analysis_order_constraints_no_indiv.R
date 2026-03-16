# setup ====================================================================================
rm(list=ls())
library(here)
library(tidyverse)
library(glue)
library(fs)
library(patchwork)
library(furrr)
library(scales)
library(multinomineq)


# whether or not to do parallel
do_parallel_bf <- T
do_parallel_post <- F

# where to store results
results_dir <- here("analysis/order_constraints/results/")

# process data ====================================================================================
# read in data, only include critical trials, recode bw cond to be more understandable
d <- here("data","cleaned","bw_all.csv") %>%
  read_csv() %>%
  mutate(bw_cond=recode(bw_cond,"bw"="best-worst","wb"="worst-best")) %>%
  filter(str_detect(effect,"attraction"))
order_levels <- c("tcd","tdc","cdt","ctd","dtc","dct")
# IMPORTANT
# TCD=p1
# TDC=p2
# CDT=p3
# CTD=p4
# DTC=p5
# DCT=p6

# critical choice props 
d_order <- d %>%
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
  ))%>% 
  rowwise() %>%
  mutate(middle_att=setdiff(c("t","c","d"),c(best_att,worst_att))) %>%
  ungroup() %>%
  mutate(order=factor(str_c(best_att,middle_att,worst_att),
                      levels=order_levels)) %>%
  ungroup()

d_order_props <- d_order %>%
  group_by(distance,order) %>%
  summarise(n=n()) %>%
  group_by(distance) %>%
  mutate(prop=n/sum(n)) %>%
  ungroup() 

d_order_counts_wide <- d_order_props %>%
  arrange(distance) %>%
  select(-prop) %>%
  pivot_wider(names_from = order, values_from = n, values_fill = 0)

# import models ====================================================================================
read_model_A <- function(f, skip, nlines,ncol=5){
  m <- scan(f,skip=skip,nlines=nlines) %>%
    matrix(nrow = nlines, ncol = ncol,byrow=T)
  colnames(m) <- c("p1","p2","p3","p4","p5")
  return(m)
}
read_model_B <- function(f, skip, nlines){
  m <- scan(f, skip=skip, nlines=nlines)
  return(m)
}

repulsion_A <- read_model_A(here("analysis/order_constraints/repulsion.txt"),
                            skip=2,nlines=10)
repulsion_B <- read_model_B(here("analysis/order_constraints/repulsion.txt"),
                            skip=13,nlines=10)
attraction_A <- read_model_A(here("analysis/order_constraints/attraction.txt"),
                             skip=2,nlines=14)
attraction_B <- read_model_B(here("analysis/order_constraints/attraction.txt"),
                             skip=17,nlines=14)
nonmonotonic_repulsion_A <- read_model_A(here("analysis/order_constraints/non-monotonic_repulsion.txt"),
                            skip=2,nlines=14)
nonmonotonic_repulsion_B <- read_model_B(here("analysis/order_constraints/non-monotonic_repulsion.txt"),
                            skip=17,nlines=14)
# nonmonotonic_repulsion_strong_A <- read_model_A(here("analysis/order_constraints/non-monotonic_repulsion_strong.txt"),
#                                          skip=2,nlines=13)
# nonmonotonic_repulsion_strong_B <- read_model_B(here("analysis/order_constraints/non-monotonic_repulsion_strong.txt"),
#                                          skip=16,nlines=13)
nonmonotonic_attraction_A <- read_model_A(here("analysis/order_constraints/non-monotonic_attraction.txt"),
                                         skip=2,nlines=8)
nonmonotonic_attraction_B <- read_model_B(here("analysis/order_constraints/non-monotonic_attraction.txt"),
                                         skip=11,nlines=8)
similarity_A <- read_model_A(here("analysis/order_constraints/similarity.txt"),
                             skip=2,nlines=15)
similarity_B <- read_model_B(here("analysis/order_constraints/similarity.txt"),
                             skip=18,nlines=15)
models <- list(
  list(attraction_A,
       attraction_B,
       "attraction"),
  list(repulsion_A,
       repulsion_B,
       "repulsion"),
  list(nonmonotonic_repulsion_A,
       nonmonotonic_repulsion_B,
       "non-monotonic_repulsion"),
  # list(nonmonotonic_repulsion_strong_A,
  #      nonmonotonic_repulsion_strong_B,
  #      "non-monotonic_repulsion_strong"),
  list(nonmonotonic_attraction_A,
       nonmonotonic_attraction_B,
       "non-monotonic_attraction"),
  list(similarity_A,
       similarity_B,
       "similarity")
)

# functions for running models ====================================================================================
run_model_bf <- function(data, model, dist, M_init=1e3){
  print(model[[3]])
  print(dist)
  
  data <- filter(data,distance==dist) %>%
    select(-distance) 
  data <- unlist(as.vector(data))

  M <- M_init
  # Number of options always 6
  options <- c(6) # always 6
  
  # name of model we are running
  model_name <- model[[3]]
  
  # limit of number of samples taken from posterior. just to avoid computer overload. hopefully won't need it
  samp_limit <- 5e20
  
  # sample from the prior
  prior_count <- count_multinom(k=0,options=c(6), A=model[[1]],b=model[[2]],M=M,progress = T)
  
  # initial sample from the posterior
  posterior <- count_multinom(k=data,options=c(6),A=model[[1]], b=model[[2]],
                              M=M,progress = T)
  
  # find bayes factor and check if not finite or if =0
  tmp_bf <- count_to_bf(posterior,prior_count)
  if(tmp_bf['bf_0u',1]==0 | !is.finite(tmp_bf['bf_0u',1])){ # re-sample as needed
    M <- 5e+19
    do_sample <- T
    while(do_sample){
      print("re-sampling")
      posterior <- count_multinom(k=data,options=c(6),A=model[[1]], b=model[[2]],
                                  M=M)
      tmp_bf <- count_to_bf(posterior,prior_count)
      if(tmp_bf['bf_0u',1]!=0 & is.finite(tmp_bf['bf_0u',1])){
        do_sample <- F
      }else if(M>samp_limit){
        print("hit sample limit")
        do_sample <- F
      }else{
        M <- M+5e+19
      }
    }
  }
  results <- tibble(
    model=character(),
    comparison=character(),
    bf=numeric(),
    se=numeric(),
    `ci.5%`=numeric(),
    `ci.95%`=numeric(),
    M=numeric()
  )
  
  tmp_bf_1 <- as_tibble(tmp_bf) %>%
    mutate(comparison=rownames(tmp_bf),
           model=model_name,
           prop_inside=attr(posterior,'proportion'),
           prop_inside_se=attr(posterior,'se'),
           M=M)
  results <- bind_rows(results,tmp_bf_1)
  
  results$distance <- dist
  return(results)
}



run_model_post <- function(data, model, dist, M=1e3){
  data <- filter(data,distance==dist) %>%
    select(-distance) 
  data <- unlist(as.vector(data))
  options <- c(6) # always 6
  model_name <- model[[3]]
  # IMPORTANT
  # TCD=p1
  # TDC=p2
  # CDT=p3
  # CTD=p4
  # DTC=p5
  # DCT=p6
  # browser()
  posterior <- sampling_multinom(k=data,options=c(6), A=model[[1]],b=model[[2]],M=M,progress=T)
  ppp <- ppp_multinom(posterior, k=data, c(6))
  results <- tibble(
    distance=dist,
    model=model_name,
    obs=ppp[1],
    pred=ppp[2],
    ppp=ppp[3]
  )
  # browser()
  return(results)
}

# run analyses ============================================================================================================================================================
# IMPORTANT - NUMBER OF SAMPLES 
M_init_post <- 50000000 # need a lot less for posterior distributions, if p=0 it's okay
M_init_bf <- 1e6


results_all_bf <- results_all_post <- vector("list")
i<-1
results_file_bf <- path(results_dir,"bf_aggregated.csv")
if(!file_exists(results_file_bf)){
  tmp_bf <- tibble(
    model=character(),
    comparison=character(),
    bf=numeric(),
    se=numeric(),
    `ci.5%`=numeric(),
    `ci.95%`=numeric(),
    M=numeric(),
    prop_inside=numeric(),
    prop_inside_se=numeric(),
    distance=numeric(),
  )
  write_csv(tmp_bf, file=results_file_bf)
}

results_file_post <- path(results_dir,"post_aggregated.csv")
if(!file_exists(results_file_post)){
  tmp_post <- tibble(
    distance=numeric(),
    model=character(),
    obs=numeric(),
    pred=numeric(),
    ppp=numeric()
  )
  write_csv(tmp_post, file=results_file_post)
}

for(model_tmp in models){
  results_bf <- read_csv(results_file_bf)
  results_post <- read_csv(results_file_post)
    model_name_tmp <- model_tmp[[3]]
    results_bf_tmp <- results_bf %>%
      filter(model==model_name_tmp)
    if(nrow(results_bf_tmp)==0){
      if(do_parallel_bf){
        print(paste("Sampling BF for model:",model_name_tmp))
        plan(multisession, workers=4)
        results_tmp <- future_map(c(2,5,9,14),~run_model_bf(d_order_counts_wide,
                                                            model_tmp,
                                                            .x,
                                                            M=M_init_bf),
                                  .options = furrr_options(seed = T,stdout=T))
        results_all_bf[[i]] <- list_rbind(results_tmp)
      }else{
        print(paste("Sampling BF for model:",model_name_tmp))
        results_tmp <- map(c(2,5,9,14),~run_model_bf(d_order_counts_wide,
                                                     model_tmp,
                                                     .x,
                                                     M=M_init_bf),
                           .options = furrr_options(seed = T,stdout=T))
        results_all_bf[[i]] <- list_rbind(results_tmp)
      }
      write_csv(results_all_bf[[i]], file=results_file_bf, append=T)
    }
    
    
    results_post_tmp <- results_post %>%
      filter(model==model_name_tmp)
    if(nrow(results_post_tmp)==0){
      if(do_parallel_post){
        print(paste("Sampling Posterior for model:",model_name_tmp))
        
        plan(multisession, workers=4)
        results_tmp <-future_map(c(2,5,9,14), ~run_model_post(d_order_counts_wide,
                                                        model_tmp, 
                                                        .x, 
                                                        M=M_init_post),
                   .options = furrr_options(seed=T,stdout=T))
        results_all_post[[i]] <- list_rbind(results_tmp)
      }else{
        print(paste("Sampling Posterior for model:",model_name_tmp))
        
        results_tmp <- map(c(2,5,9,14), ~run_model_post(d_order_counts_wide,
                                         model_tmp, 
                                         .x, 
                                         M=M_init_post))
        results_all_post[[i]] <- list_rbind(results_tmp)
      }
      write_csv(results_all_post[[i]], file=results_file_post, append=T)
    }
    i <- i+1
}


# # load results ========================================================================================================================
# model_results_bf <- here("analysis/order_constraints/results/bf.csv") %>%
#   read_csv() %>%
#   group_by(model, sub_n, distance) %>%
#   filter(M==max(M)) %>%
#   ungroup() %>%
#   filter(comparison=="bf_0u") %>%
#   select(sub_n,distance,model,bf)
# 
# model_results_post <-  here("analysis/order_constraints/results/post.csv") %>%
#   read_csv() %>%
#   select(sub_n,distance,model,ppp)
# 
# model_results_bf_post <- model_results_bf %>%
#   left_join(model_results_post) %>%
#   arrange(sub_n,distance,desc(bf))
# 
# # # # first examine bayes factors ================================================================
# model_results_bf %>%
#   group_by(distance,model) %>%
#   summarise(n_zero=sum(bf==0),
#             n_not_finite=sum(!is.finite(bf))) %>%
#   ungroup() %>%
#   arrange(model) %>%
#   print(n=nrow(.))
# 
# sort_bfs <- function(results){
#   results <- arrange(results, desc(bf))
#   results$best <- logical(nrow(results))
#   for(i in 1:nrow(results)){
#    if(results$ppp[i]>.05){
#      results$best[i] <- T
#      break
#    }
#   }
#   results_1 <- results %>%
#     mutate(model=case_when(
#       bf<=1~"none",
#       T~model)) %>%
#     filter(best) %>%
#     select(-best)
#   return(results_1)
# }
# 
# model_results_bf_post_split <- model_results_bf_post %>%
#   group_by(sub_n,distance) %>%
#   group_split() 
# 
# model_results_bf_post_split_w_max <- map(model_results_bf_post_split,sort_bfs)
# model_results_bf_post_w_max <- model_results_bf_post_split_w_max %>%
#   bind_rows()
# model_levels <- c("none",
#                   "repulsion",
#                   "non-monotonic_repulsion",
#                   "attraction",
#                   "non-monotonic_attraction",
#                   "similarity")
# model_labels <- c("none",
#                   "repulsion",
#                   "non-monotonic repulsion",
#                   "attraction",
#                   "non-monotonic attraction",
#                   "similarity")
# 
# 
# model_counts <- model_results_bf_post_w_max %>%
#   group_by(model,distance) %>%
#   summarise(n=n()) %>%
#   ungroup() %>%
#   group_by(distance) %>%
#   mutate(perc=100*(n/sum(n))) %>%
#   ungroup() %>%
#   arrange(distance, model, desc(n)) %>%
#   mutate(distance=str_glue("{distance}% TDD"),
#          distance=factor(distance,
#                          levels=c("2% TDD",
#                                   "5% TDD",
#                                   "9% TDD",
#                                   "14% TDD")),
#          model=factor(model,levels=model_levels,
#                       labels=model_labels))
# 
# 
# model_counts %>%
#   arrange(distance, desc(n)) %>%
#   ggplot(aes(model,n))+
#   geom_col(position="dodge",fill="lightblue")+
#   labs(x="model",y="N preferred")+
#   coord_flip()+
#   facet_wrap(vars(distance),nrow=2)+
#   ggthemes::theme_few()+
#   theme(text=element_text(size=14))
# ggsave(filename = path(results_dir,glue("bf_counts.pdf")),
#        width=5,height=5)
# # # #
# model_summaries_by_distance  <- model_results_bf %>%
#   mutate(distance=str_glue("{distance}% TDD"),
#          distance=factor(distance,
#                          levels=c("2% TDD",
#                                   "5% TDD",
#                                   "9% TDD",
#                                   "14% TDD"))) %>%
#   filter(model!="none") %>%
#   mutate(l=log10(bf)) %>%
#   group_by(distance, model) %>%
#   summarise(bf_joint_log=sum(log10(bf)),
#             bf_joint=10^(bf_joint_log)) %>%
#   ungroup() %>%
#   filter(is.finite(bf_joint_log)) %>%
#   arrange(distance,model,desc(bf_joint_log))
# 
# model_summaries_by_distance %>%
#   ggplot(aes(model,bf_joint_log))+
#   geom_col(position="dodge",fill="lightblue",width=.75)+
#   geom_hline(yintercept=0,alpha=.85,linetype="dashed")+
#   coord_flip()+
#   facet_wrap(vars(distance),nrow=2,scales="free_y")+
#   ggthemes::theme_few()+
#   theme(text=element_text(size=14))
# ggsave(filename = path(results_dir,glue("bf_joint_by_distance.pdf")),
#        width=8,height=5)
# 
# model_summaries <- model_results_bf %>%
#   mutate(distance=str_glue("{distance}% TDD"),
#          distance=factor(distance,
#                          levels=c("2% TDD",
#                                   "5% TDD",
#                                   "9% TDD",
#                                   "14% TDD"))) %>%
#   filter(model!="none") %>%
#   mutate(l=log10(bf)) %>%
#   group_by(model) %>%
#   summarise(bf_joint_log=sum(log10(bf)),
#             bf_joint=10^(bf_joint_log)) %>%
#   ungroup() %>%
#   filter(is.finite(bf_joint_log)) %>%
#   arrange(model,desc(bf_joint_log))
# 
# model_summaries %>%
#   ggplot(aes(model,bf_joint_log))+
#   geom_col(position="dodge",fill="lightblue",width=.75)+
#   geom_hline(yintercept=0,alpha=.85,linetype="dashed")+
#   coord_flip()+
#   ggthemes::theme_few()+
#   theme(text=element_text(size=14))
# ggsave(filename = path(results_dir,glue("bf_joint_all.pdf")),
#        width=8,height=5)
# model_summaries
