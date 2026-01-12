rm(list=ls())
library(here)
library(tidyverse)
library(glue)
library(fs)
library(patchwork)
library(furrr)
library(multinomineq)

results_dir <- here("analysis/order_constraints/results/")

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

# critical choice props ===========================================================================
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
  mutate(middle_att=setdiff(c("t","c","d"),c(best_att,worst_att)),
         order=factor(str_c(best_att,middle_att,worst_att),
                      levels=order_levels)) %>%
  ungroup()

d_order_props <- d_order %>%
  group_by(sub_n,distance,order) %>%
  summarise(n=n()) %>%
  group_by(sub_n,distance) %>%
  mutate(prop=n/sum(n)) %>%
  ungroup() 
  
d_order_counts_wide <- d_order_props %>%
  arrange(sub_n,distance) %>%
  select(-prop) %>%
  pivot_wider(names_from = order, values_from = n, values_fill = 0)

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

corr_1_A <- read_model_A(here("analysis/order_constraints/correlations_1.txt"),
                              skip=2,nlines=7)
corr_1_B <- read_model_B(here("analysis/order_constraints/correlations_1.txt"),
                              skip=10,nlines=7)

corr_2_A <- read_model_A(here("analysis/order_constraints/correlations_2.txt"),
                            skip=2,nlines=8)
corr_2_B <- read_model_B(here("analysis/order_constraints/correlations_2.txt"),
                            skip=11,nlines=8)

repulsion_A <- read_model_A(here("analysis/order_constraints/repulsion.txt"),
                            skip=2,nlines=6)
repulsion_B <- read_model_B(here("analysis/order_constraints/repulsion.txt"),
                            skip=9,nlines=6)

attraction_A <- read_model_A(here("analysis/order_constraints/attraction.txt"),
                             skip=2,nlines=7)
attraction_B <- read_model_B(here("analysis/order_constraints/attraction.txt"),
                            skip=10,nlines=7)

# null_A <- read_model_A(here("analysis/order_constraints/null.txt"),
#                        skip=2,nlines=15)
# null_B <- read_model_B(here("analysis/order_constraints/null.txt"),
#                        skip=18,nlines=15)

models <- list(
  list(corr_1_A,
       corr_1_B,
       "correlations_1"),
  list(corr_2_A,
       corr_2_B,
       "correlations_2"),
  list(attraction_A,
       attraction_B,
       "attraction"),
  list(repulsion_A,
       repulsion_B,
       "repulsion")#,
  # list(null_A,
  #      null_B,
  #      "null")
)

run_model_bf <- function(data, model, sub_n, distance, M=1e3){
  options <- c(6) # always 6
  model_name <- model[[3]]
  results <- tibble(
    model=character(),
    comparison=character(),
    bf=numeric(),
    se=numeric(),
    `ci.5%`=numeric(),
    `ci.95%`=numeric()
  )
  prior_count <- count_multinom(k=0,options=c(6), A=model[[1]],b=model[[2]],M=M,progress=T)
  posterior <- count_multinom(k=data,options=c(6),A=model[[1]], b=model[[2]],
                              M=M,progress = T)
  tmp_bf <- count_to_bf(posterior,prior_count)
  tmp_bf_1 <- as_tibble(tmp_bf) %>%
    mutate(comparison=rownames(tmp_bf),
           model=model_name,
           prop_inside=attr(posterior,'proportion'),
           prop_inside_se=attr(posterior,'se'))
  results <- bind_rows(results,tmp_bf_1)
  
  results$distance <- distance
  results$sub_n <- sub_n
  return(results)
}

run_model_bf_wrapper <- function(data_all, model, distance_cond, M=1e3){
  sub_ns <- unique(data_all$sub_n)
  n_subs <- length(sub_ns)
  data_all_filtered <- data_all %>%
    filter(distance==distance_cond)
  results <- vector("list",n_subs)
  i <- 1
  for(s in sub_ns){
    print(distance_cond)
    print(s)
    results[[i]] <- run_model_bf(filter(data_all_filtered,sub_n==s) %>%
                                 select(-c(sub_n,distance)) %>%
                                 unlist(as.vector(.)),
                               model=model,
                               sub_n=s,
                               distance=distance_cond,
                               M=M)
    i <- i+1
  }
  results <- list_rbind(results) %>%
    relocate(c(sub_n,distance,model),.before=everything())
  return(results)
}


run_model_post <- function(data, model, sub_n, distance, M=1e3){
  options <- c(6) # always 6
  model_name <- model[[3]]
  # browser()
  posterior <- sampling_multinom(k=0,options=c(6), A=model[[1]],b=model[[2]],M=M,progress=T)
  ppp <- ppp_multinom(posterior, data, c(6))
  ppp_results <- tibble(
    sub_n=sub_n,
    model=model_name,
    distance=distance,
    obs=ppp[1],
    pred=ppp[2],
    ppp=ppp[3],
  )
  return(list(posterior,ppp_results))
}

run_model_post_wrapper <- function(data_all, model, distance_cond, M=1e3){
  sub_ns <- unique(data_all$sub_n)
  n_subs <- length(sub_ns)
  data_all_filtered <- data_all %>%
    filter(distance==distance_cond)
  results <- vector("list",n_subs)
  i <- 1
  for(s in sub_ns){
    print(distance_cond)
    print(s)
    results[[i]] <- run_model_post(filter(data_all_filtered,sub_n==s) %>%
                                select(-c(sub_n,distance)) %>%
                                unlist(as.vector(.)),
                              model=model,
                              sub_n=s,
                              distance=distance_cond,
                              M=M)
    i <- i+1
  }
  # results <- list_rbind(results) %>%
  #   relocate(c(sub_n,distance,model),.before=everything())
  return(results)
}

# testing ================================================================================================================================================
# 
# run_model_post_wrapper(filter(d_order_counts_wide,
#                  sub_n %in% c(22,3922) & distance==2),
#           model=models[[1]],
#           distance=2,M=)
# 
# run_model_post_wrapper(filter(d_order_counts_wide,
#                               sub_n %in% c(22,3922) & distance==2) %>%
#                          select(-c(sub_n,distance)) %>%
#                          unlist(as.vector(.)),
#                        model=models[[1]],
#                        sub_n=22,distance=2)
# run analyses ============================================================================================================================================================

for(model_tmp in models){
  M <- 50000 # IMPORTANT - NUMBER OF SAMPLES
  model_name_tmp <- model_tmp[[3]]
  f_tmp_post <- path(results_dir,
                   glue("{model_name_tmp}_post_{M}_samples.RData"))
  f_tmp_bf <- path(results_dir,
                glue("{model_name_tmp}_bf_{M}_samples.RData"))
  if(!file_exists(f_tmp_bf)){
    plan(multisession, workers=4)
    model_results_bf <-  future_map(c(2,5,9,14),
                                ~run_model_bf_wrapper(d_order_counts_wide,
                                                    model_tmp, 
                                                    .x, 
                                                    M=M),
                                .options = furrr_options(seed = T,stdout=T))
    model_results_bf <- list_rbind(model_results_bf)
    save(model_results_bf, file=f_tmp_bf)
  }
  
  if(!file_exists(f_tmp_post)){
    plan(multisession, workers=4)
    model_results_post <-  future_map(c(2,5,9,14),
                                   ~run_model_post_wrapper(d_order_counts_wide,
                                                         model_tmp, 
                                                         .x, 
                                                         M=M),
                                   .options = furrr_options(seed = T,stdout=T))
    save(model_results_post, file=f_tmp_post)
    
   }
}

# 
load_results_bf <- function(f){
  load(f)
  return(model_results_bf)
}
models_all_bf <- results_dir %>%
  dir_ls(regexp="bf")
model_results_all_bf <- map(models_all_bf, load_results_bf) %>%
  list_rbind() %>%
  mutate(bf=case_when(
    !is.finite(bf)~median(c(`ci.5%`,`ci.95%`)),
    T~bf
  ))


load_results_post <- function(f){
  load(f)
  n_distance <- 4
  n_subs <- 369
  res <- tibble(sub_n=numeric(),
                model=character(),
                distance=numeric(),
                obs=numeric(),
                pred=numeric(),
                ppp=numeric())
  for(i in 1:n_distance){
    print(i)
    for(j in 1:n_subs){
      res <- bind_rows(res,
                       model_results_post[[i]][[j]][[2]]
      )
    }
  }
  return(res)
}
models_all_post <- results_dir %>%
  dir_ls(regexp="post")
model_results_all_post <- map(models_all_post, load_results_post) %>%
  list_rbind()
save(model_results_all_post,file=path(results_dir, glue("post_{M}_samples_clean.RData")))

#   
model_results_0u <- model_results_all_bf %>%
  filter( comparison=="bf_0u" & is.finite(bf)) 
model_results_0u_with_data <- model_results_0u %>%
  filter(is.finite(bf) & comparison=="bf_0u") %>%
  left_join(d_order_counts_wide)
# 
model_counts <- model_results_0u %>%
  group_by(sub_n,distance) %>%
  mutate(max_bf=max(bf)) %>%
  filter(bf==max(bf)) %>%
  group_by(model,distance) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  group_by(distance) %>%
  mutate(perc=100*(n/sum(n))) %>%
  ungroup() %>%
  arrange(distance, model, desc(n)) %>%
  mutate(distance=str_glue("{distance}% TDD"),
         distance=factor(distance,
                         levels=c("2% TDD",
                         "5% TDD",
                         "9% TDD",
                         "14% TDD")))

model_counts %>%
  arrange(distance, model) %>%
  ggplot(aes(model,n))+
  geom_col(position="dodge",fill="lightblue")+
  labs(x="model",y="N preferred")+
  coord_flip()+
  facet_wrap(vars(distance),nrow=2)+
  ggthemes::theme_few()+
  theme(text=element_text(size=14))
ggsave(filename = path(results_dir,glue("bf_counts_{M}_samples.jpeg")),
       width=5,height=5)
# 
model_summaries <- model_results_0u %>%
  mutate(distance=str_glue("{distance}% TDD"),
         distance=factor(distance,
                         levels=c("2% TDD",
                                  "5% TDD",
                                  "9% TDD",
                                  "14% TDD"))) %>%
  mutate(bf=case_when(
    bf==0~median(c(`ci.5%`,`ci.95%`)),
    T~bf
  )) %>%
  group_by(distance, model) %>%
  summarise(bf_joint_log10=sum(log10(bf))) %>%
  ungroup() %>%
  arrange(distance,model,desc(bf_joint_log10))
model_summaries
model_summaries$bf_joint_log10
model_summaries %>%
  ggplot(aes(model,bf_joint_log10))+
  geom_col(position="dodge",fill="lightblue",width=.75)+
  geom_hline(yintercept=0,alpha=.85,linetype="dashed")+
  labs(x="model",y="log10 joint bayes factor")+
  coord_flip()+
  facet_grid(distance~.)+
  ggthemes::theme_few()+
  theme(text=element_text(size=14))
ggsave(filename = path(results_dir,glue("bf_joint_log10_{M}_samples.jpeg")),
       width=5,height=7)


model_results_all_post %>%
  left_join(model_results_all_bf) %>%
  ggplot(aes(bf,ppp))+
  geom_point(alpha=.5)+
  geom_hline(yintercept=.05,linetype="dashed",col="red")+
  facet_wrap(vars(distance))+
  ggthemes::theme_few()
