rm(list=ls())
library(here)
library(tidyverse)
library(glue)
library(fs)
library(patchwork)
library(furrr)
library(scales)
library(multinomineq)

do_models <- T
# whether or not to do parallel
do_parallel_bf <- F
do_parallel_post <- F

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
  mutate(middle_att=setdiff(c("t","c","d"),c(best_att,worst_att))) %>%
  ungroup() %>%
  mutate(order=factor(str_c(best_att,middle_att,worst_att),
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
nonmonotonic_repulsion_strong_A <- read_model_A(here("analysis/order_constraints/non-monotonic_repulsion_strong.txt"),
                                         skip=2,nlines=13)
nonmonotonic_repulsion_strong_B <- read_model_B(here("analysis/order_constraints/non-monotonic_repulsion_strong.txt"),
                                         skip=16,nlines=13)
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

run_model_bf <- function(data, model, sub_n, distance, M_init=1e3){
  
  M <- M_init
  # Number of options always 6
  options <- c(6) # always 6
  
  # name of model we are running
  model_name <- model[[3]]
  
  # limit of number of samples taken from posterior. just to avoid computer overload. hopefully won't need it
  samp_limit <- 500000000
  
  # sample from the prior
  prior_count <- count_multinom(k=0,options=c(6), A=model[[1]],b=model[[2]],M=M,progress = T)
  
  # initial sample from the posterior
  posterior <- count_multinom(k=data,options=c(6),A=model[[1]], b=model[[2]],
                              M=M,progress = T)
  
  # find bayes factor and check if not finite or if =0
  tmp_bf <- count_to_bf(posterior,prior_count)
  if(tmp_bf['bf_0u',1]==0 | !is.finite(tmp_bf['bf_0u',1])){ # re-sample as needed
    M <- 5000000000
    do_sample <- T
    while(do_sample){
      # print("re-sampling")
      posterior <- count_multinom(k=data,options=c(6),A=model[[1]], b=model[[2]],
                                  M=M)
      tmp_bf <- count_to_bf(posterior,prior_count)
      if(tmp_bf['bf_0u',1]!=0 & is.finite(tmp_bf['bf_0u',1])){
        do_sample <- F
      }else if(M>samp_limit){
        print("hit sample limit")
        do_sample <- F
      }else{
        M <- M+100000000
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
  
  results$distance <- distance
  results$sub_n <- sub_n
  return(results)
}

run_model_bf_wrapper <- function(data_all, model_current, distance_cond, results_file, M_init=1e3){
  sub_ns <- unique(data_all$sub_n)
  n_subs <- length(sub_ns)
  data_all_filtered <- data_all %>%
    filter(distance==distance_cond)
  i <- 1
  results_df <- read.csv(results_file)
  # browser()
  for(s in sub_ns){
    tmp <- results_df %>%
      filter(sub_n==s &
             distance==distance_cond &
             model==model_current[[3]])
    if(nrow(tmp)==0){
      cat(distance_cond,"% Distance","\n")
      cat(i,"/",n_subs," Subjects\n")
      results_tmp <- run_model_bf(filter(data_all_filtered,sub_n==s) %>%
                                     select(-c(sub_n,distance)) %>%
                                     unlist(as.vector(.)),
                                   model=model_current,
                                   sub_n=s,
                                   distance=distance_cond,
                                   M_init=M_init)
      write_csv(results_tmp,
                file=results_file,
                append=T)
    }
    i <- i+1
  }
}


run_model_post <- function(data, model, sub_n, distance, M=1e3){
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
    sub_n=sub_n,
    distance=distance,
    model=model_name,
    obs=ppp[1],
    pred=ppp[2],
    ppp=ppp[3]
  )
  return(results)
}

run_model_post_wrapper <- function(data_all, model_current, distance_cond, results_file, M=1e3){
  sub_ns <- unique(data_all$sub_n)
  n_subs <- length(sub_ns)
  data_all_filtered <- data_all %>%
    filter(distance==distance_cond)
  results_df <- read.csv(results_file)
  i <- 1
  for(s in sub_ns){
    # browser()
    tmp <- results_df %>%
      filter(sub_n==s &
               distance==distance_cond &
               model==model_current[[3]])
    if(nrow(tmp)==0){
      cat(distance_cond,"% Distance","\n")
      cat(i,"/",n_subs," Subjects\n")
      results_tmp <- run_model_post(filter(data_all_filtered,sub_n==s) %>%
                                      select(-c(sub_n,distance)) %>%
                                      unlist(as.vector(.)),
                                    model=model_current,
                                    sub_n=s,
                                    distance=distance_cond,
                                    M=M)
      write_csv(results_tmp,
                file=results_file,
                append=T)
    }
    i <- i+1
  }
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
# IMPORTANT - NUMBER OF SAMPLES 
M_init_post <- 50000 # need a lot less for posterior distributions, if p=0 it's okay
M_init_bf <- 500000

results_file_bf <- path(results_dir,"bf.csv")
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
    sub_n=numeric(),
  )
  write_csv(tmp_bf, file=results_file_bf)
}

results_file_post <- path(results_dir,"post.csv")
if(!file_exists(results_file_post)){
  tmp_post <- tibble(
    sub_n=numeric(),
    distance=numeric(),
    model=character(),
    obs=numeric(),
    pred=numeric(),
    ppp=numeric()
  )
  write_csv(tmp_post, file=results_file_post)
}

if(do_models){
  for(model_tmp in models){
    
    model_name_tmp <- model_tmp[[3]]
    
    if(do_parallel_bf){
      print(paste("Sampling BF for model:",model_name_tmp))
      plan(multisession, workers=4)
      future_map(c(2,5,9,14),~run_model_bf_wrapper(d_order_counts_wide,
                                                   model_tmp,
                                                   .x,
                                                   results_file=results_file_bf,
                                                   M=M_init_bf),
                 .options = furrr_options(seed = T,stdout=T))
    }else{
      print(paste("Sampling BF for model:",model_name_tmp))
      map(c(2,5,9,14),~run_model_bf_wrapper(d_order_counts_wide,
                                            model_tmp,
                                            .x,
                                            results_file=results_file_bf,
                                            M=M_init_bf))
    }
    
    
    if(do_parallel_post){
      print(paste("Sampling Posterior for model:",model_name_tmp))
      
      plan(multisession, workers=4)
      future_map(c(2,5,9,14), ~run_model_post_wrapper(d_order_counts_wide,
                                                      model_tmp, 
                                                      .x, 
                                                      results_file=results_file_post,
                                                      M=M_init_post),
                 .options = furrr_options(seed=T,stdout=T))
    }else{
      print(paste("Sampling Posterior for model:",model_name_tmp))
      
      map(c(2,5,9,14), ~run_model_post_wrapper(d_order_counts_wide,
                                               model_tmp, 
                                               .x, 
                                               results_file=results_file_post,
                                               M=M_init_post))
    }
  }
}



# # load results ========================================================================================================================
# load_results_bf <- function(f){
#   if(!str_detect(f,"jpeg|pdf")){
#     load(f)
#     return(model_results_bf)
#   }
# }
# load_results_post <- function(f){
#   if(!str_detect(f,"jpeg|pdf")){
#     load(f)
#     model_results_post_1 <- model_results_post %>%
#       list_rbind()
#     return(model_results_post_1)
#   }
# }

# model_levels <- c("none",
#                  "repulsion",
#                  "attraction",
#                  "non-monotonic_repulsion",
#                  "non-monotonic_repulsion_strong",
#                  "non-monotonic_attraction")#,
#                  #"similarity")
# model_labels <- c(
#   "none",
#   "repulsion",
#   "attraction",
#   "non-monotonic repulsion",
#   "non-monotonic repulsion-strong",
#   "non-monotonic attraction")#,
#   #"similarity"
# #)
model_results_all_bf <- here("analysis/order_constraints/results/bf.csv") %>%
  read_csv()%>%
  filter(comparison=="bf_0u")

model_results_all_post <-  here("analysis/order_constraints/results/post.csv") %>%
  read_csv()
#
# # # first examine bayes factors ================================================================
#
model_results_all_bf_data <- model_results_all_bf %>%
  left_join(d_order_counts_wide)

model_results_all_bf_data %>%
  group_by(distance,model) %>%
  summarise(n_zero=sum(bf==0),
            n_not_finite=sum(!is.finite(bf))) %>%
  ungroup() %>%
  arrange(model) %>%
  print(n=nrow(.))
model_results_all_bf_max <- model_results_all_bf_data %>%
  group_by(sub_n,distance) %>%
  mutate(max_bf=max(bf,na.rm=T)) %>%
  ungroup() %>%
  filter(bf==max_bf) %>%
  ungroup() %>%
  left_join(model_results_all_post) %>%
  mutate(model=case_when(
    bf<1~"none",
    T~model
  ))
model_results_all_bf_max %>%
  group_by(distance)%>%
  summarise(n=n()) %>%
  ungroup()

model_counts <- model_results_all_bf_max %>%
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
  print(n=30)

model_counts %>%
  mutate(model=factor(model,
                      levels=model_levels,
                      labels=model_labels)) %>%
  arrange(distance, desc(n)) %>%
  ggplot(aes(model,n))+
  geom_col(position="dodge",fill="lightblue")+
  labs(x="model",y="N preferred")+
  coord_flip()+
  facet_wrap(vars(distance),nrow=2)+
  ggthemes::theme_few()+
  theme(text=element_text(size=14))
ggsave(filename = path(results_dir,glue("bf_counts_{M_init}_samples.pdf")),
       width=5,height=5)
# # #
model_summaries  <- model_results_all_bf_max %>%
  mutate(distance=str_glue("{distance}% TDD"),
         distance=factor(distance,
                         levels=c("2% TDD",
                                  "5% TDD",
                                  "9% TDD",
                                  "14% TDD"))) %>%
  mutate(bf=case_when(
    bf==0~1e-323,
    T~bf
  )) %>%
  filter(model!="none") %>%
  group_by(distance, model) %>%
  summarise(bf_joint_log=sum(log10(bf)),
            bf_joint=10^bf_joint_log) %>%
  ungroup() %>%
  arrange(distance,model,desc(bf_joint_log))


model_summaries %>%
  mutate(model=factor(model,
                      levels=model_levels,
                      labels=model_labels)) %>%
  ggplot(aes(model,bf_joint_log))+
  geom_col(position="dodge",fill="lightblue",width=.75)+
  geom_hline(yintercept=0,alpha=.85,linetype="dashed")+
  labs(x="model",y="log10 joint bayes factor")+
  coord_flip()+
  facet_wrap(vars(distance),nrow=2,scales="free_x")+
  ggthemes::theme_few()+
  theme(text=element_text(size=14))
ggsave(filename = path(results_dir,glue("bf_joint_{M_init}_samples.pdf")),
       width=5,height=5)

# posteriors ============================================================
model_results_all_bf_max_w_post <- model_results_all_bf_max %>%
  left_join(model_results_all_post)
model_results_all_bf_max_w_post %>%
  mutate(reject=case_when(
    ppp<.05 ~ "reject",
    ppp>=.05~"fail to reject"
  ),
  reject=factor(reject,levels=c("reject","fail to reject"))) %>%
  group_by(model,distance,reject) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  pivot_wider(names_from = reject,
              values_from = n,
              values_fill = 0) %>%
  print(n=30)

model_results_all_bf_max_w_post %>%
  ggplot(aes(ppp))+
  geom_histogram(fill="lightblue")+
  facet_grid(model~distance)+
  ggthemes::theme_few()

model_results_all_bf_max %>%
  filter(model=="none") %>%
  select(sub_n,distance,tcd,tdc,ctd,cdt,dtc,dct) %>%
  print(n=30)
