rm(list=ls())
library(here)
library(tidyverse)
library(glue)
library(fs)
library(patchwork)
library(furrr)
library(multinomineq)


# read in data, only include critical trials, recode bw cond to be more understandable
d <- here("data","cleaned","bw_all.csv") %>%
  read_csv() %>%
  mutate(bw_cond=recode(bw_cond,"bw"="best-worst","wb"="worst-best")) %>%
  filter(str_detect(effect,"attraction"))
order_levels <- c("tcd","tdc","cdt","ctd","dtc","dct")
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

strong_corr_A <- read_model_A(here("analysis/order_constraints/strong_correlations.txt"),
                              skip=2,nlines=7)
strong_corr_B <- read_model_B(here("analysis/order_constraints/strong_correlations.txt"),
                              skip=10,nlines=7)

weak_corr_A <- read_model_A(here("analysis/order_constraints/weak_correlations.txt"),
                            skip=2,nlines=8)
weak_corr_B <- read_model_B(here("analysis/order_constraints/weak_correlations.txt"),
                            skip=11,nlines=8)

repulsion_A <- read_model_A(here("analysis/order_constraints/repulsion.txt"),
                            skip=2,nlines=6)
repulsion_B <- read_model_B(here("analysis/order_constraints/repulsion.txt"),
                            skip=9,nlines=6)

attraction_A <- read_model_A(here("analysis/order_constraints/attraction.txt"),
                             skip=2,nlines=7)
attraction_B <- read_model_B(here("analysis/order_constraints/attraction.txt"),
                            skip=10,nlines=7)

null_A <- read_model_A(here("analysis/order_constraints/null.txt"),
                       skip=2,nlines=15)
null_B <- read_model_B(here("analysis/order_constraints/null.txt"),
                       skip=18,nlines=15)

models <- list(
  "strong_correlations"=list(strong_corr_A,
                             strong_corr_B),
  "weak_correlations"=list(weak_corr_A,
                           weak_corr_B),
  "attraction"=list(attraction_A,
                    attraction_B),
  "repulsion"=list(repulsion_A,
                   repulsion_B),
  "null"=list(null_A,
              null_B)
)

run_models <- function(data, models, sub_n, distance, M=1e3){
  options <- c(6) # always 6
  n_models <- length(models)
  model_names <- names(models)
  i <- 1
  # results_counts <- vector("list",n_models)
  results <- tibble(
    model=character(),
    comparison=character(),
    bf=numeric(),
    se=numeric(),
    `ci.5%`=numeric(),
    `ci.95%`=numeric()
  )
  for(m in models){
    # browser()
    prior_count <- count_multinom(k=0,options=c(6), A=m[[1]],b=m[[2]],M=M,progress=T)
    x <- 1
    do_sample <- T
    if(do_sample){
      x <- x+1
      posterior <- count_multinom(k=data,options=c(6),A=m[[1]], b=m[[2]],
                                  M=M,progress = T)
      tmp_bf <- count_to_bf(posterior,prior_count)
      if(all(is.finite(tmp_bf[,1]))|x>5){
        do_sample <- F
      }
    }
    
    tmp_bf_1 <- as_tibble(tmp_bf) %>%
      mutate(comparison=rownames(tmp_bf),
             model=model_names[i],
             prop_inside=attr(posterior,'proportion'),
             prop_inside_se=attr(posterior,'se'))
    results <- bind_rows(results,tmp_bf_1)
    # results_counts[[i]] <- posterior
    i <- i+1
  }
  results$distance <- distance
  results$sub_n <- sub_n
  return(results)
  # return(list(results,results_counts))
}
run_models_wrapper <- function(data_all, models, distance_cond, M){
  sub_ns <- unique(data_all$sub_n)
  n_subs <- length(sub_ns)
  data_all_filtered <- data_all %>%
    filter(distance==distance_cond)
  results <- vector("list",n_subs)
  i <- 1
  for(s in sub_ns){
    print(distance_cond)
    print(s)
    results[[i]] <- run_models(filter(data_all_filtered,sub_n==s) %>%
                                 select(-c(sub_n,distance)) %>%
                                 unlist(as.vector(.)),
                               models=models,
                               sub_n=s,
                               distance=distance_cond,
                               M=M)
    i <- i+1
  }
  results <- list_rbind(results) %>%
    relocate(c(sub_n,distance,model),.before=everything())
  return(results)
}
M <- 4e6 # IMPORTANT - NUMBER OF SAMPLES

f <- here("analysis/order_constraints/results",glue("bf_{M}_samples.RData"))
if(file_exists(f)){
  load(f)
}else{
  plan(multisession, workers=4)
  model_results <- future_map(c(2,5,9,14),
                              ~run_models_wrapper(d_order_counts_wide,
                                                  models, 
                                                  .x, 
                                                  M=M),
                              .options = furrr_options(seed = T,stdout=T))
  model_results <- list_rbind(model_results)
  save(model_results, file=f)
}

model_summaries <- model_results %>% 
  filter(is.finite(bf) & comparison=="bf_0u") %>%
  group_by(distance, model) %>%
  summarise(bf_multipled=sum(log(bf)))
