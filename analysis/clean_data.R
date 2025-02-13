rm(list=ls())
library(here)
library(tidyverse)
library(fs)
library(glue)

data_dir_raw <- here("data","raw")
data_dir_cleaned <- here("data","cleaned")
dir_create(data_dir_cleaned)
data_folders_raw <- dir_ls(data_dir_raw)
data_files_raw_all <- map(data_folders_raw, function(x) dir_ls(x)) %>%
  list_c()

clean <- function(f,n_trials=598){ # finished ppts should have 598 rows
  print(f)
  comp <- as.numeric(str_extract(f, "(?<=raw/c).{1,20}(?=/bw)"))
  d <- data.table::fread(f) %>%
    mutate(sub_n=str_extract(sub_n,"[:digit:]{1,3}"))
  d$sub_n <- as.numeric(glue("{unique(d$sub_n)}{comp}"))
  if(nrow(d)<n_trials){
    return(NULL)
  }else{
    dd <- d %>%
      group_by(block_n,trial_n) %>%
      mutate(choice_n=1:n()) %>%
      ungroup() %>%
      pivot_wider(names_from = choice_n,
                  values_from = c(choice,rt),
                  names_glue = "{.value}_{choice_n}") %>%
      mutate(choice_best=if_else(bw_cond=="bw",choice_1,choice_2),
             choice_worst=if_else(bw_cond=="bw",choice_2,choice_1),
             rt_best=if_else(bw_cond=="bw",rt_1,rt_2),
             rt_worst=if_else(bw_cond=="bw",rt_2,rt_1),
             a1=h1*w1,
             a2=h2*w2,
             a3=h3*w3,
             set=recode(set,"a-b-da"="h","a-b-db"="w"), # RECODING TO MAKE SURE ITS UNDERSTANDABLE
             across(c(diag,distance),function(x) na_if(x,999)),
             across(c(rt_best,rt_worst),function(x) x*1000)) %>%
      select(-c(choice_1,choice_2,rt_1,rt_2,computer_n)) %>%
      relocate(c(block_n,trial_n),.after=bw_cond) %>%
      relocate(c(a1,a2,a3),.after = w3) %>%
      filter(bw_cond=="bw" & (rt_best >=100 & rt_best<=10000) | # FILTER OUT SLOW AND FAST RTS BUT ONLY FOR FIRST CHOICE
             bw_cond=="wb" & (rt_worst >=100 & rt_worst<=10000) )
    return(dd)
  }
}

d_all <- map_dfr(data_files_raw_all,clean)

d_corr <- d_all %>%
  mutate(
    a_max=apply(cbind(d_all$a1,d_all$a2,d_all$a3),1,max),
    a_min=apply(cbind(d_all$a1,d_all$a2,d_all$a3),1,min),
    choice_best_correct=case_when(
      choice_best==1 & a1==a_max~1,
      choice_best==2 & a2==a_max~1,
      choice_best==3 & a3==a_max~1,
      T~0
    ),
    choice_worst_correct=case_when(
      choice_worst==1 & a1==a_min~1,
      choice_worst==2 & a2==a_min~1,
      choice_worst==3 & a3==a_min~1,
      T~0
    ),
    both_correct=as.integer(choice_worst_correct+choice_best_correct==2)
  )

d_corr_props <- d_corr %>%
  group_by(sub_n,effect) %>%
  summarise(prop_best=mean(choice_best_correct),
            prop_worst=mean(choice_worst_correct),
            prop_both=mean(both_correct)) %>%
  ungroup()

d_corr_props_pass_catch <- d_corr_props %>%
  filter(effect=="catch") %>%
  filter(prop_both>.8)
subs_keep <- unique(d_corr_props_pass_catch$sub_n)
d_all_cleaned <- filter(d_all, sub_n %in% subs_keep)

cat("\n=====\nremoved",length(unique(d_all$sub_n))-length(subs_keep),"participants for failing catch trials\n=====\n",sep=" ")

d_all_cleaned %>%
  distinct(bw_cond,sub_n) %>%
  group_by(bw_cond) %>%
  summarise(N=n())
write_csv(d_all_cleaned,file=path(data_dir_cleaned,"bw_all.csv"))
