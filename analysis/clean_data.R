# cleaning data from bw experiment
# also recoding set from a-b-da to just h / w
# filter out fast vs. slow RTs
# figuring out which trials are correct or not
# also remove ppts who don't get at least 80% of catch trials correct
# to get a catch trial correct had to get both best and worst correct
rm(list=ls())
library(here)
library(tidyverse)
library(fs)
library(glue)

# SETUP ==========================================================================================
# raw data files
data_dir_raw <- here("data","raw")

# directory for cleaned data files
data_dir_cleaned <- here("data","cleaned")
dir_create(data_dir_cleaned)

# sep folder for each computer
data_folders_raw <- dir_ls(data_dir_raw)

# concat all file names
data_files_raw_all <- map(data_folders_raw, function(x) dir_ls(x)) %>%
  list_c()

# CLEANING FUNCTION ==========================================================================================
clean <- function(f,n_trials=598){ # finished ppts should have 598 rows
  print(f)
  # grab computer name
  comp <- as.numeric(str_extract(f, "(?<=raw/c).{1,20}(?=/bw)"))
  
  # read in data, get subject number, convert to numeric
  d <- data.table::fread(f) %>%
    mutate(sub_n=str_extract(sub_n,"[:digit:]{1,3}")) # occasionally an RA bumped a button e.g ppt 123]
  d$sub_n <- as.numeric(glue("{unique(d$sub_n)}{comp}"))
  if(nrow(d)<n_trials){
    return(NULL)
  }else{
    # pivot data and figure out best vs. worst choice
    dd <- d %>%
      group_by(block_n,trial_n) %>%
      mutate(choice_n=1:n()) %>%
      ungroup() %>%
      pivot_wider(names_from = choice_n,
                  values_from = c(choice,rt),
                  names_glue = "{.value}_{choice_n}") %>%
      mutate(choice_best=if_else(bw_cond=="bw",choice_1,choice_2), # bw condition means best than worst
             choice_worst=if_else(bw_cond=="bw",choice_2,choice_1),# wb condition means worst than best
             rt_best=if_else(bw_cond=="bw",rt_1,rt_2),
             rt_worst=if_else(bw_cond=="bw",rt_2,rt_1),
             a1=h1*w1,
             a2=h2*w2, # figure out areas
             a3=h3*w3,
             set=recode(set,"a-b-da"="h","a-b-db"="w"), # RECODING TO MAKE SURE ITS UNDERSTANDABLE
             across(c(diag,distance),function(x) na_if(x,999)), # I called NA 999 in experiment
             across(c(rt_best,rt_worst),function(x) x*1000)) %>%
      select(-c(choice_1,choice_2,rt_1,rt_2,computer_n)) %>%
      relocate(c(block_n,trial_n),.after=bw_cond) %>%
      relocate(c(a1,a2,a3),.after = w3) %>%
      filter(bw_cond=="bw" & (rt_best >=100 & rt_best<=10000) | # FILTER OUT SLOW AND FAST RTS BUT ONLY FOR FIRST CHOICE
             bw_cond=="wb" & (rt_worst >=100 & rt_worst<=10000) ) 
  d_all <- dd %>% 
  mutate(
    a_max=apply(cbind(dd$a1,dd$a2,dd$a3),1,max), # figure out which is max and min on each trial
    a_min=apply(cbind(dd$a1,dd$a2,dd$a3),1,min),
    choice_best_correct=case_when( # figure out correct or not. had to get best and smallest
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
    return(d_all)
  }
}

# READ IN DATA ==========================================================================================
d_all <- map_dfr(data_files_raw_all,clean)

# PROP CORRECT, FILTERING OUT CATCH-FAILING PPTS ==========================================================================================
# figure out subject / effect level prop correct
d_corr_props <- d_all %>%
  group_by(sub_n,effect) %>%
  summarise(prop_best=mean(choice_best_correct),
            prop_worst=mean(choice_worst_correct),
            prop_both=mean(both_correct)) %>%
  ungroup()

# catch trials
# figuring out who "passed"
# need at least 80% correct
d_corr_props_pass_catch <- d_corr_props %>%
  filter(effect=="catch") %>%
  filter(prop_both>.8)

# vector of unique subject numbers of only those who passed the catch trials
subs_keep <- unique(d_corr_props_pass_catch$sub_n)

# FINAL CLEANED DATA ==========================================================================================
d_all_cleaned <- filter(d_all, sub_n %in% subs_keep)

cat("\n=====\nremoved",length(unique(d_all$sub_n))-length(subs_keep),"participants for failing catch trials\n=====\n",sep=" ")

# FINAL SUBJECT COUNT==========================================================================================
d_all_cleaned %>%
  distinct(bw_cond,sub_n) %>%
  group_by(bw_cond) %>%
  summarise(N=n())

# SAVE DATA =========================================================================================
write_csv(d_all_cleaned,file=path(data_dir_cleaned,"bw_all.csv"))
