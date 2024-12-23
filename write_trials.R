rm(list=ls())
library(tidyverse)
library(readxl)
library(here)
library(glue)
library(fs)

trials <- here("specs","choice_trials.xlsx") %>%
                   read_excel() %>%
                   mutate(across(c(w1,h1,w2,h2,w3,h3),~replace_na(.x,0)),
                          distance=replace_na(as.numeric(distance),0))
trial_params_dir <- c(here("specs","choice_trial_params"),
                      here("experiment_code","choice_trial_params"))
ndir <- length(trial_params_dir)
dir_create(trial_params_dir)
for(i in 1:ncol(trials)){
  print(i)
  for(j in 1:ndir){
    tmp <- unname(unlist(as.vector(trials[,i])))
    tmpdir <- trial_params_dir[j]
    tmp_file <- glue("{tmpdir}/{colnames(trials)[i]}.txt")
    fmt <- ifelse(is.numeric(tmp),"%.5f","%s")
    write_lines(sprintf(fmt,tmp),
                file = tmp_file,
                append = F,
                sep = "\n")
  }
  
}

