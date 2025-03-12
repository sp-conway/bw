# pre-setup ====================================================================================
rm(list=ls()); gc()
library(here)
library(tidyverse)
library(glue)
library(HDInterval)
library(fs)
library(rstan)
library(bayesplot)


# general setup ====================================================================================

for(which_model in c("maxdiff_1","maxdiff_2","maxdiff_3","maxdiff_separateU_1","maxdiff_separateU_2")){
  
  # paths
  model_dir <- here("analysis","bayes",which_model)
  load(path(model_dir,glue("{which_model}_fit.RData")))
  p_rank <- extract(fit,"p_rank")$p_rank
  save(p_rank,file=path(model_dir,"p_rank.RData"))
  
}
