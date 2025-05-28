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

for(which_model in c("maxdiff_1")){
  model_dir <- here("analysis","bayes",which_model)
  load(path(model_dir,glue("{which_model}_fit.RData")))
  p_rank <- extract(fit,"p_rank")$p_rank
  save(p_rank,file=path(model_dir,"p_rank.RData"))
}
