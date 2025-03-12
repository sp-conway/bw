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
# model settings
which_model <- "maxdiff_separateU_1"
debug_model <- F

# paths
model_dir <- here("analysis","bayes",which_model)
stan_file <- path(model_dir,glue("{which_model}.stan"))

load(path(model_dir,glue("{which_model}_counts_rank_rep.RData")))
load(path(model_dir,glue("{which_model}_data_for_model.RData")))
