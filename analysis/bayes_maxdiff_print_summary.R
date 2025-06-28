# pre-setup ====================================================================================
rm(list=ls())
library(here)
library(tidyverse)
library(glue)
library(fs)


# general setup ====================================================================================
# model settings
which_model <- "maxdiff_1"
debug_model <- F

# paths
model_dir <- here("analysis","bayes",which_model)
fit_file <- path(model_dir,glue("{which_model}_fit_summary.RData"))
load(fit_file)

print(fit_summary$summary[1:8,c("mean","sd","2.5%","97.5%")])
