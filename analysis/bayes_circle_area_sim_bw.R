# simulating choice predictions for bw choice, conditioned on best fitting mean params from circle exp
# ONLY simulating triangle condition
# NOW ONLY INCLUDING sigma_constant_comp_effect
# removed other models from thsi project

# SETUP ====================================================================================================
rm(list=ls())
library(tidyverse)
library(here)
library(glue)
library(mvtnorm)
library(fs)

# number of draws
N <- 1e6

# which model params
which_model <- "sigma_constant_comp_effect"
model_directory <- here("analysis/bayes_circle_area",which_model,"triangle/no_outliers")

# Whether or not to use decoy means which were varied to fit choice data from Conway & Cohen
vary_decoy_means <- T

# whether to plot actual data from current experiment
plot_data <- F

# FUNCTION TO SIMULATE MTCM ======================================================================
load_and_run_model <- function(N,model_dir,vary_decoy_means=F){
  browser()

  # load model parameters
  load(path(model_dir,"mu_avgs.RData")) # MU
  load(path(model_dir,"fit_summary.RData")) # THIS HAS S AND COR
  
  # STANDARD DEVIATIONS
  s <- fit_summary %>%
    filter(str_detect(variable,"s\\[")) %>%
    pull(mean) 
  
  # CORRELATIONS
  cors <- fit_summary %>%
    filter(str_detect(variable,"cor\\[")) %>%
    mutate(
      pair=case_when(
        str_detect(variable,"cor\\[1,2\\]")~"tc", # in order of t c d in matrix
        str_detect(variable,"cor\\[1,3\\]")~"td",
        str_detect(variable,"cor\\[2,3\\]")~"cd"
      )
    ) %>%
    filter(!is.na(pair)) %>%
    select(pair,mean)
  
  # VARIANCES
  a <- ( s %*% t(s) )
  
  # VARIANCE - COVARIANCE MATRIX
  cv <- matrix(c(1, cors$mean[1], cors$mean[2],
                 cors$mean[1], 1, cors$mean[3],
                 cors$mean[2], cors$mean[3], 1), nrow=3, ncol=3,byrow=T)*a
  
  # all TDD values
  dists <- unique(mu_avgs_fully_collapsed_w_hdis$distance)
  
  # GRABBING CONWAY AND COHEN DECOY MEANS IF NEEDED
  if(vary_decoy_means){
    mu_avgs_fully_collapsed_w_hdis <- mu_avgs_fully_collapsed_w_hdis %>%
      mutate(m=case_when(
        stim=="d" & distance==2 ~ -0.025,# from conway and cohen paper appendix e
        stim=="d" & distance==5 ~ -0.095,
        stim=="d" & distance==9 ~ -0.20,
        stim=="d" & distance==14 ~ -0.26,
        T~m
      ))
    print(mu_avgs_fully_collapsed_w_hdis %>%
            arrange(distance))
  }
  
  # DO SIMULATION
  sim <- tibble()
  
  for(d in dists){
    print(d)
    mu_tmp <- mu_avgs_fully_collapsed_w_hdis %>%
      filter(distance==d) 
    
    # vector of all mu
    mu_tmp1 <- c(mu_tmp %>%
                   filter(stim=="t") %>%
                   pull(m),
                 mu_tmp %>%
                   filter(stim=="c") %>%
                   pull(m),
                 mu_tmp %>%
                   filter(stim=="d") %>%
                   pull(m))
    
    # draw from mv gaussian
    x <- rmvnorm(N,mu_tmp1,cv)
    
    # find which max and min
    mx <- apply(x, 1, which.max)
    mn <- apply(x, 1, which.min)
    
    # find choice probs for max and min
    pmax <- c(sum(mx==1)/N,
              sum(mx==2)/N,
              sum(mx==3)/N)
    pmin <- c(sum(mn==1)/N,
              sum(mn==2)/N,
              sum(mn==3)/N)
    
    # bind to tibble
    # not efficient but with only 4 distance conditions it's fine
    sim <- bind_rows(sim,
                     tibble(option=c("t","c","d"),
                            distance=d,
                            p=pmax,
                            type="best"),
                     tibble(option=c("t","c","d"),
                            distance=d,
                            p=pmin,
                            type="worst"))
  }
  sim <- sim %>%
    pivot_wider(names_from = type, 
                values_from = p)
  return(sim)
}

# FUNCTION TO ANALYZE DATA FROM ACTUAL BW EXPERIMENT ==================================================
analyze_data <- function(){
  # read in data, only include critical trials, recode bw cond to be more understandable
  d <- here("data","cleaned","bw_all.csv") %>%
    read_csv() %>%
    mutate(bw_cond=recode(bw_cond,"bw"="best-worst","wb"="worst-best")) %>%
    filter(str_detect(effect,"attraction")) %>%
    rowwise() %>%
    mutate(min=which.min(c(a1,a2,a3))) %>%
    ungroup() %>%
    mutate(best=case_when( # figure out which option they chose. if it has the min area it was a decoy choice
      # OTHER WISE NEED TO FIGURE OUT SYSTEMATICALLY IF IT WAS H OR W
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
    # NOW LABEL T / C. if an h set it's t, if it's a w set it's c
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
    ))
  
  # NOW GET MEANS
  d %>%
    pivot_longer(contains("att"),names_to = "type",values_to = "option") %>%
    mutate(type=str_remove(type,"_att")) %>%
    group_by(sub_n,distance,type,option) %>%
    summarise(n=n()) %>%
    group_by(sub_n,distance,type) %>%
    mutate(prop=n/sum(n)) %>%
    ungroup() %>%
    select(-n) %>%
    left_join(distinct(d,sub_n,bw_cond)) %>%
    group_by(distance,type,option) %>%
    summarise(m_prop=mean(prop),
              n=n(),
              se=sd(prop)/sqrt(n),
              lower=m_prop-qt(.975,n-1)*se,
              upper=m_prop+qt(.975,n-1)*se) %>%
    ungroup() %>%
    select(-c(n,se)) %>%
    pivot_wider(names_from = type, values_from = c(m_prop,lower,upper))
}

model_output_file <- here("analysis","sim_from_bayes_circle_area",glue("bw_preds_{which_model}{dec}.csv",dec=ifelse(vary_decoy_means,"_vary_decoy_means","")))
if(file_exists(model_output_file)){
  model_sims <- read_csv(model_output_file) 
}else{
  model_sims <- load_and_run_model(N,model_directory,vary_decoy_means=vary_decoy_means) # only simming triangle cond
  write_csv(model_sims,file=model_output_file)
}

# PLOT RESULTS =============================================================================
if(plot_data){
  data <- analyze_data()
}

# INIT FILE
pdf(here("analysis",
         "sim_from_bayes_circle_area",
         glue("bw_preds_v_data_{which_model}{dec}{dat}.pdf",
              dec=ifelse(vary_decoy_means,"_vary_decoy_means",""),
              dat=ifelse(plot_data,"_with_data",""))),width=8,height=8)


par(mfrow=c(2,2),pty='s',mar = rep(4.5,4))
for(d in c(2,5,9,14)){
  plot(NA,NA,xlim=c(0,.8),ylim=c(0,.8),
       xlab="p(worst)",ylab="p(best)",
       main=glue("{d}% TDD"),
       cex.main=1.75,cex.axis=1.75,cex.lab=1.75)
  if(d==2){
    if(plot_data){
      # LEGEND FOR DATA NEEDS TO BE BIGGER
      rect(.05,.5,.4,.8)
      points(.1,.75,col='red',pch=0)
      text(.2,.75,'target',cex=1.25)
      points(.1,.7,col='blue',pch=2)
      text(.26,.7,'competitor',cex=1.25)
      points(.1,.65,col='darkgreen',pch=1)
      text(.2,.65,'decoy',cex=1.25)
      
      points(.1,.58,pch=1)
      text(.2,.58,'data',cex=1.25)
      points(.1,.53,pch=16)
      text(.2,.53,'mtcm',cex=1.25)
    }else{
      rect(.05,.6,.4,.8)
      points(.1,.75,col='red',pch=15)
      text(.2,.75,'target',cex=1.25)
      points(.1,.7,col='blue',pch=17)
      text(.26,.7,'competitor',cex=1.25)
      points(.1,.65,col='darkgreen',pch=19)
      text(.2,.65,'decoy',cex=1.25)
    }
    
  }
  m_tmp <- model_sims %>%
    filter(distance==d)
  # DATA ARE OPEN POINTS, MODEL POINTS ARE FILLED
  # TARGET RED COMPETITOR BLUE DECOY GREEN
  if(plot_data){
    d_tmp <- data %>%
      filter(distance==d)
    points(d_tmp[d_tmp$option=="t",]$m_prop_worst,
           d_tmp[d_tmp$option=="t",]$m_prop_best,
           col='red',pch=0,cex=1.5)
    points(d_tmp[d_tmp$option=="c",]$m_prop_worst,
           d_tmp[d_tmp$option=="c",]$m_prop_best,
           col='blue',pch=2,cex=1.5)
    points(d_tmp[d_tmp$option=="d",]$m_prop_worst,
           d_tmp[d_tmp$option=="d",]$m_prop_best,
           col='darkgreen',pch=1,cex=1.5)
  }
  
  points(m_tmp[m_tmp$option=="t",]$worst,
         m_tmp[m_tmp$option=="t",]$best,
         col='red',pch=15,cex=1.5)
  points(m_tmp[m_tmp$option=="c",]$worst,
         m_tmp[m_tmp$option=="c",]$best,
         col='blue',pch=17,cex=1.5)
  points(m_tmp[m_tmp$option=="d",]$worst,
         m_tmp[m_tmp$option=="d",]$best,
         col='darkgreen',pch=19,cex=1.5)
}
dev.off()
  