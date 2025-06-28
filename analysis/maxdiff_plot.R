rm(list=ls())
library(tidyverse)
library(here)
# playing around with the maxdiff model to understand how it works

maxdiff <- function(U){
  # U <- c(9,8,7) # t c d 
  # browser()
  N <- length(U)
  Ud <- outer(U,U,"-")
  p_bw <- matrix(NA,N,N)
  for(i in 1:N){
    for(j in 1:N){
      if(i==j){
        p_bw[i,j] <- 0
      }else{
        p_bw[i,j] <- exp(Ud[i,j])/sum(exp(Ud[Ud!=0]))
      }
    }
  }
  
  return(p_bw)
}

U <- cbind(seq(0,1,.01),rev(seq(0,1,.01)))
