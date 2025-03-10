rm(list=ls())

# playing around with the maxdiff model to understand how it works

do_maxdiff <- function(U){
  # U <- c(9,8,7) # t c d 
  # browser()
  Ud <- outer(U,U,"-")
  p_bw <- matrix(NA,3,3)
  for(i in 1:3){
    for(j in 1:3){
      if(i==j){
        p_bw[i,j] <- 0
      }else{
        p_bw[i,j] <- exp(Ud[i,j])/sum(exp(Ud[Ud!=0]))
      }
    }
  }
  
  p_bw[1,2] # probability target best, competitor worst, decoy middle # tcd
  p_bw[1,3] # probability target best, decoy worst, competitor middle # tdc
  p_bw[2,1] # probability competitor best, target worst, decoy middle # ctd
  p_bw[2,3] # probability competitor best, decoy worst, target middle # cdt
  p_bw[3,1] # probability decoy best, target worst, competitor middle # dtc
  p_bw[3,2] # probability decoy best, competitor worst, target middle # dct
  
  p_rank <- c(p_bw[1,2], 
              p_bw[1,3], 
              p_bw[2,1], 
              p_bw[2,3],
              p_bw[3,1],
              p_bw[3,2])
  return(p_rank)
}

do_maxdiff2 <- function(U){
  # U <- c(9,8,7) # t c d 
  # browser()
  # Ud <- outer(U,U,"-")
  # p_bw <- matrix(NA,3,3)
  # for(i in 1:3){
  #   for(j in 1:3){
  #     if(i==j){
  #       p_bw[i,j] <- 0
  #     }else{
  #       p_bw[i,j] <- exp(Ud[i,j])/sum(exp(Ud[Ud!=0]))
  #     }
  #   }
  # }
  Ub <- U
  Uw <- -U
  pb <- exp(Ub)/sum(exp(Ub))
  pw <- exp(Uw)/sum(exp(Uw))
  p_bw[1,2] <- pb[1]*pw[2]# probability target best, competitor worst, decoy middle # tcd
  p_bw[1,3] <- pb[1]*pw[3]# probability target best, decoy worst, competitor middle # tdc
  p_bw[2,1] <- pb[2]*pw[1]# probability competitor best, target worst, decoy middle # ctd
  p_bw[2,3] <- pb[2]*pw[3]# probability competitor best, decoy worst, target middle # cdt
  p_bw[3,1] <- pb[3]*pw[1]# probability decoy best, target worst, competitor middle # dtc
  p_bw[3,2] <- pb[3]*pw[2]# probability decoy best, competitor worst, target middle # dct
  
  p_rank <- c(p_bw[1,2], 
              p_bw[1,3], 
              p_bw[2,1], 
              p_bw[2,3],
              p_bw[3,1],
              p_bw[3,2])
  return(p_rank)
}

U <- c(9,8,7)
do_maxdiff(U)
do_maxdiff(U)

do_maxdiff_sep <- function(Ub, Uw){
  # U <- c(9,8,7) # t c d 
  # p_bw <- matrix(NA,3,3)
  Ud <- matrix(NA, 3, 3)
  for(i in 1:3){
    for(j in 1:3){
      if(i==j){
        Ud[i,j] <- 0
      }else{
        Ud[i,j] <- Ub[i]-Uw[j]
      }
    }
  }
  for(i in 1:3){
    for(j in 1:3){
      if(i==j){
        p_bw[i,j] <- 0
      }else{
        p_bw[i,j] <- exp(Ud[i,j])/sum(exp(Ud[Ud!=0]))
      }
    }
  }
  
  p_bw[1,2] # probability target best, competitor worst, decoy middle # tcd
  p_bw[1,3] # probability target best, decoy worst, competitor middle # tdc
  p_bw[2,1] # probability competitor best, target worst, decoy middle # ctd
  p_bw[2,3] # probability competitor best, decoy worst, target middle # cdt
  p_bw[3,1] # probability decoy best, target worst, competitor middle # dtc
  p_bw[3,2] # probability decoy best, competitor worst, target middle # dct
  
  p_rank <- c(p_bw[1,2], 
              p_bw[1,3], 
              p_bw[2,1], 
              p_bw[2,3],
              p_bw[3,1],
              p_bw[3,2])
  return(p_rank)
}


