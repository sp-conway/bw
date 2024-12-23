# create attraction effect stimuli for bw exp
# exact same code as for circle experiment
# catch & fillers are created on the fly during experiment
# note that some needed very careful to tweaking to round to correct distances
# that's why percentages change from diagonal to diagonal
# fall 2024
# setup ======================================================================
# clear environment
rm(list=ls())

# libraries
library(tidyverse)
library(here)

# CONTROLS
save_stim <- T # whether to write stim values to a file

# necessary functions
distance <- function(a,b) sqrt( (a[1]-b[1])^2 + (a[2]-b[2])^2  )
compute_tdd <- function(stim,area) round(1-prod(stim)/area,digits=2)*100

# make stim ======================================================================

# three diagonals
a_3 <- c(120,195)
b_3 <- rev(a_3)

a_2 <- c(90,165)
b_2 <- rev(a_2)

a_1 <- c(60,135)
b_1 <- rev(a_1)

# distances
distance(a_1,b_1)
distance(a_2,b_2)
distance(a_3,b_3)

# areas of each diagonal
area_1 <- prod(a_1)
area_2 <- prod(a_2)
area_3 <- prod(a_3)

# lower diagonal (1) decoys
da_1_02 <- c(round(a_1[1]*.99),round(a_1[2]*.99))
da_1_05 <- c(round(a_1[1]*.975),round(a_1[2]*.975))
da_1_09 <- c(round(a_1[1]*.955),round(a_1[2]*.955))
da_1_14 <- c(round(a_1[1]*.9275),round(a_1[2]*.9275))
db_1_02 <- rev(da_1_02)
db_1_05 <- rev(da_1_05)
db_1_09 <- rev(da_1_09)
db_1_14 <- rev(da_1_14)

# tdd
compute_tdd(da_1_02, area_1)
compute_tdd(da_1_05, area_1)
compute_tdd(da_1_09, area_1)
compute_tdd(da_1_14, area_1)
compute_tdd(db_1_02, area_1)
compute_tdd(db_1_05, area_1)
compute_tdd(db_1_09, area_1)
compute_tdd(db_1_14, area_1)

# middle diagonal (2) decoys
da_2_02 <- c(round(a_2[1]*.99),round(a_2[2]*.99))
da_2_05 <- c(round(a_2[1]*.975),round(a_2[2]*.975))
da_2_09 <- c(round(a_2[1]*.9525),round(a_2[2]*.9525))
da_2_14 <- c(round(a_2[1]*.9275),round(a_2[2]*.9275))
db_2_02 <- rev(da_2_02)
db_2_05 <- rev(da_2_05)
db_2_09 <- rev(da_2_09)
db_2_14 <- rev(da_2_14)


# tdd
compute_tdd(da_2_02, area_2)
compute_tdd(da_2_05, area_2)
compute_tdd(da_2_09, area_2)
compute_tdd(da_2_14, area_2)
compute_tdd(db_2_02, area_2)
compute_tdd(db_2_05, area_2)
compute_tdd(db_2_09, area_2)
compute_tdd(db_2_14, area_2)

# upper diagonal (3) decoys
da_3_02 <- c(round(a_3[1]*.99),round(a_3[2]*.99))
da_3_05 <- c(round(a_3[1]*.975),round(a_3[2]*.975))
da_3_09 <- c(round(a_3[1]*.9525),round(a_3[2]*.9525))
da_3_14 <- c(round(a_3[1]*.9275),round(a_3[2]*.9275))
db_3_02 <- rev(da_3_02)
db_3_05 <- rev(da_3_05)
db_3_09 <- rev(da_3_09)
db_3_14 <- rev(da_3_14)

# tdd
compute_tdd(da_3_02, area_3)
compute_tdd(da_3_05, area_3)
compute_tdd(da_3_09, area_3)
compute_tdd(da_3_14, area_3)
compute_tdd(db_3_02, area_3)
compute_tdd(db_3_05, area_3)
compute_tdd(db_3_09, area_3)
compute_tdd(db_3_14, area_3)

# combine stim ======================================================================
# combine all stim in tmp matrix
stim_tmp <- rbind(
  "a_1_0"=a_1,
  "b_1_0"=b_1,
  "da_1_02"=da_1_02,
  "da_1_05"=da_1_05,
  "da_1_09"=da_1_09,
  "da_1_14"=da_1_14,
  "db_1_02"=db_1_02,
  "db_1_05"=db_1_05,
  "db_1_09"=db_1_09,
  "db_1_14"=db_1_14,
  "a_2_0"=a_2,
  "b_2_0"=b_2,
  "da_2_02"=da_2_02,
  "da_2_05"=da_2_05,
  "da_2_09"=da_2_09,
  "da_2_14"=da_2_14,
  "db_2_02"=db_2_02,
  "db_2_05"=db_2_05,
  "db_2_09"=db_2_09,
  "db_2_14"=db_2_14,
  "a_3_0"=a_3,
  "b_3_0"=b_3,
  "da_3_02"=da_3_02,
  "da_3_05"=da_3_05,
  "da_3_09"=da_3_09,
  "da_3_14"=da_3_14,
  "db_3_02"=db_3_02,
  "db_3_05"=db_3_05,
  "db_3_09"=db_3_09,
  "db_3_14"=db_3_14
)

# combine into nice df
all_stim <- tibble(
  name=rownames(stim_tmp),
  w=unname(stim_tmp[,1]),
  h=unname(stim_tmp[,2])
) %>%
  separate(name,into=c("name","diag","distance")) %>%
  mutate(across(c(diag,distance),as.numeric),
         distance=na_if(distance,0))

# plotting ======================================================================
p <- all_stim %>%
  mutate(name=toupper(name)) %>%
  ggplot(aes(w,h,col=name))+
  geom_point(size=2.5,alpha=.5)+
  labs(x="Width",y="Height")+
  scale_color_discrete(name="Stimulus")+
  scale_x_continuous(breaks=c(0,150,300))+
  scale_y_continuous(breaks=c(0,150,300))+
  coord_fixed(xlim=c(0,300),ylim=c(0,300))+
  ggthemes::theme_few()
p

# saving stim ======================================================================
if(save_stim){
  ggsave(p,filename=here("specs","stim.jpg"),width=4,height=4)
  write_csv(all_stim,file=here("specs","all_stim.csv"))
}

# catch stim ========================================================================
d1_w <- seq(a_1[1],b_1[1])
d1_h <- lm(a_1~b_1)$coefficients["(Intercept)"]-d1_w
d2_w <- seq(a_2[1],b_2[1])
d2_h <- lm(a_2~b_2)$coefficients["(Intercept)"]-d2_w
d3_w <- seq(a_3[1],b_3[1])
d3_h <- lm(a_3~b_3)$coefficients["(Intercept)"]-d3_w

par(pty='s')
plot(NA,NA,xlim=c(0,300),ylim=c(0,300))
points(d1_w,d1_h)
points(d2_w,d2_h)
points(d3_w,d3_h)

# catch stim ========================================================================
min(all_stim$w)
max(all_stim$w)
min(all_stim$h)
max(all_stim$h)

w <- runif(10000, 56, 195)
h <- runif(10000, 56, 195)
par(pty='s')
plot(NA,NA,xlim=c(0,300),ylim=c(0,300))
points(w,h,pch=".")
a <- w*h
hist(a)
