rm(list=ls())
library(here)
library(fs)
library(tidyverse)

choicef <- here("experiment_code","data","test") %>%
  dir_ls() 
choice <- map_dfr(choicef, read_csv) %>%
  mutate(a1=h1*w1,
         a2=h2*w2,
         a3=h3*w3)
  

choice_att <- choice %>%
  filter(str_detect(effect,"att")) %>%
  mutate(a1=h1*w1,
         a2=h2*w2,
         a3=h3*w3,
         d_identity=case_when(
           a1 < a2 & a1 < a3 ~ 1,
           a2 < a1 & a2 < a3 ~ 2,
           a3 < a1 & a3 < a2 ~ 3,
         ),
         t_identity = case_when(
           set=="a-b-da" & d_identity == 1 & h2 > w2 ~ 2,
           set=="a-b-da" & d_identity == 1 & h3 > w3 ~ 3,
           set=="a-b-da" & d_identity == 2 & h3 > w3 ~ 3,
           set=="a-b-da" & d_identity == 2 & h1 > w1 ~ 1,
           set=="a-b-da" & d_identity == 3 & h1 > w1 ~ 1,
           set=="a-b-da" & d_identity == 3 & h2 > w2 ~ 2,
           
           set=="a-b-db" & d_identity == 1 & h2 < w2 ~ 2,
           set=="a-b-db" & d_identity == 1 & h3 < w3 ~ 3,
           set=="a-b-db" & d_identity == 2 & h3 < w3 ~ 3,
           set=="a-b-db" & d_identity == 2 & h1 < w1 ~ 1,
           set=="a-b-db" & d_identity == 3 & h1 < w1 ~ 1,
           set=="a-b-db" & d_identity == 3 & h2 < w2 ~ 2
         ),
         c_identity=case_when(
           d_identity==1 & t_identity==2 ~ 3,
           d_identity==1 & t_identity==3 ~ 2,
           d_identity==2 & t_identity==1 ~ 3,
           d_identity==2 & t_identity==3 ~ 1,
           d_identity==3 & t_identity==1 ~ 2,
           d_identity==3 & t_identity==2 ~ 1,
         ),
         ta=case_when(
           t_identity==1~a1,
           t_identity==2~a2,
           t_identity==3~a3,
         ),
         ca=case_when(
           c_identity==1~a1,
           c_identity==2~a2,
           c_identity==3~a3,
         ),
         da=case_when(
           d_identity==1~a1,
           d_identity==2~a2,
           d_identity==3~a3,
         ),
         tdd=round(100*(1-da/ta)))
unique(choice_att$tdd)         

mina_catch <- meda_catch <- maxa_catch <- c()
catch <- filter(choice,effect=="catch")
for(i in 1:nrow(catch)){
  print(i)
  maxa_catch <- c(maxa_catch, max(c(catch$a1[i],catch$a2[i],catch$a3[i])))
  mina_catch <- c(mina_catch, min(c(catch$a1[i],catch$a2[i],catch$a3[i])))
  meda_catch <- c(meda_catch, c(catch$a1[i],catch$a2[i],catch$a3[i])[order(c(catch$a1[i],catch$a2[i],catch$a3[i]))[2]])
}

catch_a <- bind_cols(min=mina_catch,med=meda_catch,max=maxa_catch)
catch_a %>% 
  mutate(t=1:n()) %>%
  pivot_longer(-t) %>%
  ggplot(aes(value))+
  geom_histogram()+
  facet_grid(name~.)

choice %>%
  group_by(sub_n,block_n,trial_n) %>%
  mutate(t=1:n()) %>%
  distinct(choice)

choice %>%
  filter(str_detect(effect,"attraction")) %>%
  pivot_longer(c(a1,a2,a3)) %>%
  ggplot(aes(value))+
  geom_histogram()
hist(choice_att$ta)
hist(choice_att$da)

choice_att %>%
  select(sub_n,diag,distance,set,block_n,trial_n,h1,h2,h3,w1,w2,w3,t_identity,d_identity,c_identity) %>%
  distinct(sub_n,diag,distance,set,block_n,trial_n,h1,h2,h3,w1,w2,w3,t_identity,d_identity,c_identity) %>%
  pivot_longer(c(h1,h2,h3,w1,w2,w3)) %>%
  mutate(name1=str_sub(name,1,1),
         stim=str_sub(name,2,2)) %>%
  select(-name) %>%
  pivot_wider(names_from = name1,values_from = value) %>%
  mutate(name=case_when(
    stim==1 & t_identity=="1" ~ "t",
    stim==1 & c_identity=="1" ~ "c",
    stim==1 & d_identity=="1" ~ "d",
    stim==2 & t_identity=="2" ~ "t",
    stim==2 & c_identity=="2" ~ "c",
    stim==2 & d_identity=="2" ~ "d",
    stim==3 & t_identity=="3" ~ "t",
    stim==3 & c_identity=="3" ~ "c",
    stim==3 & d_identity=="3" ~ "d",
  )) %>%
  ggplot(aes(jitter(w,amount = 10),jitter(h),col=name))+
  coord_fixed()+
  geom_text(aes(label=stim),alpha=.4)+ #shouldn't have any particular label
  facet_wrap(vars(set))+
  ggthemes::theme_few()

choice %>%
  filter(effect=="filler") %>%
  distinct(sub_n,block_n,trial_n,h1,h2,h3,w1,w2,w3) %>%
  pivot_longer(c(h1,h2,h3,w1,w2,w3)) %>%
  mutate(name1=str_sub(name,1,1),
         stim=str_sub(name,2,2)) %>%
  select(-name) %>%
  pivot_wider(names_from = name1,values_from = value) %>%
  ggplot(aes(w,h,col=stim))+
  coord_fixed()+
  geom_point(alpha=.5)+
  ggthemes::theme_few()

catch %>%
  distinct(sub_n,block_n,trial_n,h1,h2,h3,w1,w2,w3) %>%
  pivot_longer(c(h1,h2,h3,w1,w2,w3)) %>%
  mutate(name1=str_sub(name,1,1),
         stim=str_sub(name,2,2)) %>%
  select(-name) %>%
  pivot_wider(names_from = name1,values_from = value) %>%
  ggplot(aes(w,h,col=stim))+
  coord_fixed()+
  geom_point()+
  ggthemes::theme_few()

# checking to make sure the same choice was never recorded for both options
choice %>%
  select(sub_n,block_n,trial_n,choice) %>%
  group_by(sub_n,block_n,trial_n) %>%
  mutate(choice_n=1:n()) %>%
  ungroup() %>%
  pivot_longer(c(choice)) %>%
  pivot_wider(names_from = choice_n,
              values_from = value,names_prefix = "choice_") %>%
  mutate(check=choice_1==choice_2) %>%
  distinct(check)

rt <- choice %>%
  select(sub_n,block_n,trial_n,rt) %>%
  group_by(sub_n,block_n,trial_n) %>%
  mutate(choice_n=1:n()) %>%
  ungroup() %>%
  pivot_longer(c(rt)) %>%
  pivot_wider(names_from = choice_n,
              values_from = value,names_prefix = "rt_") 
hist(rt$rt_1)
hist(rt$rt_2)
table(rt$rt_1 < rt$rt_2)

choice %>%
  filter(str_detect(effect,"attraction|catch")) %>%
  distinct(sub_n,block_n,trial_n,effect,h1,h2,h3,w1,w2,w3) %>%
  pivot_longer(c(h1,h2,h3,w1,w2,w3)) %>%
  mutate(name1=str_sub(name,1,1),
         stim=str_sub(name,2,2)) %>%
  select(-name) %>%
  pivot_wider(names_from = name1,values_from = value) %>%
  ggplot(aes(w,h,col=effect))+
  coord_fixed()+
  geom_point(alpha=.5)+
  ggthemes::theme_few()
