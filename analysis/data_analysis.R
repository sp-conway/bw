rm(list=ls())
library(here)
library(tidyverse)
library(fs)

d <- here("data","cleaned","bw_all.csv") %>%
  read_csv() %>%
  mutate(bw_cond=recode(bw_cond,"bw"="best-worst","wb"="worst-best"))
d_exp <- d %>%
  filter(str_detect(effect,"practice",negate=T))
d_att <- d_exp %>%
  filter(str_detect(effect,"attraction"))
d_exp %>% 
  group_by(sub_n,effect) %>% 
  summarise(n=n())
d_exp %>%
  distinct(sub_n,bw_cond) %>%
  group_by(bw_cond) %>%
  summarise(n=n())
# prop correct ========================================================
d_corr <- d_exp %>%
  mutate(
    a_max=apply(cbind(d_exp$a1,d_exp$a2,d_exp$a3),1,max),
    a_min=apply(cbind(d_exp$a1,d_exp$a2,d_exp$a3),1,min),
    choice_best_correct=case_when(
      choice_best==1 & a1==a_max~1,
      choice_best==2 & a2==a_max~1,
      choice_best==3 & a3==a_max~1,
      T~0
    ),
    choice_worst_correct=case_when(
      choice_worst==1 & a1==a_min~1,
      choice_worst==2 & a2==a_min~1,
      choice_worst==3 & a3==a_min~1,
      T~0
    ),
    both_correct=as.integer(choice_worst_correct+choice_best_correct==2)
  )
d_corr_props <- d_corr %>%
  group_by(sub_n,effect) %>%
  summarise(prop_best=mean(choice_best_correct),
            prop_worst=mean(choice_worst_correct),
            prop_both=mean(both_correct)) %>%
  ungroup()
d_corr_props %>%
  pivot_longer(contains("prop"),names_to = "choice",values_to = "prop") %>%
  ggplot(aes(prop))+
  geom_histogram(fill="lightblue",col="black")+
  facet_grid(effect~choice,scales="free_y")+
  ggthemes::theme_few()+
  theme(text=element_text(size=15))
ggsave(filename=here("analysis","plots","prop_correct_hist.jpeg"),width=5,height=4)
d_corr_props_dist <- d_corr %>%
  filter(effect=="attraction") %>%
  group_by(sub_n,distance) %>%
  summarise(prop_best=mean(choice_best_correct),
            prop_worst=mean(choice_worst_correct),
            prop_both=mean(both_correct)) %>%
  ungroup()
d_corr_props_dist %>%
  pivot_longer(contains("prop"),names_to = "choice",values_to = "prop") %>%
  ggplot(aes(prop))+
  geom_histogram(fill="lightblue",col="black")+
  facet_grid(distance~choice,scales="free_y")+
  ggthemes::theme_few()
ggsave(filename=here("analysis","plots","prop_correct_attraction_hist.jpeg"),width=5,height=4)

m_prop_corr_by_cond <- d_corr %>%
  group_by(sub_n,effect) %>%
  summarise(prop_best=mean(choice_best_correct),
            prop_worst=mean(choice_worst_correct),
            prop_both=mean(both_correct)) %>%
  ungroup() %>%
  left_join(distinct(d_exp,sub_n,bw_cond)) %>%
  pivot_longer(contains("prop"),values_to = "prop") %>%
  group_by(name,bw_cond,effect) %>%
  summarise(m=mean(prop),
            ci_lower=m-qt(.975,n()-1)*(sd(prop)/sqrt(n())),
            ci_upper=m+qt(.975,n()-1)*(sd(prop)/sqrt(n()))) %>%
  ungroup() %>%
  mutate(name=factor(name,levels=c("prop_best","prop_worst","prop_both")))

m_prop_corr_by_cond %>%
  mutate(name=factor(name,levels=c("prop_best","prop_worst","prop_both"))) %>%
  ggplot(aes(name,m,fill=bw_cond))+
  geom_col(position=position_dodge(width=.9),width=.75,)+
  geom_point(data=left_join(pivot_longer(d_corr_props,contains("prop")),
                            distinct(d_exp,sub_n,bw_cond)),
             aes(name,value),alpha=.1,shape=20,
             position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin=ci_lower,ymax=ci_upper),position=position_dodge(.9),width=.25)+
  ggsci::scale_fill_d3(name="condition",labels=c("best-worst","worst-best"))+
  scale_x_discrete(labels=c("prop_best"="best correct","prop_worst"="worst correct","prop_both"="both correct"))+
  facet_grid(effect~.)+
  labs(x="which choice",y="mean proportion")+ 
  ggthemes::theme_few()
ggsave(filename=here("analysis","plots","prop_correct_means.jpeg"),width=6,height=4)

d_corr_att_by_dist <- d_corr %>%
  filter(effect=="attraction") %>%
  group_by(sub_n,distance) %>%
  summarise(prop_corr_best=mean(choice_best_correct),
            prop_corr_worst=mean(choice_worst_correct),
            prop_corr_both=mean(both_correct)) %>%
  ungroup() %>%
  left_join(distinct(d_corr,sub_n,bw_cond)) %>%
  pivot_longer(contains("prop_corr"),names_to = "choice",values_to = "prop") %>%
  mutate(choice=str_remove(choice,"prop_corr_"))
d_corr_att_by_dist %>%
  mutate(choice=factor(choice,levels=c("best","worst","both"))) %>%
  group_by(bw_cond,distance,choice) %>%
  summarise(m=mean(prop),
            ci_lower=m-qt(.975,n()-1)*(sd(prop)/sqrt(n())),
            ci_upper=m+qt(.975,n()-1)*(sd(prop)/sqrt(n()))) %>%
  ungroup() %>%
  ggplot(aes(distance,m,fill=choice))+
  geom_col(position=position_dodge(width=.9),width=.75,)+
  geom_point(data=d_corr_att_by_dist,
             aes(distance,prop),alpha=.4,shape=".",
             position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin=ci_lower,ymax=ci_upper),position=position_dodge(.9),width=.25)+
  ggsci::scale_fill_d3(name="")+
  scale_x_continuous(breaks=c(2,5,9,14),labels = function(x) paste0(scales::comma(x), "%")) +
  # scale_x_discrete(labels=c("prop_best"="best correct","prop_worst"="worst correct","prop_both"="both correct"))+
  facet_grid(bw_cond~.)+
  labs(x="target-decoy distance",y="mean proportion")+ 
  ggthemes::theme_few()
ggsave(filename=here("analysis","plots","att_prop_correct_means.jpeg"),width=6,height=4)


# rts ================================================================================
d_exp %>%
  pivot_longer(contains("rt")) %>%
  mutate(name=str_remove(name,"rt_")) %>%
  filter(value<10000) %>%
  ggplot(aes(value,fill=name))+
  geom_histogram(position="identity",alpha=.65)+
  facet_grid(effect~bw_cond,scales="free_y")+
  ggsci::scale_fill_cosmic(name="choice")+
  labs(x="RT (ms)")+
  ggthemes::theme_few()+
  theme(legend.position="inside",
        legend.key.size = unit(5,units = "mm"),
        legend.position.inside = c(.38,.9))
ggsave(filename=here("analysis","plots","rt_hists_all.jpeg"),width=6,height=5)

d_att %>%
  pivot_longer(contains("rt")) %>%
  mutate(name=str_remove(name,"rt_")) %>%
  filter(value<10000) %>%
  ggplot(aes(value,fill=name))+
  geom_histogram(position="identity",alpha=.65)+
  facet_grid(distance~bw_cond,scales="free_y")+
  ggsci::scale_fill_cosmic(name="choice")+
  labs(x="RT (ms)")+
  ggthemes::theme_few()+
  theme(legend.position="inside",
        legend.key.size = unit(5,units = "mm"),
        legend.position.inside = c(.38,.9))
ggsave(filename=here("analysis","plots","rt_hists_attraction.jpeg"),width=6,height=5)


d_att %>%
  group_by(distance,bw_cond) %>%
  summarise(mrt_best=mean(rt_best),
            mrst_wort=mean(rt_worst)) %>%
  ungroup() %>%
  arrange(bw_cond)
# attraction choice props ========================================================
d_att_choice <- d_att %>%
  rowwise() %>%
  mutate(min=which.min(c(a1,a2,a3))) %>%
  ungroup() %>%
  mutate(best=case_when(
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

# quick check to see if bw conditions are different

d_att_mean_by_bw_cond <- d_att_choice %>%
  select(-c(diag,h1,w1,h2,w2,h3,w3,a1,a2,a3,rt_best,rt_worst,choice_best,choice_worst,min,best,worst)) %>%
  pivot_longer(contains("att"),names_to = "type",values_to = "option") %>%
  mutate(type=str_remove(type,"_att")) %>%
  group_by(sub_n,set,distance,type,option) %>%
  summarise(n=n()) %>%
  group_by(sub_n,set,distance,type) %>%
  mutate(prop=n/sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  left_join(distinct(d_att_choice,sub_n,bw_cond)) %>%
  group_by(bw_cond,set,distance,type,option) %>%
  summarise(m_prop=mean(prop),
         n=n(),
         se=sd(prop)/sqrt(n),
         ci_lower=m_prop-qt(.975,n-1)*se,
         ci_upper=m_prop-qt(.975,n-1)*se) %>%
  ungroup() %>%
  select(-c(n,se)) %>%
  pivot_wider(names_from = type,
              values_from = c(m_prop,ci_lower,ci_upper))
d_att_mean_by_bw_cond %>%
  ggplot(aes(m_prop_worst,m_prop_best,col=option))+
  geom_point()+
  geom_errorbar(aes(ymin = ci_lower_best,ymax=ci_upper_best))+
  geom_errorbarh(aes(xmin = ci_lower_worst,xmax=ci_upper_worst))+
  coord_fixed(xlim = c(0,1),ylim=c(0,1))+
  facet_grid(bw_cond~.)+
  ggthemes::theme_few()

compute_mean_choice_props <- function(dd,groups){
  dd %>%
    select(-c(diag,h1,w1,h2,w2,h3,w3,a1,a2,a3,rt_best,rt_worst,choice_best,choice_worst,min,best,worst)) %>%
    pivot_longer(contains("att"),names_to = "type",values_to = "option") %>%
    mutate(type=str_remove(type,"_att")) %>%
    group_by(!!!groups,option,sub_n) %>%
    summarise(n=n()) %>%
    group_by(!!!groups,sub_n) %>%
    mutate(prop=n/sum(n)) %>%
    ungroup() %>%
    select(-n) %>%
    group_by(!!!groups,option) %>%
    summarise(m_prop=mean(prop),
              n=n(),
              se=sd(prop)/sqrt(n),
              ci_lower=m_prop-qt(.975,n-1)*se,
              ci_upper=m_prop+qt(.975,n-1)*se) %>%
    ungroup() %>%
    select(-c(n,se)) %>%
    pivot_wider(names_from = type,
                values_from = c(m_prop,ci_lower,ci_upper))
}

d_att_choice_mean_by_set_dist <- compute_mean_choice_props(d_att_choice,groups=vars(set,distance,type))

ggplot(d_att_choice_mean_by_set_dist,aes(m_prop_worst,m_prop_best))+
  geom_text(aes(label=option), col="darkblue",alpha=.5, size=3)+
  coord_fixed(xlim=c(0,1),ylim=c(0,1))+
  scale_x_continuous(breaks=c(0,.5,1)) +  
  scale_y_continuous(breaks=c(0,.5,1))  +  
  labs(x="p(worst)",y="p(best)",title = "data")+
  ggthemes::theme_few()+
  facet_grid(as.factor(as.numeric(distance))~set)+
  theme(text = element_text(size=14),plot.title = element_text(hjust=0.5))
ggsave(filename = here("analysis","plots","att_mean_props_by_set_dist.jpeg"),width=4,height=5)

d_att_choice_mean_by_dist <- compute_mean_choice_props(d_att_choice,groups=vars(distance,type))

ggplot(d_att_choice_mean_by_dist,aes(m_prop_worst,m_prop_best,col=option))+
  # geom_text(aes(label=option), col="darkblue",alpha=.5, size=3)+
  geom_point(shape=".")+
  geom_errorbar(aes(ymin=ci_lower_best, ymax=ci_upper_best)) +
  geom_errorbarh(aes(xmin=ci_lower_worst, xmax=ci_upper_worst), size=1) +
  coord_fixed(ratio=1, xlim=c(0, 1), ylim=c(0, 1)) +  # Enforces equal scale by fixing the ratio
  # scale_x_continuous(breaks=c(.05,.75,1)) +  
  # scale_y_continuous(breaks=c(0,.5,1))  +
  labs(x="p(worst)",y="p(best)",title = "data")+
  ggthemes::theme_few()+
  facet_grid(as.factor(as.numeric(distance))~.)+
  theme(text = element_text(size=14),plot.title = element_text(hjust=0.5))
ggsave(filename = here("analysis","plots","att_mean_props_by_dist.jpeg"),width=4,height=5)

model_preds <- here("bw_sim_from_new_bayes","bw_preds_clean.csv") %>%
  read_csv()

model_data_by_dist <- model_preds %>%
  group_by(distance,choice) %>%
  summarise(m_prop_best=mean(best),
            m_prop_worst=mean(worst)) %>%
  ungroup() %>% 
  mutate(source="model",
         option=str_sub(choice,1,1)) %>%
  bind_rows(mutate(d_att_choice_mean_by_dist,source="data"))
model_data_by_dist %>%
  ggplot(aes(m_prop_worst,m_prop_best,col=source))+
  geom_text(aes(label=option),alpha=.9, size=3)+
  coord_fixed(xlim=c(0,1),ylim=c(0,1))+
  scale_x_continuous(breaks=c(0,1,1)) +  
  scale_y_continuous(breaks=c(0,1,1))  +
  labs(x="p(worst)",y="p(best)",title = "data")+
  ggsci::scale_color_cosmic(name="")+
  ggthemes::theme_few()+
  facet_grid(as.factor(as.numeric(distance))~.)+
  theme(text = element_text(size=14),plot.title = element_text(hjust=0.5))+
  theme(legend.key = element_rect(fill = "white"))
ggsave(filename = here("analysis","plots","att_mean_props_by_dist_compare_to_model.jpeg"),width=4,height=5)

model_data_by_dist_by_set <- model_preds %>%
  group_by(set,distance,choice) %>%
  summarise(m_prop_best=mean(best),
            m_prop_worst=mean(worst)) %>%
  ungroup() %>% 
  mutate(source="model",
         option=str_sub(choice,1,1)) %>%
  bind_rows(mutate(d_att_choice_mean_by_set_dist,source="data"))
model_data_by_dist_by_set %>%
  ggplot(aes(m_prop_worst,m_prop_best,col=source))+
  geom_text(aes(label=option),alpha=.9, size=3)+
  coord_fixed(xlim=c(0,1),ylim=c(0,1))+
  scale_x_continuous(breaks=c(0,1,1)) +  
  scale_y_continuous(breaks=c(0,1,1))  +
  labs(x="p(worst)",y="p(best)",title = "data")+
  ggsci::scale_color_cosmic(name="")+
  ggthemes::theme_few()+
  facet_grid(as.factor(as.numeric(distance))~set)+
  theme(text = element_text(size=14),plot.title = element_text(hjust=0.5))+
  theme(legend.key = element_rect(fill = "white"))
ggsave(filename = here("analysis","plots","att_mean_props_by_set_dist_compare_to_model.jpeg"),width=4,height=5)
