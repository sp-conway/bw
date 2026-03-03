rm(list=ls())
library(here)
library(tidyverse)
library(fs)

# Setup ====================================================================================
d <- here("data","cleaned","bw_all.csv") %>%
  read_csv() %>%
  mutate(bw_cond=recode(bw_cond,"bw"="best-worst","wb"="worst-best")) %>%
  filter(str_detect(effect,"filler"))

# prop correct ================================================================================
prop_corr_by_cond_sub <- d %>%
  group_by(sub_n) %>%
  summarise(prop_best=mean(choice_best_correct),
            prop_worst=mean(choice_worst_correct),
            prop_both=mean(both_correct)) %>%
  ungroup() %>%
  left_join(distinct(d,sub_n,bw_cond))

m_prop_corr_by_cond <- prop_corr_by_cond_sub %>%
  pivot_longer(contains("prop"),values_to = "prop") %>%
  group_by(name,bw_cond) %>%
  summarise(m=mean(prop*100),
            s=sd(prop*100),
            ci_lower=m-qt(.975,n()-1)*(s/sqrt(n())),
            ci_upper=m+qt(.975,n()-1)*(s/sqrt(n()))) %>%
  ungroup() %>%
  mutate(name=factor(name,levels=c("prop_best","prop_worst","prop_both")))
m_prop_corr_by_cond

m_prop_corr_by_cond %>%
  mutate(name=factor(name,levels=c("prop_best","prop_worst","prop_both"))) %>%
  ggplot(aes(name,m,fill=bw_cond))+
  geom_col(position=position_dodge(width=.9),width=.75,)+
  geom_point(data=left_join(pivot_longer(prop_corr_by_cond_sub,contains("prop")), # join to subject correct
                            distinct(d,sub_n,bw_cond)),
             aes(name,value*100),alpha=.1,shape=20,
             position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin=ci_lower,ymax=ci_upper),position=position_dodge(.9),width=.25)+
  ggsci::scale_fill_d3(name="condition",labels=c("best-worst","worst-best"))+
  scale_x_discrete(labels=c("prop_best"="best correct","prop_worst"="worst correct","prop_both"="both correct"))+
  labs(x="",y="mean proportion")+ 
  ggthemes::theme_few()

# mean prop correct collapsed across condition 
options(pillar.sigfig=5)
d %>%
  pivot_longer(c(choice_best_correct,
                 choice_worst_correct,
                 both_correct)) %>%
  group_by(name) %>%
  summarise(m=mean(value*100),
            s=sd(value*100),
            ci_lower=m-qt(.975,n()-1)*(s/sqrt(n())),
            ci_upper=m+qt(.975,n()-1)*(s/sqrt(n())))
# rts ================================================================================
d %>%
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
