# analyzing catch trials 
rm(list=ls())
library(here)
library(tidyverse)
library(fs)

# Setup ====================================================================================
d <- here("data","cleaned","bw_all.csv") %>%
  read_csv() %>%
  mutate(bw_cond=recode(bw_cond,"bw"="best-worst","wb"="worst-best")) %>%
  filter(str_detect(effect,"catch"))

# prop correct ================================================================================
prop_corr_by_cond_sub <- d %>%
  group_by(sub_n,effect) %>%
  summarise(prop_best=mean(choice_best_correct),
            prop_worst=mean(choice_worst_correct),
            prop_both=mean(both_correct)) %>%
  ungroup() %>%
  left_join(distinct(d,sub_n,bw_cond))

m_prop_corr_by_cond <- prop_corr_by_cond_sub %>%
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
  geom_point(data=left_join(pivot_longer(prop_corr_by_cond_sub,contains("prop")), # join to subject correct
                            distinct(d,sub_n,bw_cond)),
             aes(name,value),alpha=.1,shape=20,
             position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin=ci_lower,ymax=ci_upper),position=position_dodge(.9),width=.25)+
  ggsci::scale_fill_d3(name="condition",labels=c("best-worst","worst-best"))+
  scale_x_discrete(labels=c("prop_best"="best correct","prop_worst"="worst correct","prop_both"="both correct"))+
  facet_grid(effect~.)+
  labs(x="which choice",y="mean proportion")+ 
  ggthemes::theme_few()
ggsave(filename=here("analysis","plots","catch_prop_correct_means.jpeg"),width=6,height=4)

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
ggsave(filename=here("analysis","plots","catch_rt_hists.jpeg"),width=6,height=5)
