rm(list=ls())
library(here)
library(tidyverse)
library(glue)
library(fs)

# read in data, only include critical trials, recode bw cond to be more understandable
d <- here("data","cleaned","bw_all.csv") %>%
  read_csv() %>%
  mutate(bw_cond=recode(bw_cond,"bw"="best-worst","wb"="worst-best")) %>%
  filter(str_detect(effect,"attraction"))

# n block
unique(d$block_n)
max(d$trial_n)
# quick check of subject counts
d %>%
  distinct(sub_n,bw_cond) %>%
  group_by(bw_cond) %>%
  summarise(n=n())

d %>%
  group_by(sub_n) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  arrange(desc(n))

length(unique(d$sub_n))
# prop correct ========================================================
# subject level prop correct
d_corr_props <- d %>%
  group_by(sub_n) %>%
  summarise(prop_best=mean(choice_best_correct),
            prop_worst=mean(choice_worst_correct),
            prop_both=mean(both_correct)) %>%
  ungroup()

# prop correct histogram
d_corr_props %>%
  pivot_longer(contains("prop"),names_to = "choice",values_to = "prop") %>%
  ggplot(aes(prop))+
  geom_histogram(fill="lightblue",col="black")+
  facet_grid(.~choice,scales="free_y")+
  ggthemes::theme_few()+
  theme(text=element_text(size=15))
ggsave(filename=here("analysis","plots","crit_prop_correct_hist.jpeg"),width=5,height=4)

# subject level prop correct, by distance
d_corr_props_dist <- d %>%
  filter(effect=="attraction") %>%
  group_by(sub_n,distance) %>%
  summarise(prop_best=mean(choice_best_correct),
            prop_worst=mean(choice_worst_correct),
            prop_both=mean(both_correct)) %>%
  ungroup()

# plotting histogram of subject level prop correct, split by distance
d_corr_props_dist %>%
  pivot_longer(contains("prop"),names_to = "choice",values_to = "prop") %>%
  ggplot(aes(prop))+
  geom_histogram(fill="lightblue",col="black")+
  facet_grid(distance~choice,scales="free_y")+
  ggthemes::theme_few()
ggsave(filename=here("analysis","plots","crit_prop_correct_by_dist_hist.jpeg"),width=5,height=4)

# mean proportion correct by condition
m_prop_corr_by_cond <- d_corr_props %>%
  left_join(distinct(d,sub_n,bw_cond)) %>% # need to join back to get bw cond values
  pivot_longer(contains("prop"),values_to = "prop") %>%
  group_by(name,bw_cond) %>%
  summarise(m=mean(prop),
            ci_lower=m-qt(.975,n()-1)*(sd(prop)/sqrt(n())),
            ci_upper=m+qt(.975,n()-1)*(sd(prop)/sqrt(n()))) %>%
  ungroup() %>%
  mutate(name=factor(name,levels=c("prop_best","prop_worst","prop_both")))

# plot mean prop correct by bw condition
m_prop_corr_by_cond %>%
  mutate(name=factor(name,levels=c("prop_best","prop_worst","prop_both"))) %>%
  ggplot(aes(name,m,fill=bw_cond))+
  geom_col(position=position_dodge(width=.9),width=.75,)+
  geom_point(data=left_join(pivot_longer(d_corr_props,contains("prop")),
                            distinct(d,sub_n,bw_cond)), #  joining to get indiv data points
             aes(name,value),alpha=.1,shape=20,
             position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin=ci_lower,ymax=ci_upper),position=position_dodge(.9),width=.25)+
  ggsci::scale_fill_d3(name="condition",labels=c("best-worst","worst-best"))+
  scale_x_discrete(labels=c("prop_best"="best correct","prop_worst"="worst correct","prop_both"="both correct"))+
  labs(x="which choice",y="mean proportion")+ 
  ggthemes::theme_few()
ggsave(filename=here("analysis","plots","crit_prop_correct_means.jpeg"),width=6,height=4)

# prop correct by condition, choice
# both meand and indiv subs
d_corr_by_dist <- d_corr_props_dist %>%
  left_join(distinct(d,sub_n,bw_cond)) %>%
  pivot_longer(contains("prop"),names_to = "choice",values_to = "prop") %>%
  mutate(choice=str_remove(choice,"prop_"))

# plot mean prop correct by condition and choice
d_corr_by_dist %>%
  mutate(choice=factor(choice,levels=c("best","worst","both"))) %>%
  group_by(bw_cond,distance,choice) %>%
  summarise(m=mean(prop),
            ci_lower=m-qt(.975,n()-1)*(sd(prop)/sqrt(n())),
            ci_upper=m+qt(.975,n()-1)*(sd(prop)/sqrt(n()))) %>%
  ungroup() %>%
  ggplot(aes(distance,m,fill=choice))+
  geom_col(position=position_dodge(width=.9),width=.75,)+
  geom_point(data=d_corr_by_dist,
             aes(distance,prop),alpha=.4,shape=".",
             position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin=ci_lower,ymax=ci_upper),position=position_dodge(.9),width=.25)+
  ggsci::scale_fill_d3(name="")+
  scale_x_continuous(breaks=c(2,5,9,14),labels = function(x) paste0(scales::comma(x), "%")) +
  # scale_x_discrete(labels=c("prop_best"="best correct","prop_worst"="worst correct","prop_both"="both correct"))+
  facet_grid(bw_cond~.)+
  labs(x="target-decoy distance",y="mean proportion")+ 
  ggthemes::theme_few()
ggsave(filename=here("analysis","plots","crit_prop_correct_means.jpeg"),width=6,height=4)


# rts ================================================================================
d %>%
  pivot_longer(contains("rt")) %>%
  mutate(name=str_remove(name,"rt_")) %>%
  filter(value<10000) %>%
  ggplot(aes(value,fill=name))+
  geom_histogram(position="identity",alpha=.65)+
  facet_grid(.~bw_cond,scales="free_y")+
  ggsci::scale_fill_cosmic(name="choice")+
  labs(x="RT (ms)")+
  ggthemes::theme_few()+
  theme(legend.position="inside",
        legend.key.size = unit(5,units = "mm"),
        legend.position.inside = c(.38,.9))
ggsave(filename=here("analysis","plots","crit_rt_hists.jpeg"),width=6,height=5)

d %>%
  group_by(distance,bw_cond) %>%
  summarise(mrt_best=mean(rt_best),
            mrst_wort=mean(rt_worst)) %>%
  ungroup() %>%
  arrange(bw_cond)

# critical choice props ===========================================================================
d_att_choice <- d %>%
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
# they're not, so we're collapsing over bw cond
compute_mean_choice_props <- function(dd,groups){
  dd %>%
    select(-c(h1,w1,h2,w2,h3,w3,a1,a2,a3,rt_best,rt_worst,choice_best,choice_worst,min,best,worst)) %>%
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
ggsave(filename = here("analysis","plots","crit_mean_props_by_set_dist.jpeg"),width=4,height=5)

d_att_choice_mean_by_dist <- compute_mean_choice_props(d_att_choice,groups=vars(distance,type))

ggplot(d_att_choice_mean_by_dist,aes(m_prop_worst,m_prop_best,col=option))+
  # geom_text(aes(label=option), col="darkblue",alpha=.5, size=3)+
  geom_point(alpha=.5)+
  coord_fixed(xlim=c(0,1),ylim=c(0,1))+
  scale_x_continuous(breaks=c(0,.5,1),labels = c("0",".5","1"))+
  scale_y_continuous(breaks=c(0,.5,1),labels = c("0",".5","1"))+
  facet_grid(distance~.)+
  labs(x="p(worst)",y="p(best)")+
  ggsci::scale_color_startrek()+
  ggthemes::theme_few()+
  theme(text=element_text(size=18))
ggsave(filename = here("analysis","plots","crit_mean_props_by_dist.jpeg"),width=4,height=6)



d_att_choice_mean_by_dist_diag <- compute_mean_choice_props(d_att_choice,groups=vars(distance,diag,type))

ggplot(d_att_choice_mean_by_dist_diag,aes(m_prop_worst,m_prop_best,col=option))+
  # geom_text(aes(label=option), col="darkblue",alpha=.5, size=3)+
  geom_point(shape=".")+
  geom_errorbar(aes(ymin=ci_lower_best, ymax=ci_upper_best)) +
  geom_errorbarh(aes(xmin=ci_lower_worst, xmax=ci_upper_worst), linewidth=1) +
  coord_fixed(ratio=1, xlim=c(0, 1), ylim=c(0, 1)) +  # Enforces equal scale by fixing the ratio
  # scale_x_continuous(breaks=c(.05,.75,1)) +  
  # scale_y_continuous(breaks=c(0,.5,1))  +
  labs(x="p(worst)",y="p(best)",title = "data")+
  ggthemes::theme_few()+
  facet_grid(as.factor(as.numeric(distance))~diag)+
  theme(text = element_text(size=14),plot.title = element_text(hjust=0.5))
ggsave(filename = here("analysis","plots","crit_mean_props_by_dist_diag.jpeg"),width=6,height=7)

d
# model predictions (all 3 variants)==============================================================================
outl <- "no_outliers"

# sigma constant - no target/comp effect
# sigma constant target effect and sigma constant comp effect are equivalent
# just which one is reference gets varied

for(which_model in c("sigma_constant","sigma_constant_target_effect","sigma_constant_comp_effect")){
  # get predictions, filter out horizontal bc only had triangle display here
  model_preds <-here("analysis","sim_from_bayes_circle_area","bayes_circle_area",which_model,glue("bw_preds_{which_model}_{outl}.csv")) %>%
    read_csv() %>%
    filter(disp_cond=="triangle")
  
  # bind with data
  model_data_by_dist <- model_preds %>%
    mutate(source="model") %>% 
    bind_rows(mutate(d_att_choice_mean_by_dist,source="data"))
  
  # plotting
  # plot labels, color labels by source
  model_data_by_dist %>%
    ggplot(aes(m_prop_worst,m_prop_best,shape=source,col=option))+
    geom_point(alpha=.5, size=1)+
    coord_fixed(xlim=c(0,1),ylim=c(0,1))+
    scale_x_continuous(breaks=c(0,1,1)) +  
    scale_y_continuous(breaks=c(0,1,1))  +
    labs(x="p(worst)",y="p(best)",title = "data")+
    ggsci::scale_color_startrek()+
    scale_shape_manual(values=c(1,4),name="")+
    ggthemes::theme_few()+
    facet_grid(as.factor(as.numeric(distance))~.)+
    theme(text = element_text(size=14),plot.title = element_text(hjust=0.5))+
    theme(legend.key = element_rect(fill = "white"))
  ggsave(filename = here("analysis","plots",glue("crit_mean_props_by_dist_compare_to_model_{which_model}_{outl}.jpeg")),width=4,height=5)
}

# indiv ppts ================================================================================================
d_att_choice %>%
  select(-c(diag,h1,w1,h2,w2,h3,w3,a1,a2,a3,rt_best,rt_worst,choice_best,choice_worst,min,best,worst)) %>%
  pivot_longer(contains("att"),names_to = "type",values_to = "option") %>%
  mutate(type=str_remove(type,"_att")) %>%
  group_by(sub_n,distance,type,option) %>%
  summarise(n=n()) %>%
  group_by(sub_n,distance,type) %>%
  mutate(prop=n/sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from = type,values_from = prop, values_fill = 0) %>%
  ggplot(aes(worst,best,col=option))+
  geom_point(shape=".")+
  coord_fixed(ratio=1, xlim=c(0, 1), ylim=c(0, 1)) +  
  ggsci::scale_color_startrek(name="")+
  labs(x="p(worst)",y="p(best)")+
  ggthemes::theme_few()+
  facet_grid(as.factor(as.numeric(distance))~.)+
  theme(text = element_text(size=16),plot.title = element_text(hjust=0.5))
ggsave(filename=here("analysis","plots","crit_props_by_dist_indiv.jpeg"),width=5,height=6)

# hawkins analysis - marginal choice props ========================================================================
d_att_choice_collapsed <- d_att_choice %>%
  mutate(best=case_when(
    best=="d" & set=="w"~glue("dw_{distance}_{diag}"),
    best=="d" & set=="h"~glue("dh_{distance}_{diag}"),
    best=="w" ~ glue("w_{diag}"),
    best=="h" ~ glue("h_{diag}"),
  ),
  worst=case_when(
    worst=="d" & set=="w"~glue("dw_{distance}_{diag}"),
    worst=="d" & set=="h"~glue("dh_{distance}_{diag}"),
    worst=="w" ~ glue("w_{diag}"),
    worst=="h" ~ glue("h_{diag}")
  )) %>%
  select(-c(best_att,worst_att)) %>%
  pivot_longer(c(best,worst),
               names_to = "type",
               values_to="option")

d_att_choice_collapsed_no_d <- d_att_choice_collapsed %>%
  filter(str_detect(option,"d",negate=T)) %>%
  group_by(sub_n,diag,option,type) %>%
  summarise(N=n()) %>%
  ungroup() %>%
  group_by(sub_n,diag,type) %>%
  mutate(prop=N/sum(N)) %>%
  group_by(option,type) %>%
  summarise(m_prop=mean(prop),
            n=n(),
            se=sd(prop)/sqrt(n),
            ci_lower=m_prop-qt(.975,n-1)*se,
            ci_upper=m_prop+qt(.975,n-1)*se) %>%
  ungroup() %>%
  select(-c(n,se)) %>%
  pivot_wider(names_from = type,
              values_from = c(m_prop,ci_lower,ci_upper))

d_att_choice_collapsed_with_d <- d_att_choice_collapsed %>%
  filter(str_detect(option,"d")) %>%
  group_by(sub_n,distance,diag,option,type) %>%
  summarise(N=n()) %>%
  ungroup() %>%
  group_by(sub_n,distance,diag,type) %>%
  mutate(prop=N/sum(N)) %>%
  group_by(option,type) %>%
  summarise(m_prop=mean(prop),
            n=n(),
            se=sd(prop)/sqrt(n),
            ci_lower=m_prop-qt(.975,n-1)*se,
            ci_upper=m_prop+qt(.975,n-1)*se) %>%
  ungroup() %>%
  select(-c(n,se)) %>%
  pivot_wider(names_from = type,
              values_from = c(m_prop,ci_lower,ci_upper))
d_att_choice_collapsed_means_all <- bind_rows(d_att_choice_collapsed_with_d,
                                              d_att_choice_collapsed_no_d)

d_att_choice_collapsed_means_all %>%
  ggplot(aes(m_prop_worst,m_prop_best))+
  # geom_point(col="lightblue",alpha=.5)+
  geom_errorbar(aes(ymin=ci_lower_best,ymax=ci_upper_best),alpha=.65,col="lightblue")+
  geom_errorbarh(aes(xmin=ci_lower_worst,xmax=ci_upper_worst),alpha=.65,col="lightblue")+
  coord_fixed(xlim = c(.3,.8),ylim = c(.3,.8))+
  labs(x="mean p(worst)",y="mean p(best)")+
  ggthemes::theme_few()+
  theme(text=element_text(size=16))
ggsave(filename=here("analysis","plots","crit_mean_props_marginal.jpeg"),width=5,height=5)  

# compute choice props by set, but plot in terms of h/w =========================================
d_mchoice_by_set_dist <- d_att_choice %>%
  mutate(set=case_when(
    set=="w"~"h-w-dw",
    set=="h"~"h-w-dh"
  )) %>%
  pivot_longer(c(best,worst),names_to = "type", values_to = "choice") %>%
  group_by(sub_n,distance,set,type,choice) %>%
  summarise(N=n()) %>%
  ungroup() %>%
  group_by(sub_n,distance,set,type) %>%
  mutate(prop=N/sum(N)) %>%
  group_by(distance,set,type,choice) %>%
  summarise(m=mean(prop),
            s=sd(prop),
            se=s/sqrt(n()),
            ci_lower=m-qt(.975,n()-1,lower.tail = T)*se,
            ci_upper=m+qt(.975,n()-1,lower.tail = T)*se) %>%
  ungroup() %>%
  select(-c(s,se)) %>%
  pivot_wider(names_from = type,
              values_from = c(m,ci_lower,ci_upper))
#Our transformation function
d_mchoice_by_set_dist %>%
  ggplot(aes(m_worst,m_best,col=choice,shape=set))+
  geom_point(alpha=.5)+
  coord_fixed(xlim=c(0,1),ylim=c(0,1))+
  facet_grid(distance~.)+
  labs(x="p(worst)",y="p(best)")+
  ggsci::scale_color_startrek()+
  scale_x_continuous(breaks=c(0,.5,1),labels = c("0",".5","1"))+
  scale_y_continuous(breaks=c(0,.5,1),labels = c("0",".5","1"))+
  ggthemes::theme_few()+
  theme(text=element_text(size=18))
ggsave(filename = here("analysis","plots","crit_mean_choice_by_set_dist_labelHW.jpeg"),
       width=4,height=6)
