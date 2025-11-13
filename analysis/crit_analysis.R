rm(list=ls())
library(here)
library(tidyverse)
library(glue)
library(fs)
library(patchwork)

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

d_att_choice_mean_by_set_dist <- compute_mean_choice_props(d_att_choice,groups=vars(set,distance,type)) %>%
  mutate(set=case_when(
    set=="h"~"target tall",
    set=="w"~"target wide"
  ))
d_att_choice_mean_by_dist <- compute_mean_choice_props(d_att_choice,groups=vars(distance,type)) %>%
  mutate(set="collapsed")
d_att_choice_mean_all <- bind_rows(d_att_choice_mean_by_set_dist,d_att_choice_mean_by_dist) %>%
  mutate(distance=factor(str_glue("{distance}%"),
                         levels=c("2%","5%","9%","14%")),
         set=factor(set,levels=c("target wide","target tall","collapsed")),
         option=factor(option,levels=c("t","c","d")))
ggplot(d_att_choice_mean_all,aes(m_prop_worst,m_prop_best,col=option))+
  geom_text(aes(label=option),alpha=.8, size=4,show.legend = F)+
  coord_fixed(xlim=c(0,1),ylim=c(0,1))+
  scale_x_continuous(breaks=c(0,.5,1)) +  
  scale_y_continuous(breaks=c(0,.5,1))  +  
  labs(x="p(worst)",y="p(best)")+
  ggsci::scale_color_startrek()+
  ggthemes::theme_few()+
  facet_grid(distance~set)+
  theme(text = element_text(size=14),plot.title = element_text(hjust=0.5))
ggsave(filename = here("analysis","plots","crit_mean_props.jpeg"),width=8,height=7)

d_att_choice_mean_by_dist %>%
  select(-set) %>%
  write_csv(file=here("analysis/collapsed_data_for_andrew/crit_collapsed_data_means.csv"))
# ranking analysis ===============================================================
order_levels <- c("tcd","ctd","tdc","cdt","dtc","dct")
d_order <- d_att_choice %>%
  rowwise() %>%
  mutate(middle_att=setdiff(c("t","c","d"),c(best_att,worst_att)),
         order=factor(str_c(best_att,middle_att,worst_att),
                      levels=order_levels)) %>%
  ungroup()

d_sub_order <- d_order %>%
  group_by(sub_n,distance,order) %>%
  summarise(N=n()) %>%
  group_by(sub_n,distance) %>%
  mutate(prop=N/sum(N)) %>%
  ungroup()

d_m_order <- d_sub_order %>%
  group_by(distance,order) %>%
  summarise(m=mean(prop),
            s=qt(.025,n()-1,lower.tail = F)*sd(prop)/sqrt(n()),
            l=m-s,
            u=m+s)

# Getting model predictions, plotting in the same way data are plotted ===============================================================
# see sim_from_bayes_circle_area for code generating these predictions
f_preds <- here("analysis/sim_from_bayes_circle_area/bayes_circle_area/sigma_constant_comp_effect/bw_preds_ordering_sigma_constant_comp_effect_no_outliers_const_tc_vary_decoy_means.csv")
preds_cor_est_vary_decoy <- f_preds %>%
  read_csv() %>%
  select(-c(n,disp_cond)) %>%
  mutate(order=factor(order,levels=order_levels),
         model="cor estimated")
f_preds_cor_equal <- here("analysis/sim_from_bayes_circle_area/bayes_circle_area/sigma_constant_comp_effect/bw_preds_ordering_sigma_constant_comp_effect_no_outliers_const_tc_eq_cors_vary_decoy_means.csv")
preds_cor_equal <- f_preds_cor_equal %>%
  read_csv() %>%
  select(-c(n,disp_cond)) %>%
  mutate(order=factor(order,levels=order_levels),
         model="cor equal")
preds_all <- bind_rows(preds,
                       preds_cor_equal)

do_plot_model <- function(m_preds, dist){
  m_preds <- m_preds %>%
    filter(distance==dist)
  p <- m_preds %>%
    ggplot(aes(order,prop,col=model))+
    geom_point(aes(order,prop,col=model,shape=model),alpha=.85,size=2)+
    geom_line(aes(group=model),alpha=.5)+
    scale_shape_manual(name="",values=c(1,4))+
    scale_y_continuous(limits=c(0,.5))+
    ggsci::scale_color_d3(name="")+
    labs(x="choice order",
         y="mean proportion of choice order",
         title=glue("TDD: {dist}%"))+
    ggthemes::theme_few()+
    theme(plot.title=element_text(hjust=0.5,size=10),text=element_text(size=10),
          legend.position="none")
  if(dist==2){
    p <- p+
      theme(legend.position="inside",
            legend.position.inside=c(.75,.76),legend.title=element_blank(),
            legend.text = element_text(size=10),legend.margin = margin(6),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
  }else if(dist==5 | dist==9){
    p <- p+theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())
  }
  return(p)  
}
p_model <- map(list(2,5,9,14),
          do_plot_model,
          m_preds=preds_all)
p_model[[1]]+p_model[[2]]+p_model[[3]]+p_model[[4]]+plot_layout(ncol=1,axis_titles = "collect")
ggsave(filename = here("analysis/plots/crit_ordering_model_preds.jpeg"),width=4,height=8)

do_plot_data <- function(means, subs, dist, plot_model=F, model_preds=NULL){
  means <- means %>%
    filter(distance==dist)
  subs <- subs %>%
    filter(distance==dist)
  p <- means %>%
    ggplot(aes(order,m))+
    geom_line(aes(group=distance),col="black")+
    geom_errorbar(aes(ymin=l,ymax=u),position = position_dodge(width=2),width=.1,col="red3")+
    geom_line(data=subs,aes(order,prop,group=sub_n),col="grey",alpha=.04,inherit.aes=F)+
    scale_y_continuous(limits=c(0,.75),breaks=seq(0,.75,.25))+
    ggsci::scale_color_d3()+
    labs(x="choice order",y="mean proportion of choice order",title=glue("TDD: {dist}%"))+
    ggthemes::theme_few()+
    theme(plot.title=element_text(hjust=0.5,size=10),
          text=element_text(size=12))
  if(plot_model){
    model_preds_1 <- model_preds %>%
      filter(distance==2) 
    p <- p + geom_line(data=model_preds_1,aes(order,prop,group=distance),col="blue",inherit.aes=F)
  }
  if(dist==2){
    p <- p+
      annotate(geom="segment",x=4,xend=4.5,y=.6,col="black")+
      annotate(geom="segment",x=4,xend=4.5,y=.45,col="grey")+
      annotate(geom="text",x=5.18,y=.6,label="Group Means",size=2.75)+
      annotate(geom="text",x=5.1,y=.45,label="Participants",size=2.75)+
      theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
  }else if(dist==5 | dist==9){
    p <- p+
      theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
  }
  return(p)
} 
p_data <- map(list(2,5,9,14),
              do_plot_data,means=d_m_order,
              subs=d_sub_order,
              plot_model=T,
              model_preds=preds) 
p_data[[1]]+p_data[[2]]+p_data[[3]]+p_data[[4]]+plot_layout(ncol=1,axis_titles = "collect")
ggsave(filename = here("analysis/plots/crit_ordering_data.jpeg"),width=4,height=8)
