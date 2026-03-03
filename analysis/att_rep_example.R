rm(list=ls())

tibble(
  pb=c(.6,.3,.1),
  pw=c(.1,.2,.7),
  option=factor(c("T","C","D"),levels=c("T","C","D"))
) %>%
  ggplot(aes(pw,pb))+
  geom_text(aes(col=option,label=option),show.legend=F)+
  labs(x="p(worst)",y="p(best)")+
  ggsci::scale_color_startrek()+
  coord_fixed(xlim=c(0,.8),ylim=c(0,.8))+
  ggthemes::theme_few()
ggsave(filename = here("analysis/plots/att_ex.jpeg"),width=3,height=3)


tibble(
  pb=c(.3,.6,.1),
  pw=c(.2,.1,.7),
  option=factor(c("T","C","D"),levels=c("T","C","D"))
) %>%
  ggplot(aes(pw,pb))+
  geom_text(aes(col=option,label=option),show.legend=F)+
  labs(x="p(worst)",y="p(best)")+
  coord_fixed(xlim=c(0,.8),ylim=c(0,.8))+
  ggsci::scale_color_startrek()+
  ggthemes::theme_few()
ggsave(filename = here("analysis/plots/rep_ex.jpeg"),width=3,height=3)
