rm(list=ls())
# simulating maxdiff model to show monotonicity
library(tidyverse)
library(here)

model <- function(index){
  U <- rnorm(3)
  Ub <- U
  Uw <- -U
  pb <- exp(Ub)/sum(exp(Ub))
  pw <- exp(Uw)/sum(exp(Uw))
  return(
    tibble(
      index=index,
      option=c(1:3),
      pb=pb,
      pw=pw
    )
  )
}

set.seed(6)
sim <- map(1:15, model) %>%
  list_rbind()

sim %>%
  ggplot(aes(pw,pb,group=index))+
  geom_point(alpha=.3,color="darkred")+
  geom_line(alpha=.6,color="darkred")+
  labs(x="p(worst)",y="p(best)")+
  coord_fixed(xlim=c(0,1),ylim=c(0,1))+
  ggthemes::theme_few()
ggsave(filename=here("analysis/plots/maxdiff_sim_monotonic.jpeg"),width=3,height=3)

