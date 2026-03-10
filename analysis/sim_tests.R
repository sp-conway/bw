rm(list=ls())
library(tidyverse)
library(mvtnorm)
library(here)
library(latex2exp)
mu_T <- mu_C <- 1
mu_D <- .75
s <- c(1,1,1)
a <- ( s %*% t(s) )
mu <- c(mu_T, mu_C, mu_D)
rho_TC <- .1
rho_CD <- .1
N <- 100000
sample <- function(N, mu, cv, rho_TD){
  X <- rmvnorm(N, mu, cv)
  do_order <- function(x){
    y <- str_replace_all(as.character(order(x,decreasing = T)),
                         c("1"="t",
                           "2"="c",
                           "3"="d"))
    z <- str_flatten(y)
    return(z)
  }
  tibble(order=apply(X, 1, do_order)) %>%
    group_by(order) %>%
    summarise(N=n()) %>%
    ungroup() %>%
    mutate(prop=N/sum(N),
           rho_TD=rho_TD) %>%
    select(-N) %>%
    arrange(desc(prop))
}

sim <- vector("list")
i <- 1
for(rho_TD in seq(.11,.99,.01)){
  print(rho_TD)
  cv <- matrix(c(1, rho_TC, rho_TD,
                 rho_TC, 1, rho_CD,
                 rho_TD, rho_CD, 1), nrow=3, ncol=3,byrow=T)*a
  sim[[i]] <- sample(N, mu, cv, rho_TD)
  i <- i+1
}

sim_clean <- list_rbind(sim)
sim_clean %>%
  pivot_wider(names_from = order, values_from = prop, values_fill = 0) %>%
  pivot_longer(-rho_TD, values_to = "prop", names_to = "order") %>%
  mutate(order=factor(order, levels=c("ctd","tdc","cdt",
                                      "tcd","dtc","dct"))) %>%
  ggplot(aes(rho_TD, prop, col=order))+
  geom_path()+
  labs(x=TeX("$\\rho_{TD}$"),
       y="choice proportion",
       caption=TeX(str_c("$\\rho_{TC}=\\rho_{CD}=",rho_TC)))+
  coord_fixed(xlim=c(.1,1),
              ylim=c(0,.5))+
  scale_x_continuous(breaks=seq(.1,1,.1))+
  ggsci::scale_color_startrek()+
  ggthemes::theme_few()+
  theme(legend.direction = "horizontal",
        legend.position = "inside",
        legend.position.inside = c(.3,.8),plot.caption=element_text(hjust=0))
ggsave(filename=here("analysis/plots/mnvorm_preds_vary_rho_td.jpeg"),width=6,height=6)


