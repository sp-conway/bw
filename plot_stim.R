rm(list=ls())
library(tidyverse)
library(here)
library(glue)
library(latex2exp)

s <- here("specs/all_stim.csv") %>%
  read_csv() %>%
  mutate(name = toupper(name),
         name = recode(name,
                       "A" = "H",
                       "B" = "W",
                       "DA" = "D_H",
                       "DB" = "D_W"),
         name = factor(name, levels = c("H", "W", "D_H", "D_W")))

ggplot(s,aes(w,h,col=name,shape=name))+
  geom_point(size = 2.5, alpha = 0.75) +
  geom_abline(intercept = 195, slope = -1, linetype = "dashed", alpha = 0.3) +
  geom_abline(intercept = 255, slope = -1, linetype = "dashed", alpha = 0.3) +
  geom_abline(intercept = 315, slope = -1, linetype = "dashed", alpha = 0.3) +
  annotate(geom="text",x=95,y=95,label="lower",angle=-45,alpha=.3)+
  annotate(geom="text",x=125,y=125,label="middle",angle=-45,alpha=.3)+
  annotate(geom="text",x=155,y=155,label="upper",angle=-45,alpha=.3)+
  labs(x = "Width (px)", y = "Height (px)") +
  coord_fixed(xlim=c(50,200),ylim=c(50,200))+
  scale_shape_manual(values=c(15,0,16,1),
                     name="Stimulus",
                     labels = list(expression(H),
                                   expression(W),
                                   expression(D[H]),
                                   expression(D[W])))+
  ggsci::scale_color_lancet(
    name = "Stimulus",
    labels = list(expression(H),
                  expression(W),
                  expression(D[H]),
                  expression(D[W]))
  )+
  ggthemes::theme_few()
ggsave(filename=here("specs/stim_plot.jpeg"),width=4,height=4)
