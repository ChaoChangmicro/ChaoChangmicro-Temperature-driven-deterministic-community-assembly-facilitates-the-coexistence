library(tidyverse)
library(gapminder)
library(ggsci)
library(ggprism)
library(rstatix)
library(ggpubr)

df <- read.csv("niche_breadth.csv",sep=",",header=T)
df$group <-factor(df$group,levels = c("low","high"))
is.na(df)
df1 <- df[!is.infinite(df$niche),]
df_p_val1 <- df1%>% 
  wilcox_test(niche ~ group) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "group", dodge =0.9)

df1%>%
  ggplot(aes(group,niche)) +
  stat_boxplot(geom="errorbar",position=position_dodge(width=0.2),width=0.1)+
  geom_boxplot(position=position_dodge(width =0.1),width=0.5)+
  geom_jitter(aes(fill=group),pch=21,width=0.25)+
  #geom_line(aes(group=paired),position = position_dodge(0.2),color="grey80") +
  # geom_point(aes(fill=group,group=group,size=MST.ij.ruzicka,alpha=MST.ij.ruzicka),pch=21,size=4,
  #            position = position_dodge(0.2))+
  stat_pvalue_manual(df_p_val1,label = "p.adj",y.position=24,label.size=6,hide.ns = F)+
  scale_size_continuous(range=c(1,3))+
  #facet_wrap(.~continent,nrow=1)+
  scale_fill_calc()+
  scale_x_discrete(guide = "prism_bracket")+
  scale_y_continuous(limits = c(0,25),guide = "prism_offset_minor")+
  labs(x=NULL,y=NULL)+
  theme_prism(base_line_size =0.5)+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        axis.line = element_line(color = "black",size = 0.4),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2,color = "#e5e5e5"),
        axis.text.y = element_text(color="black",size=15),
        axis.text.x = element_text(margin = margin(t = -5),color="black",size=15),
        legend.position = "none",
        panel.spacing = unit(0,"lines"))+
  coord_cartesian()
