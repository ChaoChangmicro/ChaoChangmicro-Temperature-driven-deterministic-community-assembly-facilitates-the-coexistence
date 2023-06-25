library(gdm)
otutab = read.delim("./asv_all_rare.txt", row.names=1)
otutab <- t(otutab )
otutab <- as.data.frame(otutab)
otutab$site <- row.names(otutab)
env = read.delim("./gdm_env.txt")
env[4:ncol(env)] <- data.frame(scale(env[,4:12]))
exFormat1a <- formatsitepair(otutab, 1, siteColumn="site", XColumn="longitude", YColumn="latitude", predData=env)
table(is.na(exFormat1a  ))
gdm.1 <- gdm(data=exFormat1a, geo=F)
summary(gdm.1)

length(gdm.1$predictors) # get ideal of number of panels
#> [1] 11
plot(gdm.1, plot.layout=c(4,3))

# plot(gdm.1$predicted,
#      gdm.1$observed,
#      lwd=3,
#      col = "blue",
#     # xlim=c(0.65,1.15),
#      ylim=c(0,0.9),
#      pch=19,
#      cex=0.5,
#      xlab="Predicted community dissimilarity",
#      ylab="Observed community dissimilarity")
#   model<-lm(gdm.1$observed~gdm.1$ecological)
#   abline(model,lty=1,cex=1)
  
  

gdm.1.splineDat <- isplineExtract(gdm.1)
str(gdm.1.splineDat)
df <- as.data.frame(gdm.1.splineDat)
write.csv(df,file = "gdm参数.csv")

library(tidyverse)
library(ggsci)
library(patchwork)
df1 <- read.csv("gdm参数1.csv",row.names = 1,header = T,sep = ",")
gp1<-ggplot(df1, aes(x, y))+
 # geom_point(size=2)+
  geom_line(aes(group =group,color=group),size=1)+
  # geom_point(size=3,aes(color=group))+
  # geom_errorbar(aes(ymin = mean-se,ymax =mean+se,
  #                   group = group,color=group),width = 0.05)+
  scale_color_brewer(palette = "Set1")+
  labs(x = "Predictor Dissimilarity", y = "Observed community dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"),
         legend.position = c(0.4,0.92),
         legend.direction = "horizontal",
         legend.key = element_blank(),
         legend.background = element_blank(),
         legend.text = element_text(size = 10, face = "bold"),
         legend.title = element_text(size = 11, face = "bold"),
         aspect.ratio = 1)
gp1
  

xxx <- as.data.frame(gdm.1$predicted)
yyy <- as.data.frame(gdm.1$observed)
zzz <- cbind(xxx,yyy)
colnames(zzz) <- c("predicted","observed")
gp2<-ggplot( zzz,aes(predicted,observed))+
  geom_point(size=1,color="black",fill="blue",shape=21)+
  geom_smooth(method = "lm", formula = y~x, color = "black", fill = "#cbc9e2")+ #颜色选自https://colorbrewer2.org/
  labs(x="Predicted community dissimilarity",y = "Observed community dissimilarity")+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold",colour = "black", size = 7), 
    axis.text.y = element_text(face = "bold", size = 7, colour = "black"), 
    axis.title= element_text(face = "bold", size = 7, colour = "black"), 
    panel.background = element_blank(), 
    panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = c(0.8,0.9),
    legend.direction = "horizontal",
    legend.text = element_text(size = 7, face = "bold"),
    legend.title = element_text(size = 7, face = "bold"),
    aspect.ratio = 1
    
  )
gp2
library(ggpmisc)
gp2+stat_poly_eq(aes(label=paste(after_stat(eq.label),after_stat(adj.rr.label),after_stat(p.value.label),sep = "~~~~")),formula = y~x,parse=T,size=2.5)
ggsave("gdm_fitting.pdf",width = 8,height = 8,units = "cm")












# dis <- vegan::vegdist(t(otutab), method="bray")
# a <- as.matrix(dis)
# a <- as.data.frame(a)
# a$site <- colnames(a)
# env = read.delim("./gdm_env.txt")
# env[4:ncol(env)] <- data.frame(scale(env[,4:12]))
# exFormat3 <- formatsitepair(a, 3, XColumn="longitude", YColumn="latitude",
#                             predData=env, siteColumn="site")
# table(is.na(exFormat3 ))
# gdm.1 <- na.omit(exFormat3)
# table(is.na(gdm.1 ))
# gdm.1 <- gdm(data=gdm.1, geo=T)
# summary(gdm.1)
# 
# length(gdm.1$predictors) # get ideal of number of panels
# #> [1] 11
# plot(gdm.1, plot.layout=c(4,3))
# gdm.1.splineDat <- isplineExtract(gdm.1)
# str(gdm.1.splineDat)
# plot(gdm.1.splineDat$x[,"NO3"], 
#      gdm.1.splineDat$y[,"NO3"], 
#      lwd=3,
#      type="l", 
#      xlab="Winter precipitation (mm)", 
#      ylab="Partial ecological distance")
