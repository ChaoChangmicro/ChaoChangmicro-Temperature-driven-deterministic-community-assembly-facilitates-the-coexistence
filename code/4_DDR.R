library(vegan)

#读取上述数据集
env<- read.csv('多样性.csv',header = T,row.names = 1,sep = ',')
abund<-read.csv('all.csv',header = T,row.names = 1,sep = ',')
abund <- t(abund)
dist.abund <- vegdist(abund, method = 'bray')

#根据环境测量指标，计算样方间的欧几里得距离
#这里只选择了其中的温度指标，期望关注物种变化与温度的相关性
# temp <- env$tem
# dist.temp <- dist(temp, method = 'euclidean')



#根据经纬度，计算样方间实际的地理距离
library(geosphere)
geo <- data.frame(env$longitude, env$latitude)
d.geo <- distm(geo, fun = distHaversine)      
dist.geo <- as.dist(d.geo)

##执行 Mantel tests，详情 ?mantel，以下为 3 个示例
#物种丰度和温度的相关性，以 spearman 相关系数为例，9999 次置换检验显著性（Mantel tests 基于随机置换的方法获取 p 值）
# abund_temp <- mantel(dist.abund, dist.temp, method = 'spearman', permutations = 9999, na.rm = TRUE)
# abund_temp

#物种丰度和地理距离的相关性，以 spearman 相关系数为例，9999 次置换检验显著性
# abund_geo <- mantel(dist.abund, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
# abund_geo

library(ggplot2)
#将上文获得的距离测度，转化为数据框，一一对应起来
aa <- as.vector(dist.abund)
# tt <- as.vector(dist.temp)
gg <- as.vector(dist.geo)
mat <- data.frame(aa,  gg)


#基于物种丰度的距离与基于温度指标的距离之间的相关性散点图，上文已知二者显著相关；同时颜色表示样方间地理距离
mm <- ggplot(mat, aes(y = aa* 100, x = gg/100000)) + 
  geom_point(size =2, alpha = 0.75, colour = "black", fill = "#8c88c0",shape=21) + 
  geom_smooth(method = "lm", colour = "black", alpha = 0.2) + 
  labs(x = "Geographic Distance (100Km)", y = "Bray-Curtis Dissimilarity (%)") + 
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
    axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
    axis.title= element_text(face = "bold", size = 14, colour = "black"), 
    panel.background = element_blank(), 
    panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = c(0.8,0.9),
    legend.direction = "horizontal",
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    aspect.ratio = 1
    
  )
library(ggpmisc)
mm+stat_poly_eq(aes(label=paste(after_stat(eq.label),after_stat(adj.rr.label),after_stat(p.value.label),sep = "~~~~")),formula = y~x,parse=T,size=2.5)

x <- lm(gg~aa,mat)
summary(x)
