##计算微生物丰度间的相关系数
library(Hmisc)

#以属水平丰度为例，“genus_table.txt” 是一个属水平的微生物丰度表
genus <- read.delim('all.txt', row.name = 1, check.names = FALSE)
genus <- t(genus)
genus <- genus /rowSums(genus)
genus <- t(genus)
#可选事先过滤一些低丰度或低频的类群
genus <- genus[which(rowSums(genus) >= 0.0005), ]    #例如只保留相对丰度总和高于 0.005 的属

genus1 <- genus
genus1[genus1>0] <- 1
 genus <- genus[which(rowSums(genus1) >= 2), ]    #例如只保留在 5 个及以上样本中出现的属

#计算两属之间是否存在丰度变化的相关性，以 spearman 相关系数为例
genus_corr <- rcorr(t(genus), type = 'spearman')

#阈值筛选
#将 spearman 相关系数低于 0.7 的关系剔除，即 r>=0.7
r <- genus_corr$r
r[abs(r) < 0.6] <- 0

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
p <- genus_corr$P
p <- p.adjust(p, method = 'fdr')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0

#根据上述筛选的 r 值和 p 值保留数据
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]

#如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
write.table(data.frame(z, check.names = FALSE), 'genus_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

##获得网络
library(igraph)

#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物属间丰度的 spearman 相关系数
g <- graph_from_adjacency_matrix(z, weighted = TRUE, mode = 'undirected')
g

#自相关也可以通过该式去除
g <- simplify(g)

#孤立节点的删除（删除度为 0 的节点）
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))

#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)

igraph <- g
##节点特征
#节点数量
length(V(igraph)$name)
#或
vcount(igraph)

#节点度（Degree）
#由于本示例是个无向网络，故无出度和入度之分
V(igraph)$degree <- degree(igraph)
V(igraph)$degree

#查看度分布
#可观察到微生物相关网络通常服从幂律分布，这个下节再讲怎样通过计算验证
degree_dist <- degree.distribution(igraph)[-1]
degree_num <- 1:max(V(igraph)$degree)

par(mfrow = c(1, 2))
hist(V(igraph)$degree, xlab = 'Degree', ylab = 'Frequency',
     main = 'Degree distribution')
plot(degree_num, degree_dist, log = 'xy', xlab = 'Log-degree',
     ylab = 'Log-intensity', main = 'Log-log degree distribution')

#查看节点度与其“邻居”的平均度的关系
#微生物网络中高度值的节点更倾向连接在一起，是普遍现象吗？
neighbor_degree <- graph.knn(igraph, V(igraph))$knn
plot(V(igraph)$degree, neighbor_degree, log = 'xy',
     xlab = 'Log degree', ylab = 'Log average neighbor degree')

#加权度（Weighted degree）
V(igraph)$weight_degree <- strength(igraph)
V(igraph)$weight_degree

#接近中心性（Closeness centrality）
V(igraph)$closeness_centrality <- closeness(igraph)
V(igraph)$closeness_centrality

#介数中心性（Betweenness centrality）
V(igraph)$betweenness_centrality <- betweenness(igraph)
V(igraph)$betweenness_centrality

#特征向量中心性（Eigenvector centrality）
V(igraph)$eigenvector_centrality <- evcent(igraph)$vector
V(igraph)$eigenvector_centrality

#探索三种描述节点中心性的特征的关系
library(car)

scatter3d(V(igraph)$closeness_centrality, V(igraph)$betweenness_centrality, V(igraph)$eigenvector_centrality,
          xlab =  'Closeness centrality', ylab = 'Betweenness centrality', zlab = 'Eigenvector centrality',
          surface = FALSE)

#探索节点度和节点中心性的关系，如与特征向量中心性的关系
plot(V(igraph)$degree, V(igraph)$eigenvector_centrality,
     xlab = 'Degree', ylab = 'Eigenvector centrality')





#为节点（微生物属）添加属性信息（界门纲目科属水平注释）
#“genus_taxonomy.txt” 记录了微生物的属性，读入该表后根据已知网络节点匹配对应的行
tax <- read.delim('taxon.txt', row.name = 1, check.names = FALSE, stringsAsFactors = FALSE)
tax <- tax[as.character(V(g)$name), ]

V(igraph)$kingdom <- tax$Kingdom
V(igraph)$phylum <- tax$Phylum
V(igraph)$class <- tax$Class
V(igraph)$order <- tax$Order
V(igraph)$family <- tax$Family
V(igraph)$genus <- tax$Genus

#查看网络图
igraph
plot(igraph)

##网络文件输出，输出特定的网络文件类型，便于后续数据分析需求
#邻接矩阵，出了上述提到的在计算相关系数后，输出筛选后的相关系数矩阵外
#还可以由 igraph 的邻接列表转换



#输出节点列表
node_list <- data.frame(
  node_id = V(igraph)$name,
  degree = V(igraph)$degree,
  weight_degree = V(igraph)$weight_degree,
  closeness_centrality = V(igraph)$closeness_centrality,
  betweenness_centrality = V(igraph)$betweenness_centrality,
  eigenvector_centrality = V(igraph)$eigenvector_centrality,
  kingdom = V(igraph)$kingdom,
  phylum = V(igraph)$phylum,
  class = V(igraph)$class,
  order = V(igraph)$order,
  family = V(igraph)$family,
  genus = V(igraph)$genus )

head(node_list)
write.table(node_list, 'node_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)


##边特征
#边的数量
ecount(igraph)

#权重（Weighted），已在数据读入时转化获得
E(igraph)$weight

#边介数中心性（Edge betweenness centrality）
E(igraph)$betweenness_centrality <- edge.betweenness(igraph)
E(igraph)$betweenness_centrality

#输出列表
edge <- data.frame(as_edgelist(igraph))    #igraph 的邻接列表转为边列表

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(igraph)$weight,
  correlation = E(igraph)$correlation,
  betweenness_centrality = E(igraph)$betweenness_centrality
)
head(edge_list)

write.table(edge_list, 'edge_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.graph(igraph, 'network.graphml', format = 'graphml')
