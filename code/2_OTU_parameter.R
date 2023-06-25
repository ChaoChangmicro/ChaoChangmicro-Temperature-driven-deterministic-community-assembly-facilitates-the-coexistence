library(tidyverse)
library(phyloseq)
library(ggClusterNet)
#--输入和设定文件#-----
ps0 = base::readRDS("ps.rds")

# # #---解决错误
# ps0 = base::readRDS("./ps.rds")
# ps0


# tax = ps0 %>% vegan_tax() %>%
#   as.data.frame()
# head(tax)
# colnames(tax) =  c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
# tax_table(ps) = as.matrix(tax)

# 基于phyloseq对象进行分析
# 设定主题，颜色，分组数量，分组排序，比较类型等内容
# 如果是微生物数据，可以定义进化树展示的微生物数量，默认150个
# 同时设定lefse用的微生物的数量


#-主题--颜色等-#-------
source("F:/2  生信代码/4 三人成师/micro/total_amplicon.R")
#---扩增子环境布置
res = theme_my(ps0)
mytheme1 = res[[1]]
mytheme2 = res[[2]]; 
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]
#--如果出现错误，设定smart = F；是由于注释信息不符合本代码标准导致的
result<- dir.amp(ps0 = ps0,smart = F)
res1path = result[[1]];res1path
id = result[[2]];id

# 微生物组数据输入#------

#--如何筛选样本
# ps_sub <- subset_samples(ps0,!ID %in% c("sample1"));ps_sub 
# 需要挑选部分微生物进行分析

if (F) {
  ps0 <- ps0 %>%
    subset_taxa(
      Kingdom %in% id
    ) ;ps0
  
}

#--最终确定的phyloseq对象定义为ps#--------
ps = ps0
ps = ps %>% filter_taxa(function(x) sum(x ) > 0 , TRUE)
#--提取有多少个分组#-----------
gnum = phyloseq::sample_data(ps)$Group %>% unique() %>% length()
gnum
# 设定排序顺序
# axis_order =  phyloseq::sample_data(ps)$Group %>%unique(fromLast = T);axis_order
axis_order <-c("Low","Med","High"  )
#--默认进化树展示的微生物数量
Top_micro = 150

# 堆叠柱状图展示Top前是个门,j 展示的水平为门水平#-----
Top = 12
# rank.names(ps0)
jj = j = "Phylum"

# -类韦恩图的二分网络设置过滤阈值#-------
ps_biost = ggClusterNet::filter_OTU_ps(ps = ps,Top = 500)

# 差异分析设定两两比对分组#------
# group1 = c("Gro1","Gro2")
# group2 = c("Gro1","Gro2")
# b= data.frame(group1,group2)
b = NULL# 如果每两个组之间都做差异，那就指定


# 热图等差异系数最大的N个OTU#--------
heatnum　=　30


#--R 语言做lefse法分析-过滤#----------
ps_Rlefse = ggClusterNet::filter_OTU_ps(ps = ps,Top = 400)

#---机器学习部分#------
# ROC是三种机器学习的ROC曲线，但是只能跑两个分组，如果两组，可以选择为T。
ROC = FALSE
# 是否做交叉检验
rfcv = FALSE
# 选择多少个重要变量
optimal = 40

#-------功能预测#----

if (is.null(ps0@refseq)) {
  Tax4Fun = FALSE
} else if(!is.null(ps0@refseq)){
  Tax4Fun = TRUE
}

ps.t = ps0 %>% ggClusterNet::filter_OTU_ps(500)
if (Tax4Fun) {
  dir.create("data")
  otu = ps.t %>% 
    # ggClusterNet::filter_OTU_ps(1000) %>%
    ggClusterNet:: vegan_otu() %>%
    t() %>%
    as.data.frame()
  # write.table(otu,"./data/otu.txt",quote = F,row.names = T,col.names = T,sep = "\t")
  rep = ps.t %>% 
    # ggClusterNet::filter_OTU_ps(1000) %>%
    phyloseq::refseq()
  rep
  # library(Biostrings)
  Biostrings::writeXStringSet(rep,"./data/otu.fa")
  ps.t = ps.t
  
  #开头空一格字符保存
  write.table("\t", "./data/otu.txt",append = F, quote = F, eol = "", row.names = F, col.names = F)
  # 保存统计结果，有waring正常
  write.table(otu, "./data/otu.txt", append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)     

}

#-----选择性功能#---------

#设置CK，用于双向柱状图绘制-目前不绘制#------
CK = unique(phyloseq::sample_data(ps)$Group)[1]

# 用python做lefse#-----
lefse.py = F
if (lefse.py) {
  lefsenum = 0
  ps_lefse <- ps %>%
    phyloseq::subset_taxa(
      # Kingdom == "Fungi"
      Kingdom == id
      # Genus  == "Genus1"
      # Species %in%c("species1") 
      # row.names(tax_table(ps0))%in%c("OTU1")
    )
  
  ps_lefse = ggClusterNet::filter_OTU_ps(ps = ps_lefse,Top = 400)
  
  # #文件预处理
  # format_input.py LEFSE_to_run_G_level.txt pri_lefse.in -c 1 -u 2 -o 1000000
  # # 注意这里 –c用来指定分组信息-u 1指定样品信息
  # #文件分析,这里-l设置LDA阈值，默认为2，我们使用4 会更加严格
  # ~/src/nsegata-lefse/run_lefse.py pri_lefse.in pri_lefse_2.res  -l 2
  # #柱状图绘制
  # plot_res.py pri_lefse_2.res lefse_barplot.pdf --format pdf
  # #树状图绘制
  # plot_cladogram.py pri_lefse_2.res lefse_tree.pdf --format pdf
  # #做每个差异的柱状图
  # mkdir biomarkers_raw_images
  # plot_features.py pri_lefse.in pri_lefse_2.res biomarkers_raw_images/
  
}






