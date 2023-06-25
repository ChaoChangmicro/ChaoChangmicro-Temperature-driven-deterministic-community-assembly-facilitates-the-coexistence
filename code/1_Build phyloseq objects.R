library(phyloseq)
library(tidyverse)
metadata = read.delim("./metadata.txt",row.names = 1)
metadata$ID = row.names(metadata)

#qiime2下机数据抽平后
otutab = read.delim("./asv100.txt", row.names=1)
head(otutab)

taxonomy = read.table("./taxonomy.txt", row.names=1,header = T)
head(taxonomy)

library(Biostrings)
rep = readDNAStringSet("./dna-sequences.fasta")
seq <- rep[names(rep) %in% row.names(otutab)]

ps = phyloseq(
  sample_data(metadata),
  otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
  tax_table(as.matrix(taxonomy)),
  refseq(rep)
)


#将otu的名字进行修改
taxa_names(ps) <- 1:ntaxa(ps)
taxa_names(ps) <- paste("ASV", taxa_names(ps), sep="")
taxa_names(ps)

a <- as.data.frame(otu_table(ps))
write.csv(a ,file = "asv_rarefied.csv")

b <- as.data.frame(tax_table(ps))
write.csv(b ,file = "tax_rarefied.csv")

c <- refseq(ps)
writeXStringSet(c , "./asv.fasta", format="fasta")



######构建进化树

# 起始文件为 result/tree目录中 otus.fa(序列)、annotation.txt(物种和相对丰度)文件
# Muscle软件进行序列对齐，3s
#cd到研究目录下
/mnt/c/bin/muscle.exe -in asv.fasta -out asv_aligned.fasta

### 方法1. 利用IQ-TREE快速构建ML进化树，2m
# rm -rf iqtree
# mkdir -p iqtree
# iqtree -s otus_aligned.fas \
# -bb 1000 -redo -alrt 1000 -nt AUTO \
# -pre iqtree/otus

## 方法2. FastTree快速建树(Linux)
# 注意FastTree软件输入文件为fasta格式的文件，而不是通常用的Phylip格式。输出文件是Newick格式。
# 该方法适合于大数据，例如几百个OTUs的系统发育树！
# Ubuntu上安装fasttree可以使用`apt install fasttree`
fasttree -gtr -nt asv_aligned.fasta > asv.nwk

####将原来的树裁剪为需要的树####
# library(picante)
# tree <- read.tree("rooted_tree.nwk")
# otu <- read.table("asv_rare.txt",sep = "\t",row.names = 1,header = T)
# otu <- t(otu)
# #裁剪时otu的名字必须为列名
# phy.tree <- prune.sample(otu,tree)
# plot(phy.tree)
# write.tree(phy.tree,file = "pruned_tree.nwk")


# library(seqinr)
# all_fasta <- read.fasta("dna-sequences.fasta")
# #选取
# sub_fasta <- all_fasta[names(all_fasta) %in% colnames(otu)]
# # 写出文件
# write.fasta(sequences = sub_fasta, names =names(sub_fasta), file.out = 'target.fasta')

# 构建phyloseq对象 ------------------------------------------------------------
library(phyloseq)
library(tidyverse)
metadata = read.delim("./metadata.txt",row.names = 1)
metadata$ID = row.names(metadata)


otutab = read.delim("./asv_rarefied.txt", row.names=1)
head(otutab)

taxonomy = read.table("./tax_rarefied.txt", row.names=1,header = T)
head(taxonomy)

tree  = read_tree("./asv.nwk")


library(Biostrings)
rep = readDNAStringSet("./asv.fasta")

ps = phyloseq(
  sample_data(metadata),
  otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
  tax_table(as.matrix(taxonomy)), 
  phy_tree(tree),
  refseq(rep)
)

saveRDS(ps,"ps.rds")

