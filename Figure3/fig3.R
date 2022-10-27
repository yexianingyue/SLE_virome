source("../scripts/zy_alpha_diversity.R")
source("../scripts/other.R")



sample_map = read.table("../data/sample.group", sep="\t", header=T)
sample_map = sample_map[sample_map$Group2 == "VLP"]

# fig 2a-c
#-------------------------
dt = read.table("../00.data/vOTU.profile.relative", quote="", sep="\t", header=T, row.names=1, check.names=F)
p1 <- zy_alpha(dt, sample_map, index="shannon")
p2 <- zy_alpha(dt, sample_map, index="simpson")
p0+p1+p2

# fig 2c
#-------------------------
dt = read.table("../data/vOTU.profile.relative", quote="",sep="\t", header=T, row.names=1, check.names=F)
p3 <- zy_pcoa(dt, sample_map=sample_map, zy_group="Group", ID="Sample", 
              sample.color=sample.color, title = "PCoA of SLE virome")
p3


# fig 2d
#-------------------------
taxo.color = c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#4e5c68","#72479f","#c5835e","#4d545c")
source("../scripts/zy_compositions.R")
p4 <- zy_group_compositions(dt, sample_map, ID="Sample", group="Group", top_N=12, 
                            taxo.color = taxo.color, order_func = "cluster", title = "Relative abundance at Genus level")+
theme(axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_blank())
p4


# fig 2e
#-------------------------
family.p <-  as.data.frame(t(read.table("../data/vOTU.family.s.profile", sep="\t", 
                                        header=T, row.names = 1,check.names = F)))
family.p <- family.p[sample_map$Sample,]
family.p <- family.p/rowSums(family.p) 
family.s <- pvalue(as.data.frame(t(family.p)), sample_map, group = "Group")
family.w <- filter(family.s,"HC","SLE",0.05,0)
family.w$label <- ifelse(family.w$pvalue < 0.001, "***",
                         ifelse(family.w$pvalue < 0.01,"**",
                                ifelse(family.w$pvalue < 0.05,"*","")))
family.w <- family.w[order(family.w$log2FC),] 
family.w$name <- factor(family.w$name,levels = unique(family.w$name))

p5 <- ggplot(family.w, aes(x = name, y = log2FC, fill = enrich))+
geom_bar(stat = "identity", position = "identity")+
geom_text(data = family.w, aes(x = name, y = log2FC, label = label))+
#scale_fill_manual(values = c("purple", "blue"), guide=FALSE)+
xlab("")+ylab("")+
coord_flip()+
theme_bw()



