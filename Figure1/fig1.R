source("../scripts/zy_alpha_diversity.R")
source("../scripts/other.R")

sample_map = read.table("sample.group", sep="\t", header=T)
sample_map = sample_map[sample_map$Group2 == "WMS"]

# fig 1a-c
#-------------------------
dt = read.table("vOTU.profile.relative", quote="", sep="\t", header=T, row.names=1, check.names=F)
p0 <- zy_nspecies(dt, sample_map)
p1 <- zy_alpha(dt, sample_map, index="shannon")
p2 <- zy_alpha(dt, sample_map, index="simpson")
p0+p1+p2

# fig 1d
#-------------------------
dt = read.table("vOTU.profile.relative", quote="",sep="\t", header=T, row.names=1, check.names=F)
p3 <- zy_pcoa(dt, sample_map=sample_map, zy_group="Group", ID="Sample", 
              sample.color=sample.color, title = "PCoA of SLE virome")
p3


# fig 1e
#-------------------------
taxo.color = c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#4e5c68","#72479f","#c5835e","#4d545c")
source("../scripts/zy_compositions.R")
p4 <- zy_group_compositions(dt, sample_map, ID="Sample", group="Group", top_N=12, 
                            taxo.color = taxo.color, order_func = "cluster", title = "Relative abundance at Genus level")+
theme(axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_blank())
p4


# fig 1f
#-------------------------
family.p <-  as.data.frame(t(read.table("vOTU.family.s.profile", sep="\t", 
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





# fig 1g
#-------------------------
dt = read.table("marker.votu", sep="\t",header=T)

p6 <- ggplot(dt, aes(x=log2FC, y=-log10(qvalue), color=enrich2))+
geom_point()+
theme_bw()+
geom_hline(yintercept=-log10(0.05), colour="red", linetype="dashed")+
geom_hline(yintercept=-log10(0.01), colour="red", linetype="dashed")+
geom_vline(xintercept=1, colour="red", linetype="dashed")+
geom_vline(xintercept=-1, colour="red", linetype="dashed")+
scale_color_manual(values=c("#ee7671","#27b9b9","gray"))

# fig 1h
#-------------------------
host <- read.table("sig_votu.host",sep="\t", header=T)
host <- host[order(-host$Count),]
host$Host <- factor(host$Host,levels = unique(host$Host))
p7 <- ggplot(host, aes(Family, Count, fill = Host)) +
geom_bar(position = 'stack', width = 1,stat="identity") +
labs(x = 'Family', y = 'sig vOTU count') + 
scale_fill_manual(values=c("#8dd3c7","#ffffb3","#bebada","#fb8072",
                           "#80b1d3","#fdb462","#b3de69","#fccde5",
                           "#d9d9d9","#bc80bd","#ccebc5","#4e5c68",
                           "#72479f","#c5835e","#4d545c"

                           ))+
facet_wrap(Enriched~.,scales = 'free_x')+
theme_bw()


