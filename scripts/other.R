pvalue <- function(dt=NA, sample_map=NA, group=NA){
    dt = dt[,sample_map$ID]
    grps = unique(sample_map[,group])
    com = t(combn(grps,2))
    nspecies = nrow(dt)
    names = rownames(dt)
    result = rbind()
    for (n in 1:nspecies){
        temp_dt = dt[n,]
        for(c in 1:nrow(com)){
            g1 = com[c,1]
            g2 = com[c,2]
            g1s = sample_map[which(sample_map[,group] == g1),'ID']
            g2s = sample_map[which(sample_map[,group] == g2),'ID']
            dt1 = as.matrix(temp_dt[,g1s])
            dt2 = as.matrix(temp_dt[,g2s])
            c1 = sum(dt1 != 0 )
            c2 = sum(dt2 != 0)
            m1 = mean(dt1)
            m2 = mean(dt2)
            p = wilcox.test(dt1,dt2)$p.value
            temp_result = data.frame(name = names[n],g1=g1, g2=g2,mean1 = m1, mean2=m2,pvalue=p, count1=c1, count2= c2)
            result = rbind(result, temp_result)
        }
    }
    result
}
filter <- function(profile, id1, id2, q, fC){
    data.list <- profile %>% subset(g1 == id1 & g2 == id2)
    data <- data.list %>% subset(count1 >= 5 | count2 >= 5)
    data$qvalue <- p.adjust(data$pvalue,method = "BH")
    data$FC <- ifelse(data$mea$mean2>1, data$mean1/data$mean2,
                      data$mea$mean1)
    data$enrich <- ifelse(data$mean1 > data$mean2, data$g1, data$g2)
    data$log2FC <- ifelse(data$enrich == id1,log2(data$FC),-log2(data$FC))
    data$sig <- ifelse(data$qvalue < q, ifelse(data$FC > fC,"sig","non-sig"),"non-sig")
    data$enrich2 <- ifelse(data$sig == "sig",data$enrich,"Stable")
    data$log2FC <- ifelse(data$log2FC > 10, 10,
                          ifelse(data$log2FC < -10, -10, data$log2FC))
    marker <- subset(data, sig == "sig")
    marker <- marker[order(marker$log2FC),]
    marker$name <- factor(marker$name,levels = unique(marker$name))
    marker$label <- ifelse(marker$qvalue < 0.001,"***",
                           ifelse(marker$qvalue < 0.01,"**",
                                  ifelse(marker$qvalue < 0.05,"*","")))
    marker
}
