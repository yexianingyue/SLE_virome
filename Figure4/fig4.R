######random_forest#######

ntree_plot <- function(out=NA,p=NA){
    pdf(out)
    plot(p)
    dev.off()
}
importance_plot <- function(out=NA,p=NA){
    pdf(out)
    varImpPlot(p)
    dev.off()
}
roc_plot <- function(out=NA,p=NA, seed=NA){
    #pdf(out)
    plot.roc(p,add=F, 
             reuse.auc=TRUE,col="blue", partial.auc=c(1, 0.8), print.auc = F,
             print.auc.cex=2, print.auc.col='Black',
             max.auc.polygeon=T,
             axes=TRUE,
             grid = c(0.2,0.2),grid.col = 'grey')
    ciauc = round(ci.auc(p), digits = 3)
    textauc = paste("AUC: ", ciauc[2], " (", ciauc[1], " ~ ", ciauc[3],") seed: ",seed, sep="")
    text(x=0.4 ,y=0.4, labels=textauc, cex=2)
    #title(main="234", xlab="x", ylab='y')
    #dev.off()
}
zy_pcoa <- function(group=NA, level=NA,color=NA, sed=NA,dt=NA, sample_map=NA){

    set.seed(sed)

    #sample_map = sample_map[which(sample_map$group2 == group),]
    dt = dt[,as.character(sample_map$ID)]
    dtt = t(dt)
    dtt = dtt[,colSums(dtt !=0)>10] # 保留物种10个样本同时有的物种
    # split data
    index = sample(nrow(dtt),nrow(dtt)*0.2)
    train_data = dtt[index, ]; train_map = sample_map[index, ]
    test_data = dtt[-index, ]; test_map = sample_map[-index, ]

    fit = randomForest(train_data,as.factor(train_map[, color]), ntree=2000,importance=T, proximity=TRUE,)

    pred = as.data.frame(predict(fit, test_data, type='prob'))
    roc_p = roc(as.character(test_map[, color]), pred[,2])
    roc_plot("roc.test.pdf", roc_p, sed)
}

zy_pcoa(color = 'Group1',dt = t(dt.wms),sample_map = map.wms,sed = 2023)

##########LOO############
Loocv <- function(dt, group, k, map, seed){
    set.seed(seed)
    flag <- 1
    for(i in 1:k){  
        fold_test <- dt[i,]   #取1个样品作为测试集  
        fold_train <- dt[-i,]   # 剩下的样品作为训练集    

        # 注意：as.factor()中的Group指分组列的列名，非参数（如果不一样需要手动修改）。
        fit = randomForest(as.factor(Group2) ~ ., data = fold_train, 
                           ntree = 1000, importance = T, proximity = TRUE) 
        temp = as.data.frame(predict(fit, fold_test, type='prob'))
        if(flag == 1){
            pred = temp
            flag = 0
        }
        else {
            pred = rbind(pred,temp) 
        }
    }  
    roc1 <- roc(as.character(map[, group]),pred[,2])
    roc_plot("roc.test.pdf", roc1, seed = seed)
}

map.f <- read.table("../map.txt",sep="\t",header=T) 

dt.sig.loo <- merge(dt.wms.s, map.wms[,c(1,3)], by.x = "row.names",by.y = "ID")
rownames(dt.sig.loo) <- dt.sig.loo[,1]
dt.sig.loo <- dt.sig.loo[,-1]
colnames(dt.sig.loo) <- sub("-", "_", colnames(dt.sig.loo))

Loocv(dt = dt.sig.loo, group = "Group2", k = 47, map = map.wms, seed = 2022)

rf.result <- rbind()
for (i in 1:1000) {
    set.seed(i)
    xx = randomForest(as.factor(Group2)~., data=dt.sig.loo, ntree=2000, importance = TRUE, proximity = TRUE)
    x1 <- varImpPlot(xx)
    x2 <- data.frame(ID = rownames(x1),imp = x1[,1])
    if (i == 1 ) {
        rf.result <- x2  
    }else{
        rf.result <- merge(rf.result, x2,by.x = "ID",by.y = "ID")
    }
}

rownames(rf.result) <- rf.result$ID
rf.result <- rf.result[,-1]

rf.result2 <- rbind()
for (i in 1:length(rownames(rf.result))){
    x3 <- t.test(rf.result[i,])
    x4 <- data.frame(ID = rownames(rf.result)[i],estimate = x3$estimate, 
                     lw = x3$conf.int[1], hw = x3$conf.int[2])
    if (i == 1 ) {
        rf.result2 <- x4  
    }else{
        rf.result2 <- rbind(rf.result2, x4)
    }
}

write.table(rf.result2,"rf.res",quote = F,sep = "\t")

dt.new <- dt.wms.s
colnames(dt.new) <- sub("-", "_", colnames(dt.wms.s))
new.list <- rf.result2 %>% subset(estimate > 4)
dt.new <- dt.new[,as.character(new.list$ID)]

dt.sig.loo <- merge(dt.new, map.wms[,c(1,3)], by.x = "row.names",by.y = "ID")
rownames(dt.sig.loo) <- dt.sig.loo[,1]
dt.sig.loo <- dt.sig.loo[,-1]
colnames(dt.sig.loo) <- sub("-", "_", colnames(dt.sig.loo))

Loocv(dt = dt.sig.loo, group = "Group2", k = 47, map = map.wms, seed = 2022)
