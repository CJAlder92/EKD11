### Cluster Expression
## This purpose of this script is to identify to clusters within a dataset of already identified DEGs and then visualise the data

rm(list = ls())
source('/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/1-R/refactor/packages.R')

nf.dir      <- "/Users/alderc/1-projects/CAMP/1-AS_timecourse/2-EKD11/1-Pipeline/1-Nextflow/"
work.dir    <- "/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/"
r.dir       <- paste(work.dir, "1-R/",sep='') #change once refactor is completed
tmp.dir     <- paste(work.dir,"tmp/",sep='')
data.dir    <- paste(work.dir,"4-data/",sep='')
results.dir <- paste(work.dir,"2-results_2020/",sep='')
genome.dir  <- "/Users/alderc/1-projects/9-Data/1-Reference_genomes/1-Mus_musculus/"

# Load data
load(paste(r.dir, "Projects/deseq.Rdata", sep = ""))
load(paste(r.dir, "Projects/results.Rdata", sep = ""))

#Extra packages
library(dynamicTreeCut)

deg.genes <- unique(unlist(sapply(res.list, function(df){
  genes <- rownames(df[df$DEG, ])
  genes
})))

deg.genes

dat <- as.matrix(assay(vsd)[deg.genes, ])



## MT and SBP

 
mt.sbp.scaled <- t(scale(t(mt.sbp.dat)))
mt.sbp.cor <- as.dist(1-cor(t(mt.sbp.scaled), method="pearson"))
hr <- hclust(mt.sbp.cor, method="complete") # Cluster Genes


clusDyn <- cutreeDynamic(hr, distM = as.matrix(mt.sbp.cor), method = "hybrid")
names(clusDyn) <- rownames(mt.sbp.scaled)

clus.centroids <- sapply(levels(factor(clusDyn)), clust.centroid, mt.sbp.scaled, clusDyn)

centroid.melt <- melt(clus.centroids)
colnames(centroid.melt) <- c("Sample", "Cluster", "value")

par(mfrow = c(7,3))

for (i in 1:length(unique(centroid.melt$Cluster))){
  centroid.df <- subset(centroid.melt, Cluster == i)
  centroid.df$day <- str_match(centroid.df$Sample, pattern= "\\.d(.+)\\.")[,2]
  centroid.df$delivery <- str_match(centroid.df$Sample, pattern="^(.*)\\..*\\.")[,2]
  # centroid.df <- centroid.df[5:nrow(centroid.df), ]
  
  p1 <- ggplot(centroid.df, aes(x = day, y = value, group= delivery, colour=as.factor(delivery))) +
    geom_point() +
    # stat_summary(aes(y = value,group=1), fun.y=mean, colour="red", geom="line",group=1) + 
    geom_smooth(se = FALSE, alpha=0.1, method = "loess", span=0.9) +
    xlab("Time") +
    ylab("Expression") +
    labs(title= paste("Cluster", i, "Expression by Time", sep = " "),color = "Cluster")
  
  plot(p1)
}


# MT vs SBP 

names(res.list)

mt.deg <- lapply(res.list[grep("^mt.+naive$", names(res.list))], function(df){
  df <- df[df$DEG, ]
})

sbp.deg <- lapply(res.list[grep("^sbp.+naive$", names(res.list))], function(df){
  df <- df[df$DEG, ]
})

mt.deg.list <- unique(unlist(sapply(mt.deg, row.names)))
sbp.deg.list <- unique(unlist(sapply(sbp.deg, row.names)))

common.list <- intersect(mt.deg.list, sbp.deg.list)

common.dat <- assay(vsd[common.list, grep("naive|^mt|sbp", colnames(vsd))])

common.scaled <- t(scale(t(common.dat)))

common.scaled.cor <- as.dist(1-cor(t(common.scaled), method = "pearson"))

hr <- hclust(common.scaled.cor, method="complete") # Cluster Genes

clusDyn <- cutreeDynamic(hr, distM = as.matrix(common.scaled.cor), method = "hybrid")
names(clusDyn) <- rownames(common.scaled)

clus.centroids <- sapply(levels(factor(clusDyn)), clust.centroid, common.scaled, clusDyn)

centroid.melt <- melt(clus.centroids)
colnames(centroid.melt) <- c("Sample", "Cluster", "value")

for (i in 1:length(unique(centroid.melt$Cluster))){
  centroid.df <- subset(centroid.melt, Cluster == i)
  centroid.df$day <- ifelse(grepl("naive", centroid.df$Sample), 0,
                            str_match(centroid.df$Sample, pattern= "\\.d(.+)\\.")[,2])
  centroid.df$delivery <- ifelse(grepl("naive", centroid.df$Sample), "Naive",
                                 str_match(centroid.df$Sample, pattern="^(.*)\\..*\\.")[,2])
  # centroid.df <- centroid.df[5:nrow(centroid.df), ]
  
  p1 <- ggplot(centroid.df, aes(x = day, y = value, group= delivery, colour=as.factor(delivery))) +
    geom_point() +
    # stat_summary(aes(y = value,group=1), fun.y=mean, colour="red", geom="line",group=1) + 
    geom_smooth(se = FALSE, alpha=0.1, method = "loess", span=0.9) +
    xlab("Time") +
    ylab("Expression") +
    labs(title= paste("Cluster", i, "Expression by Time", sep = " "),color = "Cluster")
  
  pdf(paste(results.dir, "clusters/mt_sbp2/cluster_", i, ".pdf", sep =""))
  plot(p1)
  dev.off()
}

clusDyn["ENSMUSG00000055170"]
## Cluster 1 
ind = (clusDyn == 15)
clus1 <- names(clusDyn[ind])

heat.dat <-as.matrix(assay(vsd)[clus1, ])
samples.rmv <- grep("rtmt", colnames(heat.dat))
heat.dat <- heat.dat[, -samples.rmv]
heat.dat   <- t(apply(heat.dat,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))}))
rownames(heat.dat) <- gtf.dat[clus1,"gene_name"]

#order so days are together 
heat.list <- colnames(heat.dat);
heat.list <- as.data.frame(heat.list);
heat.list$day <- ifelse(grepl("naive", heat.list$heat.list),
                        0,
                        str_match(string=heat.list$heat.list, pattern= "\\.d(.+)\\.")[,2])
heat.list$delivery <- ifelse(grepl("naive", heat.list$heat.list),
                             "naive",
                             str_match(string=heat.list$heat.list, pattern="^(.*)\\..*\\.")[,2])
heat.list$delivery <- ordered(factor(heat.list$delivery, levels= c('naive', 'sbp', 'rtmt', 'mt')))
heat.list$rep <- as.integer(str_match(string=heat.list$heat.list, pattern=".*r(.)$")[,2])
heat.list <- heat.list[order(as.numeric(heat.list$day), as.numeric(as.factor(heat.list$delivery)), heat.list$rep),]
heat.list <- as.vector(heat.list$heat.list);
heat.dat <- heat.dat[,heat.list] 

pdf(paste(results.dir, "clusters/mt_sbp2/cluster_15_hm.pdf", sep = "" ), height = 11, width = 8)
do_heatmap(heat.dat, title = "Cluster 15")
dev.off()
grep("Ifng", rownames(heat.dat))




## All conditions 

names(res.list)

mt.deg <- lapply(res.list[grep("^mt.+naive$", names(res.list))], function(df){
  df <- df[df$DEG, ]
})

sbp.deg <- lapply(res.list[grep("^sbp.+naive$", names(res.list))], function(df){
  df <- df[df$DEG, ]
})

rmt.deg <- lapply(res.list[grep("^rtmt.+naive$", names(res.list))], function(df){
  df <- df[df$DEG, ]
})

mt.deg.list <- unique(unlist(sapply(mt.deg, row.names)))
sbp.deg.list <- unique(unlist(sapply(sbp.deg, row.names)))
rmt.deg.list <- unique(unlist(sapply(rmt.deg, row.names)))

common.list <- Reduce(intersect, list(mt.deg.list, sbp.deg.list, rmt.deg.list))

common.dat <- assay(vsd[common.list,])

common.scaled <- t(scale(t(common.dat)))

common.scaled.cor <- as.dist(1-cor(t(common.scaled), method = "pearson"))

hr <- hclust(common.scaled.cor, method="complete") # Cluster Genes

clusDyn <- cutreeDynamic(hr, distM = as.matrix(common.scaled.cor), method = "hybrid")
names(clusDyn) <- rownames(common.scaled)

clus.centroids <- sapply(levels(factor(clusDyn)), clust.centroid, common.scaled, clusDyn)

centroid.melt <- melt(clus.centroids)
colnames(centroid.melt) <- c("Sample", "Cluster", "value")

for (i in 1:length(unique(centroid.melt$Cluster))){
  centroid.df <- subset(centroid.melt, Cluster == i)
  centroid.df$day <- ifelse(grepl("naive", centroid.df$Sample), 0,
                            str_match(centroid.df$Sample, pattern= "\\.d(.+)\\.")[,2])
  centroid.df$delivery <- ifelse(grepl("naive", centroid.df$Sample), "naive",
                                 str_match(centroid.df$Sample, pattern="^(.*)\\..*\\.")[,2])
  # centroid.df <- centroid.df[5:nrow(centroid.df), ]
  cols <- c("mt" = "grey50", "sbp" = "deeppink2", "naive" = "black", "rtmt" = "dodgerblue2")
  
  p1 <- ggplot(centroid.df, aes(x = day, y = value, group= delivery, colour=as.factor(delivery))) +
    geom_point() +
    # stat_summary(aes(y = value,group=1), fun.y=mean, colour="red", geom="line",group=1) + 
    geom_smooth(se = FALSE, alpha=0.1, method = "loess", span=0.9) +
    xlab("Time") +
    ylab("Expression") +
    labs(title= paste("Cluster", i, "Expression by Time", sep = " "),color = "Cluster") +
    scale_colour_manual(values = cols)
  
  pdf(paste(results.dir, "clusters/common_all/common_cluster_", i, ".pdf", sep =""))
  plot(p1)
  dev.off()
}

ifng.clus <- clusDyn["ENSMUSG00000055170"] # IFNG

ind = (clusDyn == ifng.clus)
clus1 <- names(clusDyn[ind])

heat.dat <-as.matrix(assay(vsd)[clus1, ])
heat.dat   <- t(apply(heat.dat,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))}))
rownames(heat.dat) <- gtf.dat[clus1,"gene_name"]

#order so days are together 
heat.list <- colnames(heat.dat);
heat.list <- as.data.frame(heat.list);
heat.list$day <- ifelse(grepl("naive", heat.list$heat.list),
                        0,
                        str_match(string=heat.list$heat.list, pattern= "\\.d(.+)\\.")[,2])
heat.list$delivery <- ifelse(grepl("naive", heat.list$heat.list),
                             "naive",
                             str_match(string=heat.list$heat.list, pattern="^(.*)\\..*\\.")[,2])
heat.list$delivery <- ordered(factor(heat.list$delivery, levels= c('naive', 'sbp', 'rtmt', 'mt')))
heat.list$rep <- as.integer(str_match(string=heat.list$heat.list, pattern=".*r(.)$")[,2])
heat.list <- heat.list[order(as.numeric(heat.list$day), as.numeric(as.factor(heat.list$delivery)), heat.list$rep),]
heat.list <- as.vector(heat.list$heat.list);
heat.dat <- heat.dat[,heat.list] 

pdf(paste(results.dir, "clusters/common_all/cluster_2_hm.pdf", sep = "" ), height = 20, width = 10)
do_heatmap(heat.dat, title = "Cluster 2")
dev.off()
grep("Ifng", rownames(heat.dat))


mt.rmt <- intersect(mt.deg.list, rmt.deg.list)
mt.rmt <- setdiff(mt.rmt, common.list)

common.dat <- assay(vsd[mt.rmt,])

common.scaled <- t(scale(t(common.dat)))

common.scaled.cor <- as.dist(1-cor(t(common.scaled), method = "pearson"))

hr <- hclust(common.scaled.cor, method="complete") # Cluster Genes

clusDyn <- cutreeDynamic(hr, distM = as.matrix(common.scaled.cor), method = "hybrid")
names(clusDyn) <- rownames(common.scaled)

clus.centroids <- sapply(levels(factor(clusDyn)), clust.centroid, common.scaled, clusDyn)

centroid.melt <- melt(clus.centroids)
colnames(centroid.melt) <- c("Sample", "Cluster", "value")

for (i in 1:length(unique(centroid.melt$Cluster))){
  centroid.df <- subset(centroid.melt, Cluster == i)
  centroid.df$day <- ifelse(grepl("naive", centroid.df$Sample), 0,
                            str_match(centroid.df$Sample, pattern= "\\.d(.+)\\.")[,2])
  centroid.df$delivery <- ifelse(grepl("naive", centroid.df$Sample), "naive",
                                 str_match(centroid.df$Sample, pattern="^(.*)\\..*\\.")[,2])
  # centroid.df <- centroid.df[5:nrow(centroid.df), ]
  cols <- c("mt" = "grey50", "sbp" = "deeppink2", "naive" = "black", "rtmt" = "dodgerblue2")
  
  p1 <- ggplot(centroid.df, aes(x = day, y = value, group= delivery, colour=as.factor(delivery))) +
    geom_point() +
    # stat_summary(aes(y = value,group=1), fun.y=mean, colour="red", geom="line",group=1) + 
    geom_smooth(se = FALSE, alpha=0.1, method = "loess", span=0.9) +
    xlab("Time") +
    ylab("Expression") +
    labs(title= paste("Cluster", i, "Expression by Time", sep = " "),color = "Cluster") +
    scale_colour_manual(values = cols)
  
  pdf(paste(results.dir, "clusters/mt_rmt_shared/common_cluster_", i, ".pdf", sep =""))
  plot(p1)
  dev.off()
}

heat.dat <-as.matrix(assay(vsd)[mt.rmt, ])
heat.dat   <- t(apply(heat.dat,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))}))
rownames(heat.dat) <- gtf.dat[mt.rmt,"gene_name"]

#order so days are together 
heat.list <- colnames(heat.dat);
heat.list <- as.data.frame(heat.list);
heat.list$day <- ifelse(grepl("naive", heat.list$heat.list),
                        0,
                        str_match(string=heat.list$heat.list, pattern= "\\.d(.+)\\.")[,2])
heat.list$delivery <- ifelse(grepl("naive", heat.list$heat.list),
                             "naive",
                             str_match(string=heat.list$heat.list, pattern="^(.*)\\..*\\.")[,2])
heat.list$delivery <- ordered(factor(heat.list$delivery, levels= c('naive', 'sbp', 'rtmt', 'mt')))
heat.list$rep <- as.integer(str_match(string=heat.list$heat.list, pattern=".*r(.)$")[,2])
heat.list <- heat.list[order(as.numeric(heat.list$day), as.numeric(as.factor(heat.list$delivery)), heat.list$rep),]
heat.list <- as.vector(heat.list$heat.list);
heat.dat <- heat.dat[,heat.list] 

pdf(paste(results.dir, "clusters/mt_rmt_shared/all_hm2.pdf", sep = "" ), height = 20, width = 10)
do_heatmap(heat.dat, title = "MT-RMT Shared")
dev.off()
grep("Ifng", rownames(heat.dat))

##ALL 
common.list <- Reduce(union, list(mt.deg.list, sbp.deg.list, rmt.deg.list))

common.dat <- assay(vsd[deg.genes,])

common.scaled <- t(scale(t(common.dat)))

common.scaled.cor <- as.dist(1-cor(t(common.scaled), method = "pearson"))

hr <- hclust(common.scaled.cor, method="complete") # Cluster Genes

clusDyn <- cutreeDynamic(hr, distM = as.matrix(common.scaled.cor), method = "hybrid")
names(clusDyn) <- rownames(common.scaled)

clus.centroids <- sapply(levels(factor(clusDyn)), clust.centroid, common.scaled, clusDyn)

centroid.melt <- melt(clus.centroids)
colnames(centroid.melt) <- c("Sample", "Cluster", "value")

for (i in 1:length(unique(centroid.melt$Cluster))){
  centroid.df <- subset(centroid.melt, Cluster == i)
  centroid.df$day <- ifelse(grepl("naive", centroid.df$Sample), 0,
                            str_match(centroid.df$Sample, pattern= "\\.d(.+)\\.")[,2])
  centroid.df$delivery <- ifelse(grepl("naive", centroid.df$Sample), "naive",
                                 str_match(centroid.df$Sample, pattern="^(.*)\\..*\\.")[,2])
  # centroid.df <- centroid.df[5:nrow(centroid.df), ]
  cols <- c("mt" = "grey50", "sbp" = "deeppink2", "naive" = "black", "rtmt" = "dodgerblue2")
  
  p1 <- ggplot(centroid.df, aes(x = day, y = value, group= delivery, colour=as.factor(delivery))) +
    geom_point() +
    # stat_summary(aes(y = value,group=1), fun.y=mean, colour="red", geom="line",group=1) + 
    geom_smooth(se = FALSE, alpha=0.1, method = "loess", span=0.9) +
    xlab("Time") +
    ylab("Expression") +
    labs(title= paste("Cluster", i, "Expression by Time", sep = " "),color = "Cluster") +
    scale_colour_manual(values = cols)
  
  pdf(paste(results.dir, "clusters/all_DE2/all_cluster_", i, ".pdf", sep =""))
  plot(p1)
  dev.off()
}

ifng.clus <- clusDyn["ENSMUSG00000055170"] # IFNG

ind = (clusDyn == ifng.clus)
clus1 <- names(clusDyn[ind])

heat.dat <-as.matrix(assay(vsd)[clus1, ])
heat.dat   <- t(apply(heat.dat,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))}))
rownames(heat.dat) <- gtf.dat[clus1,"gene_name"]

#order so days are together 
heat.list <- colnames(heat.dat);
heat.list <- as.data.frame(heat.list);
heat.list$day <- ifelse(grepl("naive", heat.list$heat.list),
                        0,
                        str_match(string=heat.list$heat.list, pattern= "\\.d(.+)\\.")[,2])
heat.list$delivery <- ifelse(grepl("naive", heat.list$heat.list),
                             "naive",
                             str_match(string=heat.list$heat.list, pattern="^(.*)\\..*\\.")[,2])
heat.list$delivery <- ordered(factor(heat.list$delivery, levels= c('naive', 'sbp', 'rtmt', 'mt')))
heat.list$rep <- as.integer(str_match(string=heat.list$heat.list, pattern=".*r(.)$")[,2])
heat.list <- heat.list[order(as.numeric(heat.list$day), as.numeric(as.factor(heat.list$delivery)), heat.list$rep),]
heat.list <- as.vector(heat.list$heat.list);
heat.dat <- heat.dat[,heat.list] 

pdf(paste(results.dir, "clusters/all_DE2/cluster_3_hm.pdf", sep = "" ), height = 20, width = 10)
do_heatmap(heat.dat, title = "Cluster 3")
dev.off()

plot.list <- list()

for (i in 1:length(unique(centroid.melt$Cluster))){
  centroid.df <- subset(centroid.melt, Cluster == i)
  centroid.df$day <- ifelse(grepl("naive", centroid.df$Sample), 0,
                            str_match(centroid.df$Sample, pattern= "\\.d(.+)\\.")[,2])
  centroid.df$delivery <- ifelse(grepl("naive", centroid.df$Sample), "naive",
                                 str_match(centroid.df$Sample, pattern="^(.*)\\..*\\.")[,2])
  # centroid.df <- centroid.df[5:nrow(centroid.df), ]
  cols <- c("mt" = "grey50", "sbp" = "deeppink2", "naive" = "black", "rtmt" = "dodgerblue2")
  
  p1 <- ggplot(centroid.df, aes(x = day, y = value, group= delivery, colour=as.factor(delivery))) +
    geom_point() +
    # stat_summary(aes(y = value,group=1), fun.y=mean, colour="red", geom="line",group=1) + 
    geom_smooth(se = FALSE, alpha=0.1, method = "loess", span=0.9) +
    xlab("Time") +
    ylab("Expression") +
    labs(title= paste("Cluster", i, "Expression by Time", sep = " "),color = "Cluster") +
    scale_colour_manual(values = cols)
    plot.list[[i]] <- assign(paste("clus_", i, '.plot', sep=''), p1);
}
library(cowplot)
pdf(paste(results.dir, "clusters/all_DE2/cluster_facet2.pdf", sep = "" ), height = 20, width = 20)
plot_grid(plotlist = plot.list)
dev.off()


# Cytokines
cytokines <- c("Ifng", "Il1a", "Il2", "Il6", "Il10", "Tnf", "Ccl2", "Csf2", "Csf3", "Cxcl1", "Cxcl10", "Cxcl2")
cyto.id <- sapply(cytokines, function(c){
  gene <- gtf.dat[gtf.dat$gene_name %in% c, "gene_id"]
  gene
})

clusDyn[names(clusDyn) %in% cyto.id]

# Cytpo.db

innate.db <- read.table(paste(data.dir, "InnateDB_genes.csv", sep = ""), sep = ",",
                        header = TRUE)

innate.genes <- unlist(sapply(innate.db$Query.Xref, function(x){
  gene.id <- gtf.dat[gtf.dat$gene_name %in% x, "gene_id"]
  gene.id
})) # name - gene_id, value - ensembl ID

innate.clus <- clusDyn[names(clusDyn) %in% innate.genes]
innate.tbl <- as.data.frame(table(innate.clus))
names(innate.tbl) <- c("Cluster", "Num_innate")

clusDyn.df <- as.data.frame(table(clusDyn))
names(clusDyn.df) <- c("Cluster", "Num_genes")

innate.clus.df <- merge(clusDyn.df, innate.tbl, by = "Cluster")
innate.clus.df$percent <- (innate.clus.df$Num_innate / innate.clus.df$Num_genes) * 100
innate.clus.df[order(innate.clus.df$percent, decreasing = TRUE),  ]


cluster.df <- as.data.frame(clusDyn)
cluster.df$gene_name <- gtf.dat[rownames(cluster.df), "gene_name"]
head(cluster.df)
write.table(cluster.df,
            paste(results.dir, "Cluster_gene_assignment.csv", sep = " "),
            col.names=TRUE,row.names=TRUE,sep=",",quote=F)
