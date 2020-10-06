### Tradeseq gam test
library(stringr)
library(clusterExperiment)
library(tradeSeq)
library(cowplot)

counts.mat <- counts(dds.fr, normalized = FALSE)
counts.mat <- counts.mat[rownames(deg.lrt) ,1:60]
# counts.mat <- counts.mat[rt.mt.ids, 1:60] #RT:MT shared genes
head(counts)
time <- as.numeric(str_match(string=colnames(counts.mat), pattern= "\\.d(.+)\\.")[,2])
treatment <- factor(str_match(string=colnames(counts.mat), pattern="^(.*)\\..*\\.")[,2], levels = c('sbp', 'mt', 'rtmt'))

# filter counts: at least 5 counts in 3 replicats
keep <- rowSums(counts.mat>5)>=3
counts.mat <- counts.mat[keep,]

time <- matrix(time, nrow=ncol(counts.mat), ncol=3, byrow=FALSE)
rownames(time) <- colnames(counts.mat)
weights <- matrix(0, nrow=ncol(counts.mat), ncol=3)
weights[1:20, 1] <- 1 
weights[21:40, 2] <- 1
weights[41:60, 3] <- 1

## evaluate optimal K
infMat <- evaluateK(counts.mat, pseudotime=time, cellWeights=weights, nGenes=255, k=3:5)

## fit GAM
count.lrt <- counts.mat[rownames(deg.lrt), ]
gamList <- fitGAM(count.lrt, pseudotime=time, cellWeights=weights, nknots=5)

resPat <- patternTest(gamList, nPoints=20)

oPat <- order(resPat$waldStat, decreasing = TRUE)

sum(p.adjust(resPat$pvalue, "fdr") <= 0.01, na.rm=TRUE)

deg.gams <-rownames(count.lrt)[p.adjust(resPat$pvalue, "fdr") <= 0.01]
deg.gams <- intersect(deg.gams, names(gamList))

test <- gamList[[2]]

plotSmoothers(test)

# ## test for different expression pattern
# resPat <- patternTest(gamList, nPoints=8)

#Clusters
nPointsClus <- 6
# clusPat <- clusterExpressionPatterns(gamList, genes = rownames(counts.mat), nPoints = nPointsClus) ## RTMT_genes
clusPat <- clusterExpressionPatterns(gamList, genes = rownames(count.lrt), nPoints = nPointsClus)


clusLabels <- primaryCluster(clusPat$rsec)

cUniq <- unique(clusLabels)
cUniq <- cUniq[!cUniq == -1]
cUniq <- order(cUniq)# Removes unclustered genes

# pca.cols <- c("mosquito" = "#BCBD22", "blood" = "#E377C2", "naive" = "black", 'recent_blood' = "#17BECF")
plots <- list()
for (xx in cUniq) {
  cId <- which(clusLabels == xx)
  p <- ggplot(data = data.frame(x = 1:nPointsClus,
                                y = rep(range(clusPat$yhatScaled[cId, ]),
                                        nPointsClus / 2)),
              aes(x = x, y = y)) +
    geom_point(alpha = 0) +
    labs(title = paste0("Cluster ", xx),  x = "Pseudotime", y = "Normalized expression") +
    theme_classic()
  for (ii in 1:length(cId)) {
    geneId <- rownames(clusPat$yhatScaled)[cId[ii]]
    p <- p +
      geom_line(data = data.frame(x = rep(1:nPointsClus, 3),
                                  y = clusPat$yhatScaled[geneId, ],
                                  lineage = rep(0:2, each = nPointsClus)),
                aes(col = as.character(lineage), group = lineage), lwd = 1.5)
  }
  p <- p + guides(color = FALSE) +
    scale_color_manual(values = c("#E377C2", "#BCBD22", "#17BECF"),
                       breaks = c("0", "1", "2"))  
  plots[[as.character(xx)]] <- p
}
plots$ncol <- 5
pdf(paste(results.dir, "EKD11_gam_cluster_LRT_deg_6points.pdf", sep = ""), height = 160, width = 40)
do.call(plot_grid, plots)
dev.off()

##################################################################################################################################
###RTMT: MT
rt.mt.counts <- counts.mat[rt.mt.ids, ]
## fit GAM

rt.mt.gamList <- fitGAM(rt.mt.counts, pseudotime=time, cellWeights=weights, nknots=5)

rt.mt.resPat <- patternTest(gamList, nPoints=6)

oPat <- order(rt.mt.resPat$waldStat, decreasing = TRUE)

head(rownames(rt.mt.resPat)[oPat])

sum(p.adjust(resPat$pvalue, "fdr") <= 0.01, na.rm=TRUE)

deg.gams <-rownames(rt.mt.counts)[p.adjust(rt.mt.resPat$pvalue, "fdr") <= 0.01]
deg.gams <- intersect(deg.gams, names(rt.mt.gamList))

test <- gamList[[2]]

p <- plotSmoothers(test) + ggtitle(names(test))
gam.plots <- list ()

# for (i in rt.mt.gamList){
#   name <- 
#   p <- plotSmoothers(test[[1]]) + 
#     ggtitle(names(test)) + 
#     scale_color_manual(values = c("#E377C2", "#BCBD22", "#17BECF"),
#                        name = "Transmission Route",
#                        labels=c("SBP", "MT", "RTMT"))
#   
# }
# plotSmoothers(test)

#Clusters
nPointsClus <- 6
rt.mt.clusPat <- clusterExpressionPatterns(rt.mt.gamList, genes = rownames(rt.mt.counts), nPoints = nPointsClus)


rt.mt.clusLabels <- primaryCluster(rt.mt.clusPat$rsec)
names(rt.mt.clusLabels) <- rownames(rt.mt.counts)

rt.mt.cUniq <- unique(rt.mt.clusLabels)
rt.mt.cUniq <- rt.mt.cUniq[!rt.mt.cUniq == -1] %>% sort() # Removes unclustered genes


# pca.cols <- c("mosquito" = "#BCBD22", "blood" = "#E377C2", "naive" = "black", 'recent_blood' = "#17BECF")
plots <- list()
for (xx in rt.mt.cUniq) {
  cId <- which(rt.mt.clusLabels == xx)
  p <- ggplot(data = data.frame(x = 1:nPointsClus,
                                y = rep(range(rt.mt.clusPat$yhatScaled[cId, ]),
                                        nPointsClus / 2)),
              aes(x = x, y = y)) +
    geom_point(alpha = 0) +
    labs(title = paste0("Cluster ", xx),  x = "Pseudotime", y = "Normalized expression") +
    theme_classic()
  for (ii in 1:length(cId)) {
    geneId <- rownames(rt.mt.clusPat$yhatScaled)[cId[ii]]
    p <- p +
      geom_line(data = data.frame(x = rep(1:nPointsClus, 3),
                                  y = rt.mt.clusPat$yhatScaled[geneId, ],
                                  lineage = rep(0:2, each = nPointsClus)),
                aes(col = as.character(lineage), group = lineage), lwd = 1.5)
  }
  p <- p + guides(color = FALSE) +
    scale_color_manual(values = c("#E377C2", "#BCBD22", "#17BECF"),
                       breaks = c("0", "1", "2"))  
  plots[[as.character(xx)]] <- p
}
plots$ncol <- 4
pdf(paste(results.dir, "EKD11_rt_mt_shared_cluster_ALL_6points.pdf", sep = ""), height = 20, width = 20)
do.call(plot_grid, plots)
dev.off()

rt.mt.clus.df <- as.data.frame(rt.mt.clusLabels)
rt.mt.clus.df$gene_name <- gtf.dat[rownames(rt.mt.clus.df), "gene_name"]
write.table(rt.mt.clus.df, file= paste(results.dir, "rt_mt_cluster_df.csv", sep = ""), quote = F, row.names = T, sep = ",")




###################################################
###################################################
## Modular GAMS

module.anno <- read.csv(file=paste(data.dir, '/modular_analysis/Annotation of Blood modules.csv', sep=''),
                        header=TRUE,
                        sep=',');
module.dat <- read.csv(file=paste(data.dir, '/modular_analysis/Blood modules.csv', sep=''),
                       header=TRUE,
                       sep=',');

module.num <- module.anno$X

module.geneset <- lapply(module.num, function(mod){
  dat <- module.dat[module.dat$Module == mod, 'X']; 
  dat <- dat[!is.na(dat)]
  mod <- dat
})

module.geneset.names <- lapply(module.num, function(mod){
  dat <- module.dat[module.dat$Module == mod, 'Gene.name']; 
  dat <- dat[!is.na(dat)]
  mod <- dat
})

names(module.geneset) <- paste(module.num, module.anno$Biological.process)
names(module.geneset.names) <- paste(module.num, module.anno$Biological.process)

names(module.geneset.names)
module.df <- data.frame(module = names(module.geneset.names))
module.df$genes <- module.geneset.names
module.df$genes <- vapply(module.df$genes, paste, collapse = ", ", character(1L))

# write.table(module.df, file = paste(results.dir, "module_gene_df.tsv", sep = ""), sep = "\t", row.names = F, quote = F)

counts.mat <- counts(dds, normalized = FALSE)
counts.mat <- counts.mat[ ,1:60]

ifn.il10.ids <- module.geneset[[11]]

ifn.il10.mat <- counts.mat[ifn.il10.ids, ]


## fit GAM
gamList <- fitGAM(ifn.il10.mat, pseudotime=time, cellWeights=weights, nknots=5)

resPat <- patternTest(gamList, nPoints=6)

oPat <- order(resPat$waldStat, decreasing = TRUE)

sum(p.adjust(resPat$pvalue, "fdr") <= 0.01, na.rm=TRUE)

deg.gams <-rownames(ifn.il10.mat)[p.adjust(resPat$pvalue, "fdr") <= 0.01]
deg.gams <- intersect(deg.gams, names(gamList))

test <- gamList[[2]]

plotSmoothers(test)

# ## test for different expression pattern
# resPat <- patternTest(gamList, nPoints=8)

#Clusters
nPointsClus <- 6
# clusPat <- clusterExpressionPatterns(gamList, genes = rownames(counts.mat), nPoints = nPointsClus) ## RTMT_genes
clusPat <- clusterExpressionPatterns(gamList, genes = rownames(ifn.il10.mat), nPoints = nPointsClus)


clusLabels <- primaryCluster(clusPat$rsec)

cUniq <- unique(clusLabels)
cUniq <- cUniq[!cUniq == -1]
cUniq <- sort(cUniq)# Removes unclustered genes


# pca.cols <- c("mosquito" = "#BCBD22", "blood" = "#E377C2", "naive" = "black", 'recent_blood' = "#17BECF")
plots <- list()
for (xx in cUniq) {
  cId <- which(clusLabels == xx)
  p <- ggplot(data = data.frame(x = 1:nPointsClus,
                                y = rep(range(clusPat$yhatScaled[cId, ]),
                                        nPointsClus / 2)),
              aes(x = x, y = y)) +
    geom_point(alpha = 0) +
    labs(title = paste0("Cluster ", xx),  x = "Pseudotime", y = "Normalized expression") +
    theme_classic()
  for (ii in 1:length(cId)) {
    geneId <- rownames(clusPat$yhatScaled)[cId[ii]]
    p <- p +
      geom_line(data = data.frame(x = rep(1:nPointsClus, 3),
                                  y = clusPat$yhatScaled[geneId, ],
                                  lineage = rep(0:2, each = nPointsClus)),
                aes(col = as.character(lineage), group = lineage), lwd = 1.5)
  }
  p <- p + guides(color = FALSE) +
    scale_color_manual(values = c("#E377C2", "#BCBD22", "#17BECF"),
                       breaks = c("0", "1", "2"))  
  plots[[as.character(xx)]] <- p
}
plots$ncol <- 4
do.call(plot_grid, plots)


##### Module B14

ifn.sig.ids <- module.geneset[[14]]
ifn.sig.mat <- counts.mat[ifn.sig.ids, ]

gamList <- fitGAM(ifn.sig.mat, pseudotime=time, cellWeights=weights, nknots=5)

resPat <- patternTest(gamList, nPoints=20)

oPat <- order(resPat$waldStat, decreasing = TRUE)

sum(p.adjust(resPat$pvalue, "fdr") <= 0.01, na.rm=TRUE)

deg.gams <-rownames(ifn.sig.mat)[p.adjust(resPat$pvalue, "fdr") <= 0.01]
deg.gams <- intersect(deg.gams, names(gamList))

#Clusters
nPointsClus <- 6
# clusPat <- clusterExpressionPatterns(gamList, genes = rownames(counts.mat), nPoints = nPointsClus) ## RTMT_genes
clusPat <- clusterExpressionPatterns(gamList, genes = deg.gams, nPoints = nPointsClus)


clusLabels <- primaryCluster(clusPat$rsec)

cUniq <- unique(clusLabels)
cUniq <- cUniq[!cUniq == -1]
cUniq <- sort(cUniq)# Removes unclustered genes


# pca.cols <- c("mosquito" = "#BCBD22", "blood" = "#E377C2", "naive" = "black", 'recent_blood' = "#17BECF")
plots <- list()
for (xx in cUniq) {
  cId <- which(clusLabels == xx)
  p <- ggplot(data = data.frame(x = 1:nPointsClus,
                                y = rep(range(clusPat$yhatScaled[cId, ]),
                                        nPointsClus / 2)),
              aes(x = x, y = y)) +
    geom_point(alpha = 0) +
    labs(title = paste0("Cluster ", xx),  x = "Pseudotime", y = "Normalized expression") +
    theme_classic()
  for (ii in 1:length(cId)) {
    geneId <- rownames(clusPat$yhatScaled)[cId[ii]]
    p <- p +
      geom_line(data = data.frame(x = rep(1:nPointsClus, 3),
                                  y = clusPat$yhatScaled[geneId, ],
                                  lineage = rep(0:2, each = nPointsClus)),
                aes(col = as.character(lineage), group = lineage), lwd = 1.5)
  }
  p <- p + guides(color = FALSE) +
    scale_color_manual(values = c("#E377C2", "#BCBD22", "#17BECF"),
                       breaks = c("0", "1", "2"))  
  plots[[as.character(xx)]] <- p
}
plots$ncol <- 3
do.call(plot_grid, plots)
