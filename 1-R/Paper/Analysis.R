rm(list = ls())
source('/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/1-R/Paper/packages.R')
source('/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/1-R/Paper/EKD11_functions.R')


nf.dir      <- "/Users/alderc/1-projects/CAMP/1-AS_timecourse/2-EKD11/1-Pipeline/1-Nextflow/"
work.dir    <- "/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/"
r.dir       <- paste(work.dir, "1-R/Paper",sep='') #change once refactor is completed
tmp.dir     <- paste(work.dir,"tmp/",sep='')
data.dir    <- paste(work.dir,"4-data/",sep='')
results.dir <- paste(work.dir,"2-results_2020/FINAL/",sep='')
genome.dir  <- "/Users/alderc/1-projects/9-Data/1-Reference_genomes/1-Mus_musculus/"

## Files
gtf.file <- paste(genome.dir, "Mus_musculus.GRCm38.86.gtf", sep="")
design.file <- paste(nf.dir,"design_sbp_rtmt2.csv",sep='')
design.baseline <- "naive"

## Loading GTF file
gtf.dat    <- import(gtf.file)
gtf.dat    <- as.data.frame(gtf.dat[gtf.dat$type%in% "gene",])
gtf.dat$GRCm38 <- paste("chr",gtf.dat$seqnames,":",gtf.dat$start,"-",gtf.dat$end,sep='')
rownames(gtf.dat) <- gtf.dat$gene_id


## Output names
count.output <- "EKD11.normalised_counts_with_controls.xls.gz" 
log.output   <-  "2021_EKD11.vst_counts_with_controls.xls.gz"
matrix.output <- "EKD11_sbp_rmt_counts.csv"


## Design file
design           <- read.delim(design.file, header=TRUE, sep=",", stringsAsFactors = TRUE)
rownames(design) <- design$label
design$delivery  <- relevel(design$delivery, design.baseline) # sets design.baseline as first factor (baseline)
design$grp       <- relevel(design$grp, design.baseline) # same as above but for analysis with controls 
design$day       <- factor(as.character(design$day), levels=as.character(sort(unique(design$day)))) # factors days in numberical order
design$replicate <- as.factor(design$replicate) # factors replicates (order not important)

design


## Importing counts
txi <- tximport(paste(nf.dir, "results/rsem/", design$lims.name, ".genes.results", sep=""), type="rsem")
txi$length[txi$length==0] <- 1
dat <- DESeqDataSetFromTximport(txi, design, ~ grp) # For baseline counts - full rank in other file
rowData(dat) <- gtf.dat[rownames(dat),c("gene_id","gene_name","gene_source","gene_biotype","GRCm38")]
dds <- DESeq(dat)


## Filter
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 4
dds_filter <- dds[filter,]

#Count matrices
dds_filter <-DESeq(dds_filter)
vsd <- vst(dds_filter, blind = TRUE)

vsd.d6 <- vsd[,  grep("naive|d6", colnames(vsd))]
vsd.late <- vsd[,  grep("naive|d4|d6", colnames(vsd))]

## Outputs
# Counts Dataframe (Normalised)
count.norm  <- counts(dds_filter,normalized=TRUE)
count.dBase <- data.frame(count.norm,rowData(dds_filter),stringsAsFactors=FALSE)
rownames(count.dBase) <- count.dBase$gene_id
if(!file.exists(paste(results.dir,count.output ,sep=''))){
  write.table(count.dBase,file=gzfile(paste(results.dir,count.output ,sep='')),col.names=TRUE,row.names=FALSE,sep="\t",quote=F)}


# Log counts (vst)
vst.dBase <- assay(vsd)
vst.dBase <- data.frame(vst.dBase,rowData(vsd),stringsAsFactors=FALSE)
if(!file.exists(paste(results.dir,log.output ,sep=''))){
  write.table(vst.dBase,file=gzfile(paste(results.dir,log.output,sep='')),col.names=TRUE,row.names=FALSE,sep="\t",quote=F)
  }

# Count Matrix (un-normalised)
count.mat <- counts(dds, normalized=FALSE)
# if(!file.exists(paste(results.dir, matrix.output ,sep=''))){
  # write.table(count.mat, file=gzfile(paste(results.dir, matrix.output, sep='')), col.names=TRUE, row.names=TRUE, sep=',',quote=FALSE)}

# TPM 
tpm.db <- as.data.frame(txi$abundance)
names(tpm.db) <- names(as.data.frame(count.norm))



## PCA 
rv <- rowVars(assay(vsd)) 
select <- order(rv, decreasing = TRUE)[1:1000]

pca <- prcomp(t(assay(vsd)[select, ]))
pca.sum <- summary(pca)
pca.val <- pca.sum$importance[2,1:3]
pca.df <- as.data.frame(pca$x[,1:3])
pca.df <- cbind(pca.df, design[rownames(pca.df), c("delivery", "day", "grp", "label")])
# PCA Colours
pca.cols <- c("blood" = "#E377C2", "naive" = "black", 'recent_blood' = "#17BECF")


pca.plot <- plot_ly(data = pca.df) %>%
  add_trace(type = "scatter",
            x = ~PC1, y = ~PC2,
            color = ~delivery, colors = pca.cols,
            symbols = c("0" = "triangle-up" ,"1" = "circle", "2" = "cross", "3" = "x", "4" = "square", "6" = "diamond"),
            text = ~label,
            mode = "markers",
            showlegend = F) %>% 
  add_trace(type = "scatter",
            x = ~PC1, y = ~PC2,
            symbol = ~day,
            mode = "markers",
            text = ~label,
            marker = list(color = "#17BECF", size = 8)) %>% 
  add_trace(type = "scatter",
            x = ~PC1, y = ~PC2,
            symbol = ~day,
            mode = "markers",
            text = ~label,
            marker = list(color = "#E377C2", size = 8)) %>% 
  add_trace(type = "scatter",
            x = ~PC1, y = ~PC2,
            color = ~delivery,
            symbol = ~day,
            mode = "markers",
            text = ~label,
            marker = list(size = 14),
            showlegend = F) %>% 
  layout(xaxis = list(title = paste("PC1:", round(pca.val[[1]], 2)*100, "%", sep = " "), titlefont = list(size = 25), zeroline = F, showline = T, mirror = "ticks"),
         yaxis = list(title = paste("PC2:", round(pca.val[[2]], 2)*100, "%", sep = " "), titlefont = list(size = 25), zeroline = F, showline = T, mirror = "ticks"),
         legend = list(legend_title_text = "Day"))

pca.plot
# Exporting File
Sys.setenv("PATH" = paste(Sys.getenv("PATH"), "/anaconda3/bin", sep = .Platform$path.sep)) ## To run orca
orca(pca.plot, file="Paper_PCA.pdf", format = "pdf")# Orca will append to WD (no absolute path)


## Pairwise Comparison 
res.list <- list()

res.list[["rtmt.1_vs_naive"]]      <- results(dds_filter, contrast=c("grp","rtmt.1","naive"))
res.list[["rtmt.2_vs_naive"]]      <- results(dds_filter, contrast=c("grp","rtmt.2","naive"))
res.list[["rtmt.3_vs_naive"]]      <- results(dds_filter, contrast=c("grp","rtmt.3","naive"))
res.list[["rtmt.4_vs_naive"]]      <- results(dds_filter, contrast=c("grp","rtmt.4","naive"))
res.list[["rtmt.6_vs_naive"]]      <- results(dds_filter, contrast=c("grp","rtmt.6","naive"))

res.list[["sbp.1_vs_naive"]]      <- results(dds_filter, contrast=c("grp","sbp.1","naive"))
res.list[["sbp.2_vs_naive"]]      <- results(dds_filter, contrast=c("grp","sbp.2","naive"))
res.list[["sbp.3_vs_naive"]]      <- results(dds_filter, contrast=c("grp","sbp.3","naive"))
res.list[["sbp.4_vs_naive"]]      <- results(dds_filter, contrast=c("grp","sbp.4","naive"))
res.list[["sbp.6_vs_naive"]]      <- results(dds_filter, contrast=c("grp","sbp.6","naive"))

res.list[['rtmt.1_vs_sbp.1']]      <- results(dds_filter, contrast=c("grp","rtmt.1","sbp.1"));
res.list[['rtmt.2_vs_sbp.2']]      <- results(dds_filter, contrast=c("grp","rtmt.2","sbp.2"));
res.list[['rtmt.3_vs_sbp.3']]      <- results(dds_filter, contrast=c("grp","rtmt.3","sbp.3"));
res.list[['rtmt.4_vs_sbp.4']]      <- results(dds_filter, contrast=c("grp","rtmt.4","sbp.4"));
res.list[['rtmt.6_vs_sbp.6']]      <- results(dds_filter, contrast=c("grp","rtmt.6","sbp.6"));



res.list <- lapply(res.list,function(res){
  res <- as.data.frame(res)
  colnames(res)[colnames(res) %in% "log2FoldChange"] <- "log2FC"
  colnames(res)[colnames(res) %in% "padj"] <- "FDR"
  res$DEG <- res$FDR <= 0.05 & ### Changing threshold for log2FC
    abs(res$log2FC) > log2(2) &
    res$baseMean>=30         &
    !is.na(res$log2FC)       &
    !is.na(res$FDR)
  res
})

cbind(sapply(res.list, function(x){sum(x$DEG)}))

deg.list <- lapply(res.list, function(df){
  df <- df[df$DEG, ]
  df$gene_name <- gtf.dat[rownames(df), "gene_name"]
  df
})


library(xlsx)
for (i in names(deg.list)){
  df = deg.list[[i]]
  if (nrow(df) != 0){
    df <- df[, c(8, 1:7)]
    write.xlsx(df, file="EKD11_DEGs_Table.xlsx",
               row.names = F,
               sheetName=i,
               append=T)
  }
}
# Deg by transmission
rvn.deg <- deg.list[1:5]
svn.deg <- deg.list[6:10]

rvn.genes <- unique(unlist(sapply(rvn.deg, row.names)))
svn.genes <- unique(unlist(sapply(svn.deg, row.names)))

rmt_unique <- setdiff(rvn.genes, svn.genes)
sbp_unique <- setdiff(svn.genes, rvn.genes)

## Supplementary heatmap as matrix
#RMT
rmt_unique.dat <- assay(vsd)[rmt_unique, ]
rmt_unique.dat   <- t(apply(rmt_unique.dat,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))}))
rownames(rmt_unique.dat) <- gtf.dat[rownames(rmt_unique.dat),"gene_name"]

write.table(rmt_unique.dat, "RMT_unique_z-scores.csv", sep = ",")


#SBP
sbp_unique.dat <- assay(vsd)[sbp_unique, ]
sbp_unique.dat   <- t(apply(sbp_unique.dat,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))}))
rownames(sbp_unique.dat) <- gtf.dat[rownames(sbp_unique.dat),"gene_name"]

write.table(sbp_unique.dat, "SBP_unique_z-scores.csv", sep = ",")


## Venn Diagram 
library("VennDiagram")
venn.cols <- c("#E377C2", "#17BECF")
venn <- venn.diagram(x = list(SBP = svn.genes, RMT = rvn.genes),
                     filename = NULL,
                     
                     #Circles
                     lwd = 2,
                     lty = "blank",
                     fill = venn.cols,
                     
                     # Numbers
                     fontface = "bold",
                     fontfamily = "sans",
                     cex = 2.5,
                     
                     # Set names
                     cat.fontface = "bold",
                     cat.default.pos = "outer",
                     cat.fontfamily = "sans",
                     # cat.dist = c(0.05, 0.05),
                     # cat.pos = c(20, -20),
                     cat.cex = 2.5)

pdf(file = "Venn_Diagram_Overview_Revision.pdf")
grid.newpage()
grid.draw(venn)
dev.off()

# Day by day
for (i in 1:5){
  s <- rownames(svn.deg[[i]])
  r <- rownames(rvn.deg[[i]])
  ifelse(length(r) < length(s), y <- 180, y <- 0)
  ifelse(i == 5, x <- 6, x <- i)
  venn <- venn.diagram(x = list(SBP = s, RMT = r),
                       filename = NULL,
                       
                       #Circles
                       lwd = 2,
                       lty = "blank",
                       fill = venn.cols,
                       
                       # Numbers
                       fontface = "bold",
                       fontfamily = "sans",
                       cex = 2.5,
                       rotation.degree = y,
                       
                       # Set names
                       cat.fontface = "bold",
                       cat.default.pos = "outer",
                       cat.fontfamily = "sans",
                       # cat.dist = c(0.05, 0.05),
                       # cat.pos = c(20, -20),
                       cat.cex = 0
  )
  
  pdf(file = paste("Venn_Day_Revision", x, ".pdf", sep = ""))
  grid.newpage()
  grid.draw(venn)
  dev.off()
}

# Innate DB 
innate.sym <- readxl::read_xls(paste(work.dir,"/4-data/innatedb_curated_genes.xls", sep = "")) %>% 
  subset(Species == 10090) %>% select("Gene Symbol") %>% unique()

innate.id <- unlist(sapply(innate.sym, function(x){gtf.dat[gtf.dat$gene_name %in% x, "gene_id"]}))[,1]

rvn.deg.late <- rvn.deg[4:5]
svn.deg.late <- svn.deg[4:5]


### For Paper
rvn.late.genes <- unique(unlist(sapply(rvn.deg.late, row.names)))
svn.late.genes <- unique(unlist(sapply(svn.deg.late, row.names)))

combined <- unique(c(rvn.late.genes, svn.late.genes))

combined_innate <- intersect(innate.id, combined)


pdf("innate_deg_4_6_HM_revision.pdf", height = 14, width = 5)
split.heatmap2(gene_list = combined_innate, heat.dat = assay(vsd.late), row_names = T  )
dev.off()



### TPM Plots

tpm.late <- tpm.db[ ,grep("d4|d6|naive", names(tpm.db))]

names(tpm.late)

figure.genes.sym <- c("Ifng","Cxcl9", "Cxcl10", "Ccl2",
                      "Il21","Nos2", "Lgals3", "Clec7a",
                      "Il2rb", "Tnfrsf9", "Cxcr3", "Irf1",
                      "Ifit2", "Igtp", "Gbp10", "Msr1")


figure.gene.id <- gtf.dat[gtf.dat$gene_name %in% figure.genes.sym, "gene_id"]
dir.create("./TPM_Figures/")
for (i in figure.gene.id){
  gene_name <- gtf.dat[i, "gene_name"]
  # pdf(paste("../../2-results_2020/FINAL/TPM_Figures/", gene_name, "_TPM.pdf", sep = ""))
  tpm.plot(gene_id = i, tpm.data = tpm.late)
  # dev.off()
}


## LRT
design.rvs <- design[-grep("naive", rownames(design)),]
design.rvs$delivery <- factor(design.rvs$delivery)
design.rvs$day       <- factor(as.character(design.rvs$day), levels=as.character(sort(unique(design.rvs$day)))) # factors days in numberical order
design.rvs$replicate <- as.factor(design.rvs$replicate)
design.rvs$grp <- factor(design.rvs$grp)# factors replicates (order not important)

## Importing counts
txi.rvs <- tximport(paste(nf.dir, "results/rsem/", design.rvs$lims.name, ".genes.results", sep=""), type="rsem")
txi.rvs$length[txi.rvs$length==0] <- 1


dat.rvs <- DESeqDataSetFromTximport(txi.rvs, design.rvs, ~ delivery + day + delivery:day) # For baseline counts - full rank in other file
rowData(dat.rvs) <- gtf.dat[rownames(dat.rvs),c("gene_id","gene_name","gene_source","gene_biotype","GRCm38")]
dds.rvs <- DESeq(dat.rvs)


nc <- counts(dds.rvs, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 4
dds.rvs <- dds.rvs[filter,]
dds.rvs <-DESeq(dds.rvs)

vsd.rvs <- vst(dat.rvs, blind = TRUE)

rvs.lrt <- DESeq(dds.rvs, test="LRT", reduced = ~ delivery + day) ## Showing genes with any different across time

lrt.res <- as.data.frame(results(rvs.lrt, alpha = 0.05))
lrt.sig <- lrt.res[lrt.res$padj < 0.05, ] %>% drop_na()
lrt.sig$gene_name <- gtf.dat[rownames(lrt.sig), "gene_name"]
lrt.sig$gene_name
lrt.sig

#HM 
pdf(paste(results.dir, "HM_Supplementary7.pdf", sep = ""), height = 14, width = 7)
split.heatmap2(gene_list = rownames(lrt.sig), heat.dat = assay(vsd.d6), row_names = TRUE)
dev.off()


for(g in rownames(lrt.sig)){
  gene_name <- lrt.sig[g, "gene_name"]
  pdf(paste("LRT_TPM/", gene_name, "_TPM.pdf", sep = ""))
  tpm.plot(gene_id = g, tpm.db = tpm.db)
  dev.off()
}


rvs3.res <- res.list[[13]]
rvs4.res <- res.list[[14]]
rvs6.res <- res.list[[15]]


rvs3.fdr <- rvs3.res[rvs3.res$FDR < 0.05, ] %>% drop_na())
rvs4.fdr <- rvs4.res[rvs4.res$FDR < 0.05 & abs(rvs4.res$log2FC) > 0.6 , ] %>% drop_na()
rvs6.fdr <- rvs6.res[rvs6.res$FDR < 0.05 & abs(rvs6.res$log2FC) > 0.6 , ] %>% drop_na()

rvs.fdr <- unique(c(rownames(rvs4.fdr), rownames(rvs6.fdr)))

gtf.dat[rvs.fdr, "gene_name"]
tpm.plot(rvs.fdr[50])
split.heatmap2(gene_list = rvs.f, heat.dat = assay(vsd.late), row_names = TRUE)


### OLD DATASET - Justification for removing d4 samples
# load("../../1-R/Paper/deseq.Rdata")

library(dendextend)
# Changing the labels for figure purposes
c.assay <- assay(vsd)
labels <- vsd$label.short #change to vsd tomorrow
labels <- toupper(substring(labels, 3))
colnames(c.assay) <- labels

# Dissimilarity matrix (pearsons)
c <- cor(c.assay[select, ], method = 'pearson')
c.dist <- as.dist(1-c)
# dist.mat <- dist(t(assay(vsd)), method = "euclidean")
hr.cor <- hclust(c.dist, method = "ward.D2", members=NULL)

# Dendrogram
dend.cor <- (as.dendrogram(hr.cor))
delivery <- vsd$delivery
cols <- c("black", "deeppink2", "dodgerblue") #order in levels of delivery
col_grp <- cols[delivery]
col_grp <- col_grp[order.dendrogram(dend.cor)]
pdf("sample_correlation_dendrogram_colour.pdf", width = 10, height = 10)
par(mar = c(3,1,1,7))
dend.cor  %>% 
  set("labels_colors", col_grp) %>% #change label colors to GROUP
  plot(main = "Dendrogram of sample correlation", horiz = TRUE)
legend("topleft", legend = levels(delivery), fill = cols, cex = 0.75)
dev.off()



