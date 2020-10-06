rm(list = ls())
source('/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/1-R/refactor/packages.R')


### Functions
count.plot <- function(dds_object, gene_id){
  gene_name <- gtf.dat[gene_id, "gene_name"]
  cnt.cols <- c("mosquito" = "#BCBD22", "blood" = "#E377C2", "naive" = "Black", 'recent_blood' = "#17BECF")
  cnt.dat <- plotCounts(dds_object, gene_id,
                        intgroup = c("day","delivery"), returnData = TRUE)
  # print(cnt.dat)
  cnt.trans <- cnt.dat[1:40,]
  cnt.naive <- cnt.dat[41:44, ]
  naive.mean <- mean(cnt.naive$count)
  cnt.plot <- ggplot(cnt.trans,
                     aes(x = day, y = count, color = delivery, group = delivery)) + 
    ggtitle(gene_name) + 
    geom_point() + stat_summary(fun=mean, geom="line") + 
    geom_hline(yintercept = naive.mean, col = "black", linetype="dashed") +
    scale_x_discrete(limits=factor(1:6)) + 
    scale_colour_manual(values = cnt.cols)
  return(cnt.plot)
}

tpm.plot <- function(gene_id){
  gene_name <- gtf.dat[gene_id, "gene_name"]
  cnt.cols <- c("mosquito" = "#BCBD22", "blood" = "#E377C2", "naive" = "Black", 'recent_blood' = "#17BECF")
  tpm.df <- tpm.db[gene_id, ]
  tpm.df <- melt(as.matrix(tpm.df))
  tpm.df<- cbind(tpm.df ,design[tpm.df$Var2, c("delivery", "day")])
  tpm.plot <- ggplot(tpm.df,
                     aes(x = day, y = value, group = delivery, fill = delivery)) + 
    ggtitle(gene_name) + 
    geom_bar(position =  position_dodge2(preserve = "single", padding = 0), stat = "summary", fun= "mean", width = 0.75) + 
    # geom_hline(yintercept = naive.mean, col = "black", linetype="dashed") +
    # scale_x_discrete(limits=factor(1,2,3,4,6)) + 
    scale_fill_manual("Delivery", values = cnt.cols, labels = c("Naive", "SBP", "RMT")) + 
    scale_x_discrete("Day", labels= c("Control", 1,2,3,4,6)) +
    scale_y_continuous("TPM")
  plot(tpm.plot)
}

heatmap.draw <- function(gene_list, heat.dat = assay(vsd), row_names = FALSE, cluster_samples = FALSE){
  sample.order <- design[colnames(heat.dat), ] %>% 
    arrange(delivery, day) %>% 
    pull(label) %>% as.vector()
  
  heat.genes <- gene_list
  heat.dat <- heat.dat[heat.genes,sample.order ]
  heat.dat   <- t(apply(heat.dat,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))})) %>% data.frame() %>% drop_na() %>% as.matrix()
  rownames(heat.dat) <- gtf.dat[rownames(heat.dat),"gene_name"]
  

  heat.dat <- heat.dat[,sample.order]
  max <- max(heat.dat)
  min <- min(heat.dat)
  col.scale <- colorRamp2(seq(min, max, length.out = 100), rev(colorRampPalette(brewer.pal(7, "RdYlBu"))(100)))
  n <- ncol(heat.dat)
  solid.lines <- c(0,4,24,44) # Solid lines for decorate heatmap
  ht <- Heatmap(heat.dat,
                name = "zscore",
                show_row_names          = row_names,
                show_column_names       = TRUE,
                column_names_side       = "bottom",
                show_column_dend        = cluster_samples,
                cluster_columns         = cluster_samples,
                cluster_rows            = TRUE, 
                rect_gp                 = gpar(col = "white", lwd = 1),
                row_dend_width          = unit(30, "mm"),
                column_dend_height      = unit(30, "mm"),
                show_heatmap_legend     = TRUE,
                clustering_distance_rows = 'euclidean',
                # clustering_method_rows = "complete",
                column_names_max_height = unit(8, "cm"),
                row_names_gp            = gpar(fontsize = 8),
                column_names_gp         = gpar(fontsize = 12),
                col                     = col.scale);
                # col                     = colorRamp2(c(min, 0, max), rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(3))));
  
  draw(ht)
  
  decorate_heatmap_body("zscore", c(for(d in seq(0, n, 1)){
    ifelse(d %in% solid.lines, l.type <- 'solid', l.type <- 'dashed')
    grid.lines(x = c(d/n, d/n), y = c(0,1), gp= gpar(col = "black", lty='solid', lwd=2))},
    grid.lines(x = c(0,1), y = c(0,0), gp= gpar(col = "black", lty="solid", lwd=2)),
    grid.lines(x = c(0,1), y = c(1,1), gp= gpar(col = "black", lty="solid", lwd=2))))
}

split.heatmap <- function(gene_list, heat.dat = assay(vsd), z.both = TRUE, row_names = FALSE, cluster_samples = FALSE){
  #order so days are together 
  sample.order <- design[colnames(heat.dat), ] %>% 
    arrange(delivery, day) %>% 
    pull(label) %>% as.vector()
  
  heat.genes <- intersect(gene_list, rownames(heat.dat))
  heat.dat <- heat.dat[heat.genes,sample.order ]
  if(z.both == TRUE){
    heat.dat   <- t(apply(heat.dat,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))})) %>% data.frame() %>% drop_na() %>% as.matrix()}
  rownames(heat.dat) <- gtf.dat[rownames(heat.dat),"gene_name"]

  heat.rmt <- heat.dat[, -grep("sbp", colnames(heat.dat))]
  heat.sbp <- heat.dat[, -grep("rtmt", colnames(heat.dat))]
  # heat.sbp <- heat.sbp - rowMeans(heat.sbp)
  if(z.both == FALSE){
    heat.sbp   <- t(apply(heat.sbp,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))})) %>% data.frame() %>% drop_na() %>% as.matrix()
    heat.rmt   <- t(apply(heat.rmt,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))})) %>% data.frame() %>% drop_na() %>% as.matrix()}
  

  common.max <- max(c(heat.rmt, heat.sbp))
  common.min <- min(c(heat.rmt, heat.sbp))
  
  col.scale <- colorRamp2(seq(common.min, common.max, length.out = 100), rev(colorRampPalette(brewer.pal(7, "RdYlBu"))(100)))
  # solid.lines <- c(0,4,24,44) # Solid lines for decorate heatmap
  ht.rmt <- Heatmap(heat.rmt,
                name = "zscore",
                show_row_names          = FALSE,
                show_column_names       = TRUE,
                column_names_side       = "bottom",
                cluster_columns         = cluster_samples,
                show_column_dend        = cluster_samples,
                cluster_rows            = TRUE, 
                rect_gp                 = gpar(col = "white", lwd = 1),
                row_dend_width          = unit(30, "mm"),
                column_dend_height      = unit(30, "mm"),
                show_heatmap_legend     = TRUE,
                clustering_distance_rows = 'euclidean',
                # clustering_method_rows = "complete",
                column_names_max_height = unit(8, "cm"),
                row_names_gp            = gpar(fontsize = 8),
                column_names_gp         = gpar(fontsize = 12),
                # split = kclus$cluster,
                col                     = col.scale);
                # col                     = colorRamp2(c(min, 0, max), rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(3))));
  
  
  
  

  ht.sbp <- Heatmap(heat.sbp,
                    name = "sbp",
                    show_row_names          = row_names,
                    show_column_names       = TRUE,
                    column_names_side       = "bottom",
                    cluster_columns         = cluster_samples,
                    show_column_dend        = cluster_samples,
                    cluster_rows            = FALSE, 
                    rect_gp                 = gpar(col = "white", lwd = 1),
                    row_dend_width          = unit(30, "mm"),
                    column_dend_height      = unit(30, "mm"),
                    show_heatmap_legend     = FALSE,
                    clustering_distance_rows = 'euclidean',
                    # clustering_method_rows = "complete",
                    column_names_max_height = unit(8, "cm"),
                    row_names_gp            = gpar(fontsize = 8),
                    column_names_gp         = gpar(fontsize = 12),
                    # split = kclus$cluster,
                    col                     = col.scale);
                    # col                     = colorRamp2(c(min, 0, max), rev(colorRampPalette(brewer.pal(7, "RdYlBu"))(3))));
  
  
  
  ht_list = ht.rmt + ht.sbp
  ht_list = draw(ht_list, main_heatmap = "zscore") # main_heatmap - dictates row order of other heatmaps based on declared heatmap

}

nf.dir      <- "/Users/alderc/1-projects/CAMP/1-AS_timecourse/2-EKD11/1-Pipeline/1-Nextflow/"
work.dir    <- "/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/"
r.dir       <- paste(work.dir, "1-R/",sep='') #change once refactor is completed
tmp.dir     <- paste(work.dir,"tmp/",sep='')
data.dir    <- paste(work.dir,"4-data/",sep='')
results.dir <- paste(work.dir,"2-results_2020/",sep='')
genome.dir  <- "/Users/alderc/1-projects/9-Data/1-Reference_genomes/1-Mus_musculus/"



## Files
gtf.file <- paste(genome.dir, "Mus_musculus.GRCm38.86.gtf", sep="")
design.file <- paste(nf.dir,"design_sbp_rtmt.csv",sep='')
design.baseline <- "naive"

## Output names
count.output <- "sbp_rmt.normalised_counts_with_controls.xls.gz" 
log.output   <-  "sbp_rmt.vst_counts_with_controls.xls.gz"
matrix.output <- "EKD11_sbp_rmt_counts.csv"

## Loading GTF file
gtf.dat    <- import(gtf.file)
gtf.dat    <- as.data.frame(gtf.dat[gtf.dat$type%in% "gene",])f
gtf.dat$GRCm38 <- paste("chr",gtf.dat$seqnames,":",gtf.dat$start,"-",gtf.dat$end,sep='')
rownames(gtf.dat) <- gtf.dat$gene_id

## Contam File 
contam.df <- read.table(paste(data.dir, "Pancreas_contam.txt", sep = ""), col.names = "gene_name")
contam.df$gene_id <- sapply(contam.df$gene_name, function(x){gtf.dat[gtf.dat$gene_name %in% x, "gene_id"]})

## Design file
design           <- read.delim(design.file, header=TRUE, sep=",", stringsAsFactors = TRUE)
rownames(design) <- design$label
design$delivery  <- relevel(design$delivery, design.baseline) # sets design.baseline as first factor (baseline)
design$grp       <- relevel(design$grp, design.baseline)# same as above but for analysis with controls ]
design$day       <- factor(as.character(design$day))
design$replicate <- as.factor(design$replicate) # factors replicates (order not important)



## Importing counts
txi <- tximport(paste(nf.dir, "results/rsem/", design$lims.name, ".genes.results", sep=""), type="rsem")
txi$length[txi$length==0] <- 1

## DESeq2
if (file.exists(paste(r.dir, "Paper/deseq.Rdata", sep = ""))){
  load(paste(r.dir, "Paper/deseq.Rdata", sep = ""))
}else{
  mm <- model.matrix(~delivery + day + delivery:day, design)
  all.zero <- apply(mm, 2, function(x) all(x==0))
  idx <- which(all.zero)
  # delete the column with all zeros from the matrix
  mm = mm[, -idx]
  dat <- DESeqDataSetFromTximport(txi, design, design = ~delivery + day) # For baseline counts - full rank in other file
  rowData(dat) <- gtf.dat[rownames(dat),c("gene_id","gene_name","gene_source","gene_biotype","GRCm38")]
  dds <- DESeq(dat, full = mm)
  vsd <- vst(dat, blind = TRUE)
  vsd_filter <- vst(dds_filter, blind = TRUE)
  save(dds, vsd, file = paste(r.dir, "Paper/deseq.Rdata", sep = ""))}

## Filtering and cleaning
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 4
dds_filter <- dds[filter,]
dds_filter <- dds_filter[-which(rownames(dds_filter) %in% contam.df$gene_id), ]
dds_filter <-DESeq(dds_filter)
# dds_clean <- dds[which(mcols(dds)$betaConv),]
# dds_clean <-DESeq(dds_clean)

vsd_filter <- vst(dds_filter, blind = TRUE)



### Test
dds <- estimateSizeFactors(dat)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, maxit=500)
vsd_test <- vst(dds, blind = F)


## Outputs
# Counts Dataframe (Normalised)
count.norm  <- counts(dds,normalized=TRUE)
count.dBase <- data.frame(count.norm,rowData(dds),stringsAsFactors=FALSE)
rownames(count.dBase) <- count.dBase$gene_id
if(!file.exists(paste(results.dir,count.output ,sep=''))){
  write.table(count.dBase,file=gzfile(paste(results.dir,count.output ,sep='')),col.names=TRUE,row.names=FALSE,sep="\t",quote=F)}


# Log counts (vst)
vst.dBase <- assay(vsd)
vst.dBase <- data.frame(vst.dBase,rowData(dds),stringsAsFactors=FALSE)
if(!file.exists(paste(results.dir,log.output ,sep=''))){
  write.table(vst.dBase,file=gzfile(paste(results.dir,log.output,sep='')),col.names=TRUE,row.names=FALSE,sep="\t",quote=F)}

# Count Matrix (un-normalised)
count.mat <- counts(dds, normalized=FALSE)
if(!file.exists(paste(results.dir, matrix.output ,sep=''))){
  write.table(count.mat, file=gzfile(paste(results.dir, matrix.output, sep='')), col.names=TRUE, row.names=TRUE, sep=',',quote=FALSE)}

### TPM 
tpm.db <- as.data.frame(txi$abundance)
names(tpm.db) <- names(as.data.frame(count.norm))

##################################################################################################################################
##################################################################################################################################

### PCA ###
rv <- rowVars(assay(vsd_filter)) 
select <- order(rv, decreasing = TRUE)[1:1000]

pca <- prcomp(t(assay(vsd_filter)[select, ]))
pca.sum <- summary(pca)
pca.val <- pca.sum$importance[2,1:3]
pca.df <- as.data.frame(pca$x[,1:3])
pca.df <- cbind(pca.df, design[rownames(pca.df), c("delivery", "day", "grp", "label")])
pca.cols <- c("blood" = "#E377C2", "naive" = "black", 'recent_blood' = "#17BECF")


p.1 <- plot_ly(data = pca.df) %>%
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
         # legend = list(legend_title_text = "Day", xanchor = "center", x = 0.925))

p.1
## Exporting File
Sys.setenv("PATH" = paste(Sys.getenv("PATH"), "/anaconda3/bin", sep = .Platform$path.sep)) ## To run orca
orca(p.1, file="../../2-results_2020/PC1_PC2_RMT_SBP_september.pdf", format = "pdf")# Orca will append to WD (no absolute path)


### GGPLOT
pca.plot <- ggplot(pca.df, aes(x=PC1, y=PC2, colour=delivery, shape=day)) + geom_point(size =3);
pca.plot <- pca.plot + xlab(paste0("PC1: ", perc.var[1], "% variance"));
pca.plot <- pca.plot + ylab(paste0("PC2: ", perc.var[2], "% variance"));
pca.plot <- pca.plot + geom_text(aes(label=label.short),hjust=0, vjust=0, size=3, nudge_x = 0.2, nudge_y = 0.4);
pca.cols <- c("mosquito" = "firebrick", "blood" = "darkblue", "naive" = "#E69F00", 'recent_blood' = 'springgreen')
pca.plot <- pca.plot + scale_colour_manual(values = pca.cols)

##################################################################################################################################
##################################################################################################################################

### Pairwise comparison ###

if(file.exists(paste(r.dir, "Paper/results.Rdata", sep = ""))){
  load(paste(r.dir, "Paper/results.Rdata", sep = ""))
}else{
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
      abs(res$log2FC) > 2 &
      res$baseMean>=30       &
      !is.na(res$log2FC)       &
      !is.na(res$FDR)
    res
  })
  
  save(res.list, file = paste(r.dir, "Paper/results.Rdata", sep = "")) ##### SORT LATER 
}

# res.outlier <- list()
# res.outlier[['rtmt.1_vs_sbp.1']]      <- results(dds.outlier, contrast=c("grp","rtmt.1","sbp.1"));
# res.outlier[['rtmt.2_vs_sbp.2']]      <- results(dds.outlier, contrast=c("grp","rtmt.2","sbp.2"));
# res.outlier[['rtmt.3_vs_sbp.3']]      <- results(dds.outlier, contrast=c("grp","rtmt.3","sbp.3"));
# res.outlier[['rtmt.4_vs_sbp.4']]      <- results(dds.outlier, contrast=c("grp","rtmt.4","sbp.4"));
# res.outlier[['rtmt.6_vs_sbp.6']]      <- results(dds.outlier, contrast=c("grp","rtmt.6","sbp.6"));
# 
# res.filter <- list()
# res.filter[['rtmt.1_vs_sbp.1']]      <- results(dds_filter, contrast=c("grp","rtmt.1","sbp.1"));
# res.filter[['rtmt.2_vs_sbp.2']]      <- results(dds_filter, contrast=c("grp","rtmt.2","sbp.2"));
# res.filter[['rtmt.3_vs_sbp.3']]      <- results(dds_filter, contrast=c("grp","rtmt.3","sbp.3"));
# res.filter[['rtmt.4_vs_sbp.4']]      <- results(dds_filter, contrast=c("grp","rtmt.4","sbp.4"));
# res.filter[['rtmt.6_vs_sbp.6']]      <- results(dds_filter, contrast=c("grp","rtmt.6","sbp.6"));
# 
# res.cooks <- list()
# res.cooks[['rtmt.1_vs_sbp.1']]      <- results(dds_filter, contrast=c("grp","rtmt.1","sbp.1"));
# res.cooks[['rtmt.3_vs_sbp.3']]      <- results(dds_filter, contrast=c("grp","rtmt.3","sbp.3"));
# res.cooks[['rtmt.2_vs_sbp.2']]      <- results(dds_filter, contrast=c("grp","rtmt.2","sbp.2"));
# res.cooks[['rtmt.4_vs_sbp.4']]      <- results(dds_filter, contrast=c("grp","rtmt.4","sbp.4"));
# res.cooks[['rtmt.6_vs_sbp.6']]      <- results(dds_filter, contrast=c("grp","rtmt.6","sbp.6"));


deg.list <- lapply(res.list, function(df){
  df <- df[df$DEG, ]
  df$gene_name <- gtf.dat[rownames(df), "gene_name"]
  df
})

cbind(sapply(res.list, function(x){sum(x$DEG)}))

for (n in names(deg.list)){
  tmp.df <- deg.list[[n]]
  write.table(tmp.df, paste(results.dir, n, '.csv', sep = ''),
              row.names = TRUE,
              col.names = NA,
              sep = ',',
              quote = F)
}


rvn.deg <- deg.list[1:5]
svn.deg <- deg.list[6:10]
rvs.deg <- deg.list[11:15]

rvn.genes <- unique(unlist(sapply(rvn.deg, row.names)))
svn.genes <- unique(unlist(sapply(svn.deg, row.names)))
rvs.genes <- unique(unlist(sapply(rvs.deg, row.names)))


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

pdf(file = paste(results.dir, "Paper/Venn_RMT_SBP_September_FC.pdf", sep = ""))
grid.newpage()
grid.draw(venn)
dev.off()


#### Individual Days
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
  
  pdf(file = paste(results.dir, "Paper/Venn/RMT_SBP_Day_", x, "_FC_Scaled_September.pdf", sep = ""))
  grid.newpage()
  grid.draw(venn)
  dev.off()
}


for(i in 1:5){
  s <- rownames(svn.deg[[i]])
  r <- rownames(rvn.deg[[i]])
  ifelse(i == 5, x <- 6, x <- i)
  r.unique <- setdiff(r,s)
  for (g in r.unique){
    gene_name <- gtf.dat[g, "gene_name"]
    pdf(paste(results.dir, "RMT_SBP_gene_profiles/RMT_Unique/Day_", x, "/", gene_name, ".pdf", sep = ""))
    plot(count.plot(dds, g))
    dev.off()
  }
}

for(i in 1:5){
  s <- rownames(svn.deg[[i]])
  r <- rownames(rvn.deg[[i]])
  ifelse(i == 5, x <- 6, x <- i)
  s.unique <- setdiff(s,r)
  for (g in s.unique){
    gene_name <- gtf.dat[g, "gene_name"]
    pdf(paste(results.dir, "RMT_SBP_gene_profiles/SBP_Unique/Day_", x, "/", gene_name, ".pdf", sep = ""))
    plot(count.plot(dds, g))
    dev.off()
  }
}

# ### Looking into time differences of gene activation
# s.day <- lapply(svn.deg, row.names)
# r.day <- lapply(rvn.deg, row.names)
# 
# r.same <- list()
# r.early <- list()
# r.delayed <- list()
# r.unique <- list()
# max.day <- 5
# for (n in 1:max.day){
#   r = r.day[[n]]
#   for (g in r){
#     for (m in 1:max.day){
#       if (g %in% s.day[[m]]){
#         if (m == n){
#           r.same <- unlist(c(r.same, g))
#         } else if (m > n){
#           r.delayed <- unlist(c(r.delayed, g))
#         } else if (m < n){
#           r.early <- unlist(c(r.early, g))
#         }
#       }
#     }
#   }
# } 

#########TESTING
cbind(sapply(deg.list,function(x){sum(x$DEG)}));

r <- unique(c(rownames(deg.list[[4]]), rownames(deg.list[[5]])))
s <- unique(c(rownames(deg.list[[9]]), rownames(deg.list[[10]])))

d4_6.common <- intersect(r,s)
r.common <- deg.list[[4]][d4.common, c("log2FC", "FDR")]
s.common <- deg.list[[9]][d4.common, c("log2FC", "FDR")]
sig <- which(abs(r.common$log2FC - s.common$log2FC)  > 0.3)
d4_6.sig <- d4_6.common[sig]


### Looking at DEGs and there assosciated difference in FC from rvn to svn
rvn.res <- res.list[1:5]
svn.res <- res.list[6:10]
s.day <- lapply(svn.deg, row.names)
r.day <- lapply(rvn.deg, row.names)

test.list <- list()
for (i in 1:5){
deg.day <- unique(c(r.day[[i]], s.day[[i]]))
r.res <- rvn.res[[i]][deg.day, c("log2FC", "FDR")]
s.res <- svn.res[[i]][deg.day, c("log2FC", "FDR")]
sig <- which(abs(r.res$log2FC - s.res$log2FC)  > 0.5)
test.list[[i]] <- deg.day[sig]
print(gtf.dat[deg.day[sig], "gene_name"])
}

test.genes <- unique(unlist(test.list))

dir1 = paste(results.dir, "Paper/", sep = "")
dir.create(dir1)
for (n in 1:5){
  ifelse(n == 5, day <- "Day_6/", day <-  paste("Day_", n, "/", sep = ""))
  g.list = test.list[[n]]
  dir2 = paste(dir1, day, sep = "")
  dir.create(dir2)
  for (g in g.list){
    name = gtf.dat[g, "gene_name"]
    pdf(paste(dir2, name, ".pdf", sep = ""))
    tpm.plot(g)
    dev.off()
    }
}




## Venn Diagram with test gene list

rvn.fc <- intersect(rvn.genes, test.genes)
svn.fc <- intersect(svn.genes, test.genes)
venn.cols <- c("#E377C2", "#17BECF")
venn <- venn.diagram(x = list(SBP = svn.fc, RMT = rvn.fc),
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

pdf(file = paste(results.dir, "August/Venn_RMT_SBP_FC_August.pdf", sep = ""))
grid.newpage()
grid.draw(venn)
dev.off()

rvn.fc.unique <- setdiff(rvn.fc, svn.fc)
svn.fc.unique <- setdiff(svn.fc, rvn.fc)
gtf.dat[svn.fc.unique, "gene_name"]

for (g in rvn.fc.unique){
  gene_name <- gtf.dat[g, "gene_name"]
  pdf(paste(results.dir, "August/rvn_fc_unqiue_TPM/", gene_name,"_TPM.pdf", sep =""))
  tpm.plot(g)
  dev.off()
}

for (g in svn.fc.unique){
  gene_name <- gtf.dat[g, "gene_name"]
  pdf(paste(results.dir, "August/svn_fc_unqiue_TPM/", gene_name,"_TPM.pdf", sep =""))
  tpm.plot(g)
  dev.off()
}
#### Cytokine Registry
cytokine.db <- xlsx::read.xlsx2(file = paste(data.dir, "CytokineRegistry.November_2015.xls", sep = ""), sheetIndex = 1 )
cytokine.sym <- cytokine.db[,12] %>% unique()
cytokine.sym <- cytokine.sym[-which(cytokine.sym == "")]
cytokine.sym <- gsub(x = cytokine.sym, pattern = " ", replacement = "")
cytokine.id <- gtf.dat[gtf.dat$gene_name %in% cytokine.sym,  "gene_id"]

rvn.cyto <- intersect(rvn.genes,cytokine.id)
svn.cyto <- intersect(svn.genes, cytokine.id)
rvs.cyto <- intersect(rvs.genes, cytokine.id)

cyto.sig <-c(setdiff(svn.cyto, rvn.cyto), setdiff(rvn.cyto, svn.cyto))



for (g in c(cyto.sig, rvs.cyto)){
  gene_name <- gtf.dat[g, "gene_name"]
  pdf(paste(results.dir, "August/", gene_name,"_TPM.pdf", sep =""))
  tpm.plot(g)
  dev.off()
}
setdiff(svn.cyto, rvn.cyto)

rvs.genes

cyto.id


## IG genes
ig.genes <- gtf.dat[grep("Ig", gtf.dat$gene_name), "gene_id"]

all.deg <- unique(unlist(sapply(deg.list, row.names)))

ig.deg <- intersect(ig.genes, all.deg)

count.plot(dds, ig.deg[1])

#### Heatmap
heat.genes <- svn.fc.unique
heat.dat <- assay(vsd)[heat.genes, ]
heat.dat   <- t(apply(heat.dat,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))})) %>% data.frame() %>% drop_na() %>% as.matrix()
rownames(heat.dat) <- gtf.dat[rownames(heat.dat),"gene_name"]


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
max <- max(heat.dat)
min <- min(heat.dat)
n <- ncol(heat.dat)
ht <- Heatmap(heat.dat,
              name = "zscore",
              show_row_names          = TRUE,
              show_column_names       = TRUE,
              column_names_side       = "bottom",
              show_column_dend        = FALSE,
              cluster_columns         = FALSE,
              cluster_rows            = TRUE, 
              rect_gp                 = gpar(col = "white", lwd = 1),
              row_dend_width          = unit(30, "mm"),
              column_dend_height      = unit(30, "mm"),
              show_heatmap_legend     = TRUE,
              clustering_distance_rows = 'euclidean',
              clustering_method_rows = "complete",
              column_names_max_height = unit(8, "cm"),
              row_names_gp            = gpar(fontsize = 8),
              column_names_gp         = gpar(fontsize = 12),
              col                     = colorRamp2(c(min, 0, max), rev(colorRampPalette(brewer.pal(9, "RdBu"))(3))));

pdf(paste(results.dir, "August/svn_fc_unique_heatmap.pdf", sep=""), height = 25, width = 15)
ht
decorate_heatmap_body("zscore", c(for(d in seq(0, n, 4)){
  grid.lines(x = c(d/n, d/n), y = c(0,1), gp= gpar(col = "black", lty="solid", lwd=2))},
  grid.lines(x = c(0,1), y = c(0,0), gp= gpar(col = "black", lty="solid", lwd=2)),
  grid.lines(x = c(0,1), y = c(1,1), gp= gpar(col = "black", lty="solid", lwd=2))))
dev.off()


######
test <- res.list[[6]]
hist(test$pvalue, breaks = seq(0,1,0.01))
#################################################################################



common.genes <- intersect(rvn.genes, svn.genes)

rvn.unique <- setdiff(rvn.genes, svn.genes)
rvn.pc <- intersect(rvn.unique, gtf.dat$gene_id[gtf.dat$gene_biotype %in% "protein_coding"])

rvn.pc.df <- data.frame(gene_id = rvn.pc, gene_name = gtf.dat[rvn.pc, "gene_name"])

write.table(rvn.pc.df, file = paste(results.dir, "rvn_pc_unique_df.tsv", sep = ""), sep = "\t", row.names = F, quote = F)

rvn.deg.upreg <- unique(unlist(sapply(rvn.deg, function(df){
  rownames(df[df$log2FC > 0, ])
})))
rvn.deg.downreg <- unique(unlist(sapply(rvn.deg, function(df){
  rownames(df[df$log2FC < 0, ])
})))

svn.deg.upreg <- unique(unlist(sapply(svn.deg, function(df){
  rownames(df[df$log2FC > 0, ])
})))
svn.deg.downreg <- unique(unlist(sapply(svn.deg, function(df){
  rownames(df[df$log2FC < 0, ])
})))

common.upreg <- intersect(rvn.deg.upreg, svn.deg.upreg)
common.downreg <- intersect(rvn.deg.downreg, svn.deg.downreg)

rvn.upreg.unique <- setdiff(rvn.deg.upreg, svn.deg.upreg)
rvn.downreg.unique <- setdiff(rvn.deg.downreg, svn.deg.downreg)

rvn.df <- data.frame(gene_id = c(rvn.upreg.unique, rvn.downreg.unique))
rvn.df$up_down <- sapply(rvn.df$gene_id, function(x){ifelse(x %in% rvn.upreg.unique, "Up", "Down")})
rvn.df$gene_name <- gtf.dat[rvn.df$gene_id, "gene_name"]
rvn.df.pc <- rvn.df[rvn.df$gene_id %in% gtf.dat$gene_id[gtf.dat$gene_biotype %in% "protein_coding"],  ]


count.plot(dds, "ENSMUSG00000032754")



rvn.d4 <- rvn.deg[[4]]
svn.d4 <- svn.deg[[4]]

r_unique.d4 <- setdiff(rownames(rvn.d4), rownames(svn.d4))


r_unique.d4.df <- rvn.d4[r_unique.d4, ]
r_unique.d4.df[order(r_unique.d4.df$log2FC, decreasing = T), ]



##### Day 4 and 6 only
cbind(sapply(res.list, function(x){sum(x$DEG)}))

rvn.deg.late <- rvn.deg[4:5]
svn.deg.late <- svn.deg[4:5]

rvn.deg.late.2 <- lapply(rvn.deg.late, function(deg){
  deg %>% subset(abs(log2FC) > 2)
})

svn.deg.late.2 <- lapply(svn.deg.late, function(deg){
  deg %>% subset(abs(log2FC) > 2)
})

rvn.late.2 <- unique(unlist(sapply(rvn.deg.late.2, row.names)))
svn.late.2 <- unique(unlist(sapply(svn.deg.late.2, row.names)))

deg.late.2 <- unique(c(rvn.late.2, svn.late.2))
deg.late.2.sym <- gtf.dat[deg.late.2, "gene_name"]
deg.late.2.sym <- deg.late.2.sym[-grep("Igh|Igk", deg.late.2.sym)]
deg.late.id <- unlist(sapply(deg.late.2.sym, function(x){gtf.dat[gtf.dat$gene_name %in% x, "gene_id"]}))


ig.rmv <- c("Ighv3-6", "Igkv5-43", "Lmo7", "Gm12504", "Pde2a", "Ighv2-3", "Ighv1-4", "Igkv1-99")
ig.rmv <- unlist(sapply(ig.rmv, function(x){gtf.dat[gtf.dat$gene_name %in% x, "gene_id"]}))

pdf("./Paper/d4_6_DEG_sample_clusters_HM.pdf", height = 20, width = 10)
split.heatmap(gene_list = deg.late.id, heat.dat = vsd.late, z.both = TRUE, cluster_samples = TRUE)
dev.off()

pdf("./Paper/d4_6_DEG_sample_clusters_non_split_HM.pdf", height = 20, width = 10)
heatmap.draw(gene_list = deg.late.id, heat.dat = vsd.late, cluster_samples = TRUE)
dev.off()

####Cytokines
pdf("./Paper/d4_6_cytokines_HM.pdf", height = 10, width = 10)
split.heatmap(gene_list = cytokine.id, heat.dat = vsd.late)
dev.off()


## Module_genes

module.genes <- read.table("../4-data/blood_module_gene_list.tsv", sep = "\t", header = T) %>% column_to_rownames("module")
 
str_split(module.genes[1,2], pattern = ", ")[[1]]

module.list <- lapply(rownames(module.genes), function(module){
  module.string <- module.genes[module, 1]
  module.sym <- str_split(module.string, pattern = ", ")[[1]]
  module.id <- unlist(sapply(module.sym, function(x){gtf.dat[gtf.dat$gene_name %in% x, "gene_id"]}))
  return(module.id)
})
names(module.list) <- rownames(module.genes)


pdf("./Paper/modules_b31_HM.pdf", height = 15, width = 10)
split.heatmap(module.list[[31]], heat.dat = vsd.late)
dev.off()

## Innate DB 

innate.sym <- readxl::read_xls("../4-data/innatedb_curated_genes.xls") %>% 
  subset(Species == 10090) %>% select("Gene Symbol") %>% unique()

innate.id <- unlist(sapply(innate.sym, function(x){gtf.dat[gtf.dat$gene_name %in% x, "gene_id"]}))[,1]

innate_deg.id <- intersect(innate.id, deg.late.id)

pdf("./Paper/innate_db_HM.pdf", height = 20, width = 10)
split.heatmap(innate.id, heat.dat = vsd.late)
dev.off()
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################



vsd.dat <- assay(vsd)
vsd.late <- vsd.dat[, grep("naive|d4|d6", colnames(vsd.dat))]

split.heatmap(d4_6.sig, heat.dat = vsd.late)


##########LRT
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
filter <- rowSums(nc >= 10) >= 2
dds.rvs <- dds.rvs[filter,]
dds.rvs <-DESeq(dds.rvs)

vsd.rvs <- vst(dat.rvs, blind = TRUE)

rvs.lrt <- DESeq(dds.rvs, test="LRT", reduced = ~ delivery + day) ## Showing genes with any different across time

results(dds.rvs,)

lrt.res <- as.data.frame(results(rvs.lrt, alpha = 0.1))
lrt.sig <- lrt.res[lrt.res$padj < 0.1, ] %>% drop_na()
lrt.sig$gene_name <- gtf.dat[rownames(lrt.sig), "gene_name"]
lrt.sig$gene_name
lrt.sig

dat.rvs2 <- DESeqDataSetFromTximport(txi.rvs, design.rvs, ~day + delivery + day:delivery)
dds.rvs2 <- DESeq(dat.rvs2)

rvs.lrt2 <- DESeq(dds.rvs2, test="LRT", reduced = ~ day)

lrt.res2 <- as.data.frame(results(rvs.lrt2, alpha = 0.1))
lrt.sig2 <- lrt.res2[lrt.res2$padj < 0.1, ] %>% drop_na()

dir.create(paste(dir1, "LRT/", sep = ""))

# for (g in rownames(lrt.sig)){
#   gene_name <- gtf.dat[g, "gene_name"]
#   pdf(paste(results.dir, "RMT_SBP_LRT_gene_profiles/", gene_name, ".pdf", sep = ""))
#   plot(count.plot(dds,g))
#   dev.off()
# }

for(g in rownames(lrt.sig)){
  gene_name <- lrt.sig[g, "gene_name"]
  pdf(paste(results.dir, "August/LRT_TPM/", gene_name, "_TPM.pdf", sep = ""))
  tpm.plot(g)
  dev.off()
}

design(dds.rvs) <- ~grp
dds.rvs <- DESeq(dds.rvs)


design.rvs$grp
rvs.list <- list()
rvs.list[['day1']]      <- results(dds.rvs, contrast=c("grp","rtmt.1","sbp.1"));
rvs.list[['day2']]      <- results(dds.rvs, contrast=c("grp","rtmt.2","sbp.2"));
rvs.list[['day3']]      <- results(dds.rvs, contrast=c("grp","rtmt.3","sbp.3"));
rvs.list[['day4']]      <- results(dds.rvs, contrast=c("grp","rtmt.4","sbp.4"));
rvs.list[['day6']]      <- results(dds.rvs, contrast=c("grp","rtmt.6","sbp.6"));


rvs.list <- lapply(rvs.list,function(res){
  res <- as.data.frame(res);
  colnames(res)[colnames(res) %in% "log2FoldChange"] <- "log2FC";
  colnames(res)[colnames(res) %in% "padj"] <- "FDR";
  res$DEG <- res$FDR <= 0.1 &
    # abs(res$log2FC)>=2 & #only filters out logFC ±1
    res$baseMean>=30         &
    !is.na(res$log2FC)       &
    !is.na(res$FDR);
  res;
})
cbind(sapply(rvs.list,function(x){sum(x$DEG)}));

d3 <- rvs.list[[3]]
d3[d3$DEG,]
design(dds.rvs)

for (g in rownames(rvs.d5)){ #should be d6
  gene_name <- gtf.dat[g, "gene_name"]
  pdf(paste(results.dir, "RMT_SBP_gene_profiles/RMT_vs_SBP/Day6/", gene_name, ".pdf", sep = ""))
  plot(count.plot(dds, g))
  dev.off()
}





library(DEGreport)
lrt.rlog <- assay(vsd.rvs[rownames(lrt.sig), ])
nrow(lrt.rlog)

clusters <- degPatterns(lrt.rlog, meta = design.rvs, time = "day", col="delivery")

cluster_groups <- clusters$df
cluster_groups$gene_name <- gtf.dat[rownames(cluster_groups), "gene_name"]
cluster_groups

cluster.plot <- degPlotCluster(clusters$normalized, time = "day", color = "delivery")
cluster.plot <- cluster.plot + scale_colour_manual(values = cnt.cols)

##### RVS DEG Heatmaps - d4, 6 and control
rvs.genes <- unique(unlist(sapply(rvs.deg, row.names)))
rvs.remove <- c("ENSMUSG00000023433", "ENSMUSG00000031957", "ENSMUSG00000093931", "ENSMUSG00000096569", "ENSMUSG00000096770", "ENSMUSG00000024225",
                "ENSMUSG00000036938", "ENSMUSG00000042179", "ENSMUSG00000046008", "ENSMUSG00000054106", "ENSMUSG00000054446", "ENSMUSG00000074268")
rvs.genes <- setdiff(rvs.genes, rvs.remove)
heat.genes <- rvs.genes
heat.dat <-as.matrix(assay(vsd)[heat.genes, ])
# sample.rmv <- grep("d1|d2|d3", colnames(heat.dat))
# heat.dat <- heat.dat[, -sample.rmv]
heat.dat   <- t(apply(heat.dat,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))}))
rownames(heat.dat) <- gtf.dat[heat.genes,"gene_name"]

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
max <- max(heat.dat)
min <- min(heat.dat)
n <- ncol(heat.dat)
ht <- Heatmap(heat.dat,
              name = "zscore",
              show_row_names          = TRUE,
              show_column_names       = TRUE,
              column_names_side       = "bottom",
              show_column_dend        = FALSE,
              cluster_columns         = FALSE,
              rect_gp                 = gpar(col = "white", lwd = 1),
              row_dend_width          = unit(30, "mm"),
              column_dend_height      = unit(30, "mm"),
              show_heatmap_legend     = TRUE,
              clustering_distance_rows = 'pearson',
              column_names_max_height = unit(8, "cm"),
              row_names_gp            = gpar(fontsize = 6),
              column_names_gp         = gpar(fontsize = 10),
              col                     = colorRamp2(c(min, 0, max), rev(colorRampPalette(brewer.pal(9, "RdBu"))(3))));

pdf(paste(results.dir, "RMTvsSBP_DEG_HM_July.pdf", sep=""), height = 20, width = 15)
ht
decorate_heatmap_body("zscore", c(for(d in seq(0, n, 4)){
  grid.lines(x = c(d/n, d/n), y = c(0,1), gp= gpar(col = "black", lty="solid", lwd=2))},
  grid.lines(x = c(0,1), y = c(0,0), gp= gpar(col = "black", lty="solid", lwd=2)),
  grid.lines(x = c(0,1), y = c(1,1), gp= gpar(col = "black", lty="solid", lwd=2))))
dev.off()

tpm.plot("ENSMUSG00000042066")

count.plot(dds, "ENSMUSG00000046008")
