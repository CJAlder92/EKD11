rm(list = ls())
source('/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/1-R/refactor/packages.R')

nf.dir      <- "/Users/alderc/1-projects/CAMP/1-AS_timecourse/2-EKD11/1-Pipeline/1-Nextflow/"
work.dir    <- "/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/"
r.dir       <- paste(work.dir, "1-R/",sep='') #change once refactor is completed
tmp.dir     <- paste(work.dir,"tmp/",sep='')
data.dir    <- paste(work.dir,"4-data/",sep='')
results.dir <- paste(work.dir,"2-results_2020/",sep='')
genome.dir  <- "/Users/alderc/1-projects/9-Data/1-Reference_genomes/1-Mus_musculus/"

## Files
gtf.file <- paste(genome.dir, "Mus_musculus.GRCm38.86.gtf", sep="")
design.file <- paste(nf.dir,"design.csv",sep='')
design.baseline <- "naive"

## Output names
count.output <- "genes.normalised_counts_with_controls.xls.gz" 
log.output   <-  "genes.vst_counts_with_controls.xls.gz"
matrix.output <- "EKD11_counts.csv"

## Loading GTF file
gtf.dat    <- import(gtf.file)
gtf.dat    <- as.data.frame(gtf.dat[gtf.dat$type%in% "gene",])
gtf.dat$GRCm38 <- paste("chr",gtf.dat$seqnames,":",gtf.dat$start,"-",gtf.dat$end,sep='')
rownames(gtf.dat) <- gtf.dat$gene_id

## Design file
design.file      <- paste(nf.dir, "design_2020.csv", sep='') # For 2020 analysis
design           <- read.delim(design.file, header=TRUE, sep=",", stringsAsFactors = TRUE)
rownames(design) <- design$label
design$delivery  <- relevel(design$delivery, design.baseline) # sets design.baseline as first factor (baseline)
design$grp       <- relevel(design$grp, design.baseline) # same as above but for analysis with controls 
design$day       <- factor(as.character(design$day), levels=as.character(sort(unique(design$day)))) # factors days in numberical order
design$replicate <- as.factor(design$replicate) # factors replicates (order not important)

## Importing counts
txi <- tximport(paste(nf.dir, "results/rsem/", design$lims.name, ".genes.results", sep=""), type="rsem")
txi$length[txi$length==0] <- 1

## DESeq2
if (file.exists(paste(r.dir, "Projects/deseq.Rdata", sep = ""))){
  load(paste(r.dir, "Projects/deseq.Rdata", sep = ""))
}else{
  dat <- DESeqDataSetFromTximport(txi, design, ~ grp) # For baseline counts - full rank in other file
  rowData(dat) <- gtf.dat[rownames(dat),c("gene_id","gene_name","gene_source","gene_biotype","GRCm38")]
  dds <- DESeq(dat)
  vsd <- vst(dat, blind = TRUE)
  save(dds, vsd, file = paste(r.dir, "Projects/deseq.Rdata", sep = ""))}

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


##################################################################################################################################
##################################################################################################################################

### PCA ###
rv <- rowVars(assay(vsd)) 
select <- order(rv, decreasing = TRUE)[1:1000]

pca <- prcomp(t(assay(vsd)[select, ]))
pca.sum <- summary(pca)
pca.val <- pca.sum$importance[2,1:3]
pca.df <- as.data.frame(pca$x[,1:3])
pca.df <- cbind(pca.df, design[rownames(pca.df), c("delivery", "day", "grp", "label")])
pca.cols <- c("mosquito" = "#BCBD22", "blood" = "#E377C2", "naive" = "black", 'recent_blood' = "#17BECF")


p.1 <- plot_ly(data = pca.df) %>%
  add_trace(type = "scatter",
            x = ~PC1, y = ~PC2,
            color = ~delivery, colors = pca.cols,
            symbols = c("0" = "triangle-up" ,"1" = "circle", "2" = "cross", "3" = "x", "4" = "square", "6" = "diamond"),
            text = ~label,
            mode = "markers") %>% 
  add_trace(type = "scatter",
            x = ~PC1, y = ~PC2,
            symbol = ~day,
            mode = "markers",
            text = ~label,
            marker = list(color = "grey", size = 8)) %>% 
  add_trace(type = "scatter",
            x = ~PC1, y = ~PC2,
            color = ~delivery,
            symbol = ~day,
            mode = "markers",
            text = ~label,
            marker = list(size = 10),
            showlegend = F) %>% 
  layout(xaxis = list(title = paste("PC1:", round(pca.val[[1]], 2)*100, "%", sep = " ")),
         yaxis = list(title = paste("PC2:", round(pca.val[[2]], 2)*100, "%", sep = " ")))

p.1

##################################################################################################################################
##################################################################################################################################

### Pairwise comparison ###

if(file.exists(paste(r.dir, "Projects/results.Rdata", sep = ""))){
  load(paste(r.dir, "Projects/results.Rdata", sep = ""))
}else{
  res.list <- list()
  
  res.list[["mt.1_vs_naive"]]      <- results(dds, contrast=c("grp","mt.1","naive"))
  res.list[["mt.2_vs_naive"]]      <- results(dds, contrast=c("grp","mt.2","naive"))
  res.list[["mt.3_vs_naive"]]      <- results(dds, contrast=c("grp","mt.3","naive"))
  res.list[["mt.4_vs_naive"]]      <- results(dds, contrast=c("grp","mt.4","naive"))
  res.list[["mt.6_vs_naive"]]      <- results(dds, contrast=c("grp","mt.6","naive"))
  
  res.list[["rtmt.1_vs_naive"]]      <- results(dds, contrast=c("grp","rtmt.1","naive"))
  res.list[["rtmt.2_vs_naive"]]      <- results(dds, contrast=c("grp","rtmt.2","naive"))
  res.list[["rtmt.3_vs_naive"]]      <- results(dds, contrast=c("grp","rtmt.3","naive"))
  res.list[["rtmt.4_vs_naive"]]      <- results(dds, contrast=c("grp","rtmt.4","naive"))
  res.list[["rtmt.6_vs_naive"]]      <- results(dds, contrast=c("grp","rtmt.6","naive"))
  
  res.list[["sbp.1_vs_naive"]]      <- results(dds, contrast=c("grp","sbp.1","naive"))
  res.list[["sbp.2_vs_naive"]]      <- results(dds, contrast=c("grp","sbp.2","naive"))
  res.list[["sbp.3_vs_naive"]]      <- results(dds, contrast=c("grp","sbp.3","naive"))
  res.list[["sbp.4_vs_naive"]]      <- results(dds, contrast=c("grp","sbp.4","naive"))
  res.list[["sbp.6_vs_naive"]]      <- results(dds, contrast=c("grp","sbp.6","naive"))
  
  res.list[['mt.1_vs_sbp.1']]      <- results(dds, contrast=c("grp","mt.1","sbp.1"));
  res.list[['mt.2_vs_sbp.2']]      <- results(dds, contrast=c("grp","mt.2","sbp.2"));
  res.list[['mt.3_vs_sbp.3']]      <- results(dds, contrast=c("grp","mt.3","sbp.3"));
  res.list[['mt.4_vs_sbp.4']]      <- results(dds, contrast=c("grp","mt.4","sbp.4"));
  res.list[['mt.6_vs_sbp.6']]      <- results(dds, contrast=c("grp","mt.6","sbp.6"));
  
  res.list[['mt.1_vs_rtmt.1']]      <- results(dds, contrast=c("grp","mt.1","rtmt.1"));
  res.list[['mt.2_vs_rtmt.2']]      <- results(dds, contrast=c("grp","mt.2","rtmt.2"));
  res.list[['mt.3_vs_rtmt.3']]      <- results(dds, contrast=c("grp","mt.3","rtmt.3"));
  res.list[['mt.4_vs_rtmt.4']]      <- results(dds, contrast=c("grp","mt.4","rtmt.4"));
  res.list[['mt.6_vs_rtmt.6']]      <- results(dds, contrast=c("grp","mt.6","rtmt.6"));
  
  res.list[['rtmt.1_vs_sbp.1']]      <- results(dds, contrast=c("grp","rtmt.1","sbp.1"));
  res.list[['rtmt.2_vs_sbp.2']]      <- results(dds, contrast=c("grp","rtmt.2","sbp.2"));
  res.list[['rtmt.3_vs_sbp.3']]      <- results(dds, contrast=c("grp","rtmt.3","sbp.3"));
  res.list[['rtmt.4_vs_sbp.4']]      <- results(dds, contrast=c("grp","rtmt.4","sbp.4"));
  res.list[['rtmt.6_vs_sbp.6']]      <- results(dds, contrast=c("grp","rtmt.6","sbp.6"));
  
  save(res.list, file = paste(r.dir, "Projects/results.Rdata", sep = "")) ##### SORT LATER 
}

res.list <- lapply(res.list,function(res){
  res <- as.data.frame(res)
  colnames(res)[colnames(res) %in% "log2FoldChange"] <- "log2FC"
  colnames(res)[colnames(res) %in% "padj"] <- "FDR"
  res$DEG <- res$FDR <= 0.05 &
    abs(res$log2FC)>= 1 &  #only filters out logFC ±1 N.B CHANGE FOR RMT vs SBP COMPARISONS
    res$baseMean>=30         &
    !is.na(res$log2FC)       &
    !is.na(res$FDR)
  res
})

deg.list <- lapply(res.list, function(df){
  df <- df[df$DEG, ]
})

##################################################################################################################################
##################################################################################################################################







