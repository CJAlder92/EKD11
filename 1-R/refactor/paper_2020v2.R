rm(list = ls())
source('/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/1-R/refactor/packages.R')

nf.dir      <- "/Users/alderc/1-projects/CAMP/1-AS_timecourse/2-EKD11/1-Pipeline/1-Nextflow/"
work.dir    <- "/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/"
r.dir       <- paste(work.dir, "1-R/",sep='') #change once refactor is completed
tmp.dir     <- paste(work.dir,"tmp/",sep='')
data.dir    <- paste(work.dir,"4-data/",sep='')
results.dir <- paste(work.dir,"2-results_2020/",sep='')
genome.dir  <- "/Users/alderc/1-projects/9-Data/1-Reference_genomes/1-Mus_musculus/"

# Checks for directories before creating (to prevent overwriting)
for (dir in grep("\\.dir$",ls(),value=T)) {
  if (!file.exists(get(dir))) { dir.create(get(dir),recursive=TRUE, mode="0755"); }
}

## Files
gtf.file <- paste(genome.dir, "Mus_musculus.GRCm38.86.gtf", sep="")
design.file <- paste(nf.dir,"design.csv",sep='')
design.baseline <- "naive" # This is the treatment group you want to use for baseline measurements and test others against - Usually control

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
# design.file <- paste(nf.dir, "design.csv", sep = "")
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
# save(object = vsd,  file = paste(r.dir, "Projects/deseq.Rdata", sep = ""))}

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
# write.table(vst.dBase,file=gzfile(paste(results.dir,log.output,sep='')),col.names=TRUE,row.names=FALSE,sep="\t",quote=F)

# Count Matrix (un-normalised)
count.mat <- counts(dds, normalized=FALSE)
# write.table(count.mat, file=gzfile(paste(results.dir, matrix.output, sep='')), col.names=TRUE, row.names=TRUE, sep=',',quote=FALSE)

### Put into separate R file at some point
ensembl.decode <- data.frame(gene_id = rownames(gtf.dat),
                             gene_name = gtf.dat$gene_name)

# write.table(ensembl.decode, file = paste(results.dir, "ensembl_decode.csv", sep=""), col.names = TRUE, row.names = TRUE, sep=",", quote=FALSE)


##################################################################################################################################
##################################################################################################################################

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
## Exporting File
# Sys.setenv("PATH" = paste(Sys.getenv("PATH"), "/anaconda3/bin", sep = .Platform$path.sep)) ## To run orca
# orca(p.1, file="../../2-results_2020/PC1_PC2.pdf", format = "pdf")# Orca will append to WD (no absolute path)

p.1

##################################################################################################################################
##################################################################################################################################


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

## Creating DEG table for each condition (compared against naive)
mvn.deg <- deg.list[grep('^mt.._vs_naive', names(deg.list))]
rvn.deg <- deg.list[grep('^rtmt.._vs_naive', names(deg.list))]
svn.deg <- deg.list[grep('^sbp.._vs_naive', names(deg.list))]


## Create DEG gene lists for each condition
mvn.genes <- unique(unlist(sapply(mvn.deg, row.names)))
rvn.genes <- unique(unlist(sapply(rvn.deg, row.names))) 
svn.genes <- unique(unlist(sapply(svn.deg, row.names)))

common.genes <- Reduce(intersect, list(mvn.genes, rvn.genes, svn.genes)) # list of genes found in all conditions

print(paste(length(mvn.genes), length(rvn.genes), length(svn.genes), sep =","))

length(common.genes)

## Venn Diagram 
library("VennDiagram")
venn.cols <- c("#BCBD22", "#E377C2", "#17BECF")
venn <- venn.diagram(x = list(mvn.genes, svn.genes, rvn.genes),
                     category.names = c("MT", "SBP", "RTMT"),
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
                     cat.pos = c(-27, 27, 135),
                     cat.dist = c(0.055, 0.055, 0.085),
                     cat.fontfamily = "sans",
                     cat.cex = 2.5,
                     rotation = 1);

# pdf(file = paste(results.dir, "venn/All_newcols.pdf", sep = ""))
grid.newpage()
grid.draw(venn)
dev.off()

all.genes <- unique(unlist(sapply(deg.list, row.names)))
length(all.genes)

all.pc.genes <- intersect(all.genes,gtf.dat$gene_id[gtf.dat$gene_biotype %in% "protein_coding"])
length(all.pc.genes)

mvn.pc <- intersect(mvn.genes, all.pc.genes)
rvn.pc <- intersect(rvn.genes, all.pc.genes)
svn.pc <- intersect(svn.genes, all.pc.genes)

common.pc <- Reduce(intersect, list(mvn.pc, rvn.pc, svn.pc))

venn <- venn.diagram(x = list(mvn.pc, svn.pc, rvn.pc),
                     category.names = c("MT", "SBP", "RTMT"),
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
                     cat.pos = c(-27, 27, 135),
                     cat.dist = c(0.055, 0.055, 0.085),
                     cat.fontfamily = "sans",
                     cat.cex = 2.5,
                     rotation = 1)
pdf(file = paste(results.dir, "venn/All_newcols_PC.pdf", sep = ""))
grid.newpage()
grid.draw(venn)
dev.off()


##################################################################################################################################
##################################################################################################################################

rt.mt.unique <- intersect(mvn.pc, rvn.pc)
rt.mt.unique <- setdiff(rt.mt.unique, common.pc)

rt.mt.df <- data.frame(gene_id = rt.mt.unique)
rt.mt.df$gene_name <- gtf.dat[rt.mt.df$gene_id, "gene_name"]

write.table(rt.mt.df, file= paste(results.dir, "rt_mt_df.csv", sep = ""), quote = F, row.names = F, sep = ",")

rt.mt.vsd <- as.data.frame(assay(vsd)[rt.mt.ids, ])
rt.mt.vsd$gene_name <- gtf.dat[rownames(rt.mt.vsd), "gene_name"]
write.table(rt.mt.vsd, file= paste(results.dir, "rt_mt_log_counts.csv", sep = ""), quote = F, row.names = F, sep = ",")


########## Gene Profiles 

gene.id <- rt.mt.unique
for (gid in gene.id){
  sym <- gtf.dat[gid, "gene_name"]
  plot.file <- paste(results.dir,"rmt_mt_shared_gene_profiles",sym,".png",sep='');
  if (!file.exists(plot.file)) {
    cnt.dat  <- plotCounts(dds, which(rownames(dds) %in% gid), intgroup = c("day","delivery"), returnData = TRUE);
    cnt.trans <- cnt.dat[1:60,]
    cnt.naive <- cnt.dat[61:64, ]
    naive.mean <- mean(cnt.naive$count)
    #cnt.dat$day <- as.numeric(as.character(cnt.dat$day))
    cnt.plot <- ggplot(cnt.trans, aes(x = day, y = count, color = delivery, fill=delivery, group = delivery)) + geom_smooth(se = FALSE, alpha=0.1, method = "loess", span=0.75);
    cnt.plot <- cnt.plot + geom_hline(yintercept = naive.mean, col = "black", linetype="dashed")
    cnt.plot <- cnt.plot + ggtitle(sym) + ylab("normalised count") + xlab("day")
    cnt.plot <- cnt.plot + scale_x_discrete(limits=1:6)
    #cnt.plot <- cnt.plot + scale_y_log10();
    cnt.cols <- c("mosquito" = "#BCBD22", "blood" = "#E377C2", "naive" = "Black", 'recent_blood' = "#17BECF")
    cnt.plot <- cnt.plot + scale_colour_manual(values = cnt.cols)
    png(file=paste(results.dir, "rmt_mt_shared_gene_profiles/" ,sym,".png",sep=''),height=500,width=700);
    print(cnt.plot);
    dev.off();
  }
}

##################################################################################################################################
# RTMT:MT Venn Day by Day

# Individual days

for (i in 1:5){
  m <- rownames(mvn.deg[[i]]) %>% intersect(rt.mt.ids)
  s <- rownames(svn.deg[[i]]) %>% intersect(rt.mt.ids)
  r <- rownames(rvn.deg[[i]]) %>% intersect(rt.mt.ids)
  ifelse(i == 5, x <- 6, x <- i)
  venn <- venn.diagram(x = list(m, s, r),
                       main = paste("Day", x, "DEG genes", sep = " "),
                       main.fontface = "bold",
                       main.fontfamily = "sans",
                       category.names = c("MT", "SBP", "RTMT"),
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
                       cat.pos = c(-27, 27, 135),
                       cat.dist = c(0.055, 0.055, 0.085),
                       cat.fontfamily = "sans",
                       cat.cex = 2.5,
                       rotation = 1);
  pdf(file = paste(results.dir, "venn/rt_mt_shared_Day_", x, ".pdf", sep = ""))
  grid.newpage()
  grid.draw(venn)
  dev.off()
}


##################################################################################################################################
# RT:MT Heatmap 
heat.genes <- rt.mt.ids
heat.dat <-as.matrix(assay(vsd)[heat.genes, ])
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
              clustering_distance_rows = 'euclidean',
              column_names_max_height = unit(8, "cm"),
              row_names_gp            = gpar(fontsize = 6),
              column_names_gp         = gpar(fontsize = 10),
              # split                   = split.vec, #This would use vector I created above to split and cluster the genes
              # top_annotation          = ha1,
              # bottom_annotation       = anno_bottom,
              col                     = colorRamp2(c(-3.5, 0, 3.5), rev(colorRampPalette(brewer.pal(9, "RdBu"))(3))));

pdf(paste(results.dir, "rt_mt_shared_HM_April.pdf", sep=""), height = 20, width = 15)
ht
dev.off()


## Investigation SBP results for RT:MT shared genes

sbp.res <- res.list[grep('^sbp.._vs_naive', names(res.list))]

names(sbp.res)

sbp.res.rt.mt.genes <- lapply(sbp.res, function(df){
  df <- df[rt.mt.ids, ]
})

par(mfrow=c(2,3))
for (res in names(sbp.res.rt.mt.genes)){
  data <- sbp.res.rt.mt.genes[[res]]
  p <- hist(data$log2FC, breaks = seq(-2,2, 0.25), main = paste("Log2FC Histogram of", res))
  }

##################################################################################################################################
# RT:MT Shared pathway

sbp.pathways <- read.delim(file = paste(results.dir, "sbp_unique_metacore_pathways.txt", sep = ""),
                             sep = "\t",
                             header = TRUE)

sbp.object_ids <- read.table(file = paste(results.dir, "sbp_unique_objects.csv", sep = ""),
                               sep = ",", 
                               header = TRUE)


sbp.path.genes <- sbp.pathways$Network.Objects.from.Active.Data

names(sbp.path.genes) <- sbp.pathways$Maps

names(sbp.path.genes)


sbp.path.symbols <- list()
for (n in names(sbp.path.genes)){
  gene.string <- sbp.path.genes[[n]]
  gene.list <- str_split(gene.string, ", ")[[1]]
  gene.names <- sapply(gene.list, function(x){
    id <- sbp.object_ids[sbp.object_ids$Network.Object.Name %in% x, "Gene.Symbol"]
    id[[1]]
  })
  sbp.path.symbols[[as.character(n)]] <- gene.names
}

sbp.pathway.df <- data.frame(pathway = names(sbp.path.symbols))
sbp.pathway.df$num_genes_pathway <- sbp.pathways$Total
sbp.pathway.df$FDR <- sbp.pathways$FDR
sbp.pathway.df$p_value <- sbp.pathways$pValue
sbp.pathway.df$genes <- sbp.path.symbols
sbp.pathway.df$genes <- vapply(sbp.pathway.df$genes, paste, collapse = ", ", character(1L))


write.table(sbp.pathway.df, file = paste(results.dir, "sbp_unique_pathway_results.txt", sep = ""), sep = "\t", quote = F, row.names = F)



##################################################################################################################################
# SBP unique genes in MT and RMT

# MT 
mt.res <- res.list[grep('^mt.._vs_naive', names(res.list))]
sbp.unique.mt.res <- lapply(mt.res, function(res){
  res <- res[sbp.unique, ] %>% drop_na()
})
##################################################################################################################################
# SBP Unique

sbp.unique <- setdiff(svn.pc, common.pc)
length(sbp.unique)
sbp.unique <- setdiff(sbp.unique, rvn.pc)
length(sbp.unique)
sbp.unique <- setdiff(sbp.unique, mvn.pc)

sbp.unique.names <- gtf.dat[sbp.unique, "gene_name"]

sbp.unique.df <- data.frame(sbp.unique.names)

write.table(sbp.unique.df, file = paste(results.dir, "sbp_unique_df.tsv", sep = ""), sep = "\t", quote = F)

#########
#SBP Unique Heatmap

heat.genes <- sbp.unique
heat.dat <-as.matrix(assay(vsd)[heat.genes, ])
sample.keep <- grep("d6", colnames(heat.dat))
heat.dat <- heat.dat[,sample.keep]
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
              clustering_distance_rows = 'euclidean',
              column_names_max_height = unit(8, "cm"),
              row_names_gp            = gpar(fontsize = 6),
              column_names_gp         = gpar(fontsize = 10),
              # split                   = split.vec, #This would use vector I created above to split and cluster the genes
              # top_annotation          = ha1,
              # bottom_annotation       = anno_bottom,
              col                     = colorRamp2(c(-3.5, 0, 3.5), rev(colorRampPalette(brewer.pal(9, "RdBu"))(3))));

ht

## Gene Profiles

gene.id <- sbp.unique
for (gid in gene.id){
  sym <- gtf.dat[gid, "gene_name"]
  plot.file <- paste(results.dir,"2-gene_profiles/rmt_sbp/",sym,".png",sep='');
  if (!file.exists(plot.file)) {
    cnt.dat  <- plotCounts(dds, which(rownames(dds) %in% gid), intgroup = c("day","delivery"), returnData = TRUE);
    cnt.trans <- cnt.dat[1:60,]
    cnt.naive <- cnt.dat[61:64, ]
    naive.mean <- mean(cnt.naive$count)
    #cnt.dat$day <- as.numeric(as.character(cnt.dat$day))
    cnt.plot <- ggplot(cnt.trans, aes(x = day, y = count, color = delivery, fill=delivery, group = delivery)) + geom_smooth(se = FALSE, alpha=0.1, method = "loess", span=0.75);
    cnt.plot <- cnt.plot + geom_hline(yintercept = naive.mean, col = "black", linetype="dashed")
    cnt.plot <- cnt.plot + ggtitle(sym) + ylab("normalised count") + xlab("day")
    cnt.plot <- cnt.plot + scale_x_discrete(limits=1:6)
    #cnt.plot <- cnt.plot + scale_y_log10();
    cnt.cols <- c("mosquito" = "#BCBD22", "blood" = "#E377C2", "naive" = "Black", 'recent_blood' = "#17BECF")
    cnt.plot <- cnt.plot + scale_colour_manual(values = cnt.cols)
    png(file=paste(results.dir, "sbp_unique_gene_profiles/" ,sym,".png",sep=''),height=500,width=700);
    print(cnt.plot);
    dev.off();
  }
}

##################################################################################################################################

### SBP vs RMT


rvs.res <- res.list[grep("rtmt.\\d_vs_sbp", names(res.list))]

rvs.res <- lapply(rvs.res,function(res){
  res <- as.data.frame(res)
  colnames(res)[colnames(res) %in% "log2FoldChange"] <- "log2FC"
  colnames(res)[colnames(res) %in% "padj"] <- "FDR"
  res$DEG <- res$FDR <= 0.05 &
    res$baseMean>=30         &
    !is.na(res$log2FC)       &
    !is.na(res$FDR)
  res
})

rvs.deg <- lapply(rvs.res, function(df){
  df <- df[df$DEG, ]
})

cbind(sapply(rvs.deg, function(x){sum(x$DEG)}))

rvs.gene.id <- unique(unlist(sapply(rvs.deg, row.names)))

rvs.pc.id <- intersect(rvs.gene.id,gtf.dat$gene_id[gtf.dat$gene_biotype %in% "protein_coding"])

rvs.gene.names <- gtf.dat[rvs.gene.id, "gene_name"]

### SBP vs RMT Gene profiles - Protein Coding
gene.id <- rvs.pc.id
for (gid in gene.id){
  sym <- gtf.dat[gid, "gene_name"]
  plot.file <- paste(results.dir,"2-gene_profiles/rmt_sbp/",sym,".png",sep='');
  if (!file.exists(plot.file)) {
    cnt.dat  <- plotCounts(dds, which(rownames(dds) %in% gid), intgroup = c("day","delivery"), returnData = TRUE);
    cnt.trans <- cnt.dat[1:60,]
    cnt.naive <- cnt.dat[61:64, ]
    naive.mean <- mean(cnt.naive$count)
    #cnt.dat$day <- as.numeric(as.character(cnt.dat$day))
    cnt.plot <- ggplot(cnt.trans, aes(x = day, y = count, color = delivery, fill=delivery, group = delivery)) + geom_smooth(se = FALSE, alpha=0.1, method = "loess", span=0.75);
    cnt.plot <- cnt.plot + geom_hline(yintercept = naive.mean, col = "black", linetype="dashed")
    cnt.plot <- cnt.plot + ggtitle(sym) + ylab("normalised count") + xlab("day")
    cnt.plot <- cnt.plot + scale_x_discrete(limits=1:6)
    #cnt.plot <- cnt.plot + scale_y_log10();
    cnt.cols <- c("mosquito" = "#BCBD22", "blood" = "#E377C2", "naive" = "Black", 'recent_blood' = "#17BECF")
    cnt.plot <- cnt.plot + scale_colour_manual(values = cnt.cols)
    png(file=paste(results.dir, "rmt_vs_sbp_gene_profiles/" ,sym,".png",sep=''),height=500,width=700);
    print(cnt.plot);
    dev.off();
  }
}


rvs.df <- data.frame(rvs.pc.id, gene_name = gtf.dat[rvs.pc.id, "gene_name"])

write.table(rvs.df, paste(results.dir, "rmt_vs_sbp_genes.tsv", sep = ""), sep = "\t", quote = F, row.names = F)


##### RMT vs SBP pathways

rvs.pathways <- read.delim(file = paste(results.dir, "rmt_sbp_metacore_pathways.txt", sep = ""),
                           sep = "\t",
                           header = TRUE)

rvs.object_ids <- read.table(file = paste(results.dir, "rmt_sbp_object_ids.txt", sep = ""),
                             sep = "\t", 
                             header = TRUE)


rvs.path.genes <- rvs.pathways$Network.Objects.from.Active.Data

names(rvs.path.genes) <- rvs.pathways$Maps

names(rvs.path.genes)


rvs.path.symbols <- list()
for (n in names(rvs.path.genes)){
  gene.string <- rvs.path.genes[[n]]
  gene.list <- str_split(gene.string, ", ")[[1]]
  gene.names <- sapply(gene.list, function(x){
    id <- rvs.object_ids[rvs.object_ids$Network.Object.Name %in% x, "Gene.Symbol"]
    id[[1]]
  })
  rvs.path.symbols[[as.character(n)]] <- gene.names
}

rvs.pathway.df <- data.frame(pathway = names(rvs.path.symbols))
rvs.pathway.df$num_genes_pathway <- rvs.pathways$Total
rvs.pathway.df$FDR <- rvs.pathways$FDR
rvs.pathway.df$p_value <- rvs.pathways$pValue
rvs.pathway.df$genes <- rvs.path.symbols
rvs.pathway.df$genes <- vapply(rvs.pathway.df$genes, paste, collapse = ", ", character(1L))


write.table(rvs.pathway.df, file = paste(results.dir, "rmt_vs_sbp_pathway_results.tsv", sep = ""), sep = "\t", quote = F, row.names = F)


###########################################################################################################################################
#### Doing direct comparisons of each transmission 

deg.list2 <- lapply(res.list, function(df){
  df <- df[df$DEG, ]
})

mvs.deg <- deg.list2[grep('^mt.._vs_sbp', names(deg.list2))]
mvr.deg <- deg.list2[grep('^mt.._vs_rtmt', names(deg.list2))]
rvs.deg <- deg.list2[grep('^rtmt.._vs_sbp', names(deg.list2))]

mvs.genes <- unique(unlist(sapply(mvs.deg, row.names)))
mvr.genes <- unique(unlist(sapply(mvr.deg, row.names)))
rvs.genes <- unique(unlist(sapply(rvs.deg, row.names)))

length(mvs.genes)
length(mvr.genes)
length(rvs.genes)

mvs.pc <- intersect(mvs.genes, gtf.dat$gene_id[gtf.dat$gene_biotype %in% "protein_coding"])
mvr.pc <- intersect(mvr.genes, gtf.dat$gene_id[gtf.dat$gene_biotype %in% "protein_coding"])
rvs.pc <- intersect(rvs.genes, gtf.dat$gene_id[gtf.dat$gene_biotype %in% "protein_coding"])

length(mvs.pc)
length(mvr.pc)
length(rvs.pc)

library("VennDiagram")
venn.cols <- c("#BCBD22", "#E377C2", "#17BECF")
venn <- venn.diagram(x = list(mvs.pc, mvr.pc, rvs.pc),
                     category.names = c("MT vs SBP", "MT vs RMT", "RMT vs SBP"),
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
                     cat.pos = c(-27, 27, 135),
                     cat.dist = c(0.055, 0.055, 0.085),
                     cat.fontfamily = "sans",
                     cat.cex = 2.5,
                     rotation = 1);

grid.newpage()
grid.draw(venn)

cross.common.pc <-Reduce(intersect, list(mvs.pc, mvr.pc, rvs.pc)) 

mvs.rvs.shared <- intersect(mvs.pc, rvs.pc)
# rvs.unique <- setdiff(mvs.rvs.shared, cross.common.pc)


gene.id <- cross.common.pc
for (gid in gene.id){
  sym <- gtf.dat[gid, "gene_name"]
  plot.file <- paste(results.dir,"cross_comparison/mvs_rvs_shared/",sym,".png",sep='');
  if (!file.exists(plot.file)) {
    cnt.dat  <- plotCounts(dds, which(rownames(dds) %in% gid), intgroup = c("day","delivery"), returnData = TRUE);
    cnt.trans <- cnt.dat[1:60,]
    cnt.naive <- cnt.dat[61:64, ]
    naive.mean <- mean(cnt.naive$count)
    #cnt.dat$day <- as.numeric(as.character(cnt.dat$day))
    cnt.plot <- ggplot(cnt.trans, aes(x = day, y = count, color = delivery, fill=delivery, group = delivery)) + geom_smooth(se = FALSE, alpha=0.1, method = "loess", span=0.75);
    cnt.plot <- cnt.plot + geom_hline(yintercept = naive.mean, col = "black", linetype="dashed")
    cnt.plot <- cnt.plot + ggtitle(sym) + ylab("normalised count") + xlab("day")
    cnt.plot <- cnt.plot + scale_x_discrete(limits=1:6)
    #cnt.plot <- cnt.plot + scale_y_log10();
    cnt.cols <- c("mosquito" = "#BCBD22", "blood" = "#E377C2", "naive" = "Black", 'recent_blood' = "#17BECF")
    cnt.plot <- cnt.plot + scale_colour_manual(values = cnt.cols)
    png(file=paste(results.dir, "cross_comparison/mvs_rvs_shared/" ,sym,".png",sep=''),height=500,width=700);
    print(cnt.plot);
    dev.off();
  }
}

mvs.rvs.df <- data.frame(id = mvs.rvs.shared, gene_name = gtf.dat[mvs.rvs.shared, "gene_name"])

write.table(mvs.rvs.df, file = paste(results.dir, "mvs_rvs_shared_genes.tsv", sep = ""), sep = "\t", row.names = F, quote = F)


### Metacore results

mvs_rvs.pathways <- read.delim(file = paste(results.dir, "mvs_rvs_metacore_pathways.txt", sep = ""),
                           sep = "\t",
                           header = TRUE)

mvs_rvs.object_ids <- read.table(file = paste(results.dir, "mvs_rvs_object_ids.txt", sep = ""),
                             sep = "\t", 
                             header = TRUE)


mvs_rvs.path.genes <- mvs_rvs.pathways$Network.Objects.from.Active.Data

names(mvs_rvs.path.genes) <- mvs_rvs.pathways$Maps

names(mvs_rvs.path.genes)


mvs_rvs.path.symbols <- list()
for (n in names(mvs_rvs.path.genes)){
  gene.string <- mvs_rvs.path.genes[[n]]
  gene.list <- str_split(gene.string, ", ")[[1]]
  gene.names <- sapply(gene.list, function(x){
    id <- mvs_rvs.object_ids[mvs_rvs.object_ids$Network.Object.Name %in% x, "Gene.Symbol"]
    id[[1]]
  })
  mvs_rvs.path.symbols[[as.character(n)]] <- gene.names
}

mvs_rvs.pathway.df <- data.frame(pathway = names(mvs_rvs.path.symbols))
mvs_rvs.pathway.df$num_genes_pathway <- mvs_rvs.pathways$Total
mvs_rvs.pathway.df$FDR <- mvs_rvs.pathways$FDR
mvs_rvs.pathway.df$p_value <- mvs_rvs.pathways$pValue
mvs_rvs.pathway.df$genes <- mvs_rvs.path.symbols
mvs_rvs.pathway.df$genes <- vapply(mvs_rvs.pathway.df$genes, paste, collapse = ", ", character(1L))

write.table(mvs_rvs.pathway.df, file = paste(results.dir, "mvs_rvs_pathway_results.tsv", sep = ""), sep = "\t", row.names = F, quote = F)


#### Heatmap
heat.genes <- mvs.rvs.shared
heat.dat <-as.matrix(assay(vsd)[heat.genes, ])
# samples.rmv <- grep("d1|d2|d3|d4", x = colnames(heat.dat))
# heat.dat <- heat.dat[ ,-samples.rmv]
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
              clustering_distance_rows = 'euclidean',
              column_names_max_height = unit(8, "cm"),
              row_names_gp            = gpar(fontsize = 10),
              column_names_gp         = gpar(fontsize = 10),
              # split                   = split.vec, #This would use vector I created above to split and cluster the genes
              # top_annotation          = ha1,
              # bottom_annotation       = anno_bottom,
              col                     = colorRamp2(c(-3.5, 0, 3.5), rev(colorRampPalette(brewer.pal(9, "RdBu"))(3))));

pdf(paste(results.dir, "mvs_rvs_heatmap.pdf", sep=""), height = 10, width = 15)
ht
dev.off()

######## RVS Unique

rvs_unique <- setdiff(rvs.pc, mvs.rvs.shared)

rvs_unique.df <- data.frame(id = rvs_unique, gene_name = gtf.dat[rvs_unique, "gene_name"])

write.table(rvs_unique.df, file = paste(results.dir, "rvs_unique_df.tsv", sep = ""), sep = "\t", row.names = F, quote = F)

### Metacore results

rvs_unique.pathways <- read.delim(file = paste(results.dir, "rvs_unique_metacore_pathways.txt", sep = ""),
                               sep = "\t",
                               header = TRUE)

rvs_unique.object_ids <- read.table(file = paste(results.dir, "rvs_unique_object_ids.txt", sep = ""),
                                 sep = "\t", 
                                 header = TRUE)


rvs_unique.path.genes <- rvs_unique.pathways$Network.Objects.from.Active.Data

names(rvs_unique.path.genes) <- rvs_unique.pathways$Maps

names(rvs_unique.path.genes)


rvs_unique.path.symbols <- list()
for (n in names(rvs_unique.path.genes)){
  gene.string <- rvs_unique.path.genes[[n]]
  gene.list <- str_split(gene.string, ", ")[[1]]
  gene.names <- sapply(gene.list, function(x){
    id <- rvs_unique.object_ids[rvs_unique.object_ids$Network.Object.Name %in% x, "Gene.Symbol"]
    id[[1]]
  })
  rvs_unique.path.symbols[[as.character(n)]] <- gene.names
}

rvs_unique.pathway.df <- data.frame(pathway = names(rvs_unique.path.symbols))
rvs_unique.pathway.df$num_genes_pathway <- rvs_unique.pathways$Total
rvs_unique.pathway.df$FDR <- rvs_unique.pathways$FDR
rvs_unique.pathway.df$p_value <- rvs_unique.pathways$pValue
rvs_unique.pathway.df$genes <- rvs_unique.path.symbols
rvs_unique.pathway.df$genes <- vapply(rvs_unique.pathway.df$genes, paste, collapse = ", ", character(1L))

write.table(rvs_unique.pathway.df, file = paste(results.dir, "rvs_unique_pathway_results.tsv", sep = ""), sep = "\t", row.names = F, quote = F)


#### gene_profiles

gene.id <- rvs_unique
for (gid in gene.id){
  sym <- gtf.dat[gid, "gene_name"]
  cnt.dat  <- plotCounts(dds, which(rownames(dds) %in% gid), intgroup = c("day","delivery"), returnData = TRUE);
  cnt.trans <- cnt.dat[1:60,]
  cnt.naive <- cnt.dat[61:64, ]
  naive.mean <- mean(cnt.naive$count)
  #cnt.dat$day <- as.numeric(as.character(cnt.dat$day))
  cnt.plot <- ggplot(cnt.trans, aes(x = day, y = count, color = delivery, fill=delivery, group = delivery)) + geom_smooth(se = FALSE, alpha=0.1, method = "loess", span=0.75);
  cnt.plot <- cnt.plot + geom_hline(yintercept = naive.mean, col = "black", linetype="dashed")
  cnt.plot <- cnt.plot + ggtitle(sym) + ylab("normalised count") + xlab("day")
  cnt.plot <- cnt.plot + scale_x_discrete(limits=1:6)
  #cnt.plot <- cnt.plot + scale_y_log10();
  cnt.cols <- c("mosquito" = "#BCBD22", "blood" = "#E377C2", "naive" = "Black", 'recent_blood' = "#17BECF")
  cnt.plot <- cnt.plot + scale_colour_manual(values = cnt.cols)
  png(file=paste(results.dir, "rvs_unique_gene_profiles/" ,sym,".png",sep=''),height=500,width=700);
  print(cnt.plot);
  dev.off();
}



#### Heatmap
heat.genes <- rvs_unique
heat.dat <-as.matrix(assay(vsd)[heat.genes, ])
# samples.rmv <- grep("d1|d2|d3|d4", x = colnames(heat.dat))
# heat.dat <- heat.dat[ ,-samples.rmv]
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
              clustering_distance_rows = 'euclidean',
              column_names_max_height = unit(8, "cm"),
              row_names_gp            = gpar(fontsize = 10),
              column_names_gp         = gpar(fontsize = 10),
              # split                   = split.vec, #This would use vector I created above to split and cluster the genes
              # top_annotation          = ha1,
              # bottom_annotation       = anno_bottom,
              col                     = colorRamp2(c(-3.5, 0, 3.5), rev(colorRampPalette(brewer.pal(9, "RdBu"))(3))));

pdf(paste(results.dir, "rvs_unique_heatmap.pdf", sep=""), height = 10, width = 15)
ht
dev.off()







library("VennDiagram")
venn.cols <- c("#BCBD22", "#E377C2")
venn <- venn.diagram(x = list(mvs.pc, rvs.pc),
                     category.names = c("MT vs SBP", "RMT vs SBP"),
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
                     cat.pos = c(-27, 27),
                     cat.dist = c(0.055, 0.055),
                     cat.fontfamily = "sans",
                     cat.cex = 2.5,
                     scaled = FALSE);

pdf(file = paste(results.dir, "venn/mvs_rvs_venn.pdf", sep = ""))
grid.newpage()
grid.draw(venn)
dev.off()


 mvs.unique <- setdiff(mvs.pc, mvr.pc)

mvs.rvs.shared2 <- intersect(mvs.pc, rvs.pc)

gene.id <- mvs.rvs.shared2
for (gid in gene.id){
  sym <- gtf.dat[gid, "gene_name"]
  plot.file <- paste(results.dir,"cross_comparison/mvs_rvs_shared2/",sym,".png",sep='');
  if (!file.exists(plot.file)) {
    cnt.dat  <- plotCounts(dds, which(rownames(dds) %in% gid), intgroup = c("day","delivery"), returnData = TRUE);
    cnt.trans <- cnt.dat[1:60,]
    cnt.naive <- cnt.dat[61:64, ]
    naive.mean <- mean(cnt.naive$count)
    #cnt.dat$day <- as.numeric(as.character(cnt.dat$day))
    cnt.plot <- ggplot(cnt.trans, aes(x = day, y = count, color = delivery, fill=delivery, group = delivery)) + geom_smooth(se = FALSE, alpha=0.1, method = "loess", span=0.75);
    cnt.plot <- cnt.plot + geom_hline(yintercept = naive.mean, col = "black", linetype="dashed")
    cnt.plot <- cnt.plot + ggtitle(sym) + ylab("normalised count") + xlab("day")
    cnt.plot <- cnt.plot + scale_x_discrete(limits=1:6)
    #cnt.plot <- cnt.plot + scale_y_log10();
    cnt.cols <- c("mosquito" = "#BCBD22", "blood" = "#E377C2", "naive" = "Black", 'recent_blood' = "#17BECF")
    cnt.plot <- cnt.plot + scale_colour_manual(values = cnt.cols)
    png(file=paste(results.dir, "cross_comparison/mvs_rvs_shared2/" ,sym,".png",sep=''),height=500,width=700);
    print(cnt.plot);
    dev.off();
  }
}




############# FULL RANK DESEQ

design.fr.file <- paste(nf.dir,"design_full_rank.csv",sep="")
design.fr.baseline <- "blood" # This is the treatment group you want to use for baseline measurements and test others against - Usually control

## Design file
design.fr           <- read.delim(design.fr.file, header=TRUE, sep=",", stringsAsFactors = TRUE)
rownames(design.fr) <- design.fr$label
design.fr$delivery  <- relevel(design.fr$delivery, design.fr.baseline) # sets design.baseline as first factor (baseline)
design.fr$day       <- factor(as.character(design.fr$day), levels=as.character(sort(unique(design.fr$day)))) # factors days in numberical order
design.fr$replicate <- as.factor(design.fr$replicate) # factors replicates (order not important)

## Importing counts
txi <- tximport(paste(nf.dir, "results/rsem/", design.fr$lims.name, ".genes.results", sep=""), type="rsem")
txi$length[txi$length==0] <- 1


dat.fr <- DESeqDataSetFromTximport(txi, design.fr, ~ delivery + day + delivery:day) # For baseline counts - full rank in other file
rowData(dat.fr) <- gtf.dat[rownames(dat.fr),c("gene_id","gene_name","gene_source","gene_biotype","GRCm38")]
dds.fr <- DESeq(dat.fr) # non-converging genes likely contamination from pancreatic tissue (remove man)

dds_lrt <- DESeq(dds.fr, test="LRT", reduced = ~ delivery + day)
## DEseq vignette
res_lrt <- results(dds_lrt)
topGenes <- head(order(res_lrt$padj),20)
mat <- betas[topGenes, -c(1,2)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)

res.lrt <- as.data.frame(results(dds_lrt))
colnames(res.lrt)[colnames(res.lrt) %in% "log2FoldChange"] <- "log2FC"
colnames(res.lrt)[colnames(res.lrt) %in% "padj"] <- "FDR"
res.lrt$DEG <- res.lrt$FDR <= 0.05 &
  res.lrt$baseMean>=30         &
  !is.na(res.lrt$log2FC)       &
  !is.na(res.lrt$FDR)

deg.lrt <- res.lrt[res.lrt$DEG, ]
deg.lrt <- deg.lrt[order(deg.lrt$log2FC, decreasing = TRUE), ]
deg.lrt$gene_sym <- gtf.dat[rownames(deg.lrt), "gene_name"]
deg.lrt$gene_sym


deg.lrt <- deg.lrt[order(deg.lrt$FDR), ]

count.plot <- function(gene_id){
  x <- plotCounts(dds.fr, gene_id,
                  intgroup = c("day","delivery"), returnData = TRUE) 
  
  ggplot(x,
         aes(x = day, y = count, color = delivery, group = delivery)) + 
    geom_point() + stat_summary(fun=mean, geom="line")
  
}

gene.id <- rownames(deg.lrt)
for (gid in gene.id){
  sym <- gtf.dat[gid, "gene_name"]
  plot.file <- paste(results.dir,"lrt/",sym,".png",sep='');
  if (!file.exists(plot.file)) {
    cnt.dat  <- plotCounts(dds, which(rownames(dds) %in% gid), intgroup = c("day","delivery"), returnData = TRUE);
    cnt.trans <- cnt.dat[1:60,]
    cnt.naive <- cnt.dat[61:64, ]
    naive.mean <- mean(cnt.naive$count)
    #cnt.dat$day <- as.numeric(as.character(cnt.dat$day))
    cnt.plot <- ggplot(cnt.trans, aes(x = day, y = count, color = delivery, fill=delivery, group = delivery)) + geom_smooth(se = FALSE, alpha=0.1, method = "loess", span=0.75);
    cnt.plot <- cnt.plot + geom_hline(yintercept = naive.mean, col = "black", linetype="dashed")
    cnt.plot <- cnt.plot + ggtitle(sym) + ylab("normalised count") + xlab("day")
    cnt.plot <- cnt.plot + scale_x_discrete(limits=1:6)
    #cnt.plot <- cnt.plot + scale_y_log10();
    cnt.cols <- c("mosquito" = "#BCBD22", "blood" = "#E377C2", "naive" = "Black", 'recent_blood' = "#17BECF")
    cnt.plot <- cnt.plot + scale_colour_manual(values = cnt.cols)
    png(file=paste(results.dir, "lrt/" ,sym,".png",sep=''),height=500,width=700);
    print(cnt.plot);
    dev.off();
  }
}



#### Clustering R Log 
library(DEGreport)
deg_lrt.rlog <- assay(vsd.fr[rownames(deg.lrt), ])
nrow(deg_lrt.rlog)
# clusters.test <- degPatterns(deg_lrt.rlog, meta = design.fr, time = "day", col="delivery",
                             # consensusCluster = TRUE)

clusters <- degPatterns(deg_lrt.rlog, meta = design.fr, time = "day", col="delivery")

clusters$hr
cluster_groups <- clusters$df

cluster.plot <- degPlotCluster(clusters$normalized, time = "day", color = "delivery")

cluster.plot <- cluster.plot + scale_colour_manual(values = cnt.cols)

pdf(paste(results.dir, "DEG_report_cluster_facet.pdf", sep = ""), height = 14, width = 14)
plot(cluster.plot)
dev.off()

group25 <- cluster_groups %>% 
  filter(cluster == 25)

group12 <- cluster_groups %>% 
  filter(cluster == 12)

group14 <- cluster_groups %>% 
  filter(cluster == 14)

group4 <- cluster_groups %>% 
  filter(cluster == 4)

# grp25.plots <- list()
# for (g in group25$genes){
#   plot <- count.plot(dds, g)
#   grp25.plots[[as.character(g)]] <- plot
# }
# grp25.plots$ncol <- 5
# do.call(plot_grid, grp25.plots)  
# 
# grp4.plots <- list()
#  for (g in group4$genes){
#   plot <- count.plot(dds, g)
#   grp4.plots[[as.character(g)]] <- plot
# }
# grp4.plots$ncol <- 5
# do.call(plot_grid, grp4.plots)  
# 
# grp14.plots <- list()
# for (g in group14$genes){
#   plot <- count.plot(dds, g)
#   grp14.plots[[as.character(g)]] <- plot
# }
# grp14.plots$ncol <- 5
# do.call(plot_grid, grp14.plots)  
# 
# grp14.plots <- list()
# for (g in group14$genes){
#   plot <- count.plot(dds, g)
#   grp14.plots[[as.character(g)]] <- plot
# }
# grp14.plots$ncol <- 5
# do.call(plot_grid, grp14.plots)  

group.deg <- cluster_groups %>% 
  filter(cluster == 4 | cluster == 12 | cluster == 14 | cluster == 25)

group.deg$gene_name <- sapply(group.deg$genes, function(x){gtf.dat[x, "gene_name"]})

write.table(group.deg, file = paste(results.dir, "DEGReport_group_df.tsv", sep = ""), sep = "\t", row.names = F, quote = F)

apply(group.deg, 1, function(id){
  gene_name <- id[3]
  pdf(paste(results.dir, "DEG_Report_profiles/", gene_name, ".pdf", sep = ""))
  print(count.plot(dds, id[1]))
  dev.off()
})





# group25$gene_name <- gtf.dat[group25$genes, "gene_name"]
# group25
# sig.groups <- rbind(group25, group4)
# mvs.rvs.shared %in% sig.groups$genes
# rownames(cluster_groups) <- cluster_groups$genes
# 
# cluster_groups[mvs.rvs.shared, ] %>% drop_na()
# 
# test <- count.plot(dds, il18)
# 
# ## DEG Report, Protein Coding
# lrt.deg.pc <- intersect(rownames(deg.lrt),gtf.dat$gene_id[gtf.dat$gene_biotype %in% "protein_coding"])
# lrt.pc.rlog <- assay(vsd.fr[lrt.deg.pc, ])
# clusters.pc<- degPatterns(lrt.pc.rlog, meta = design.fr, time = "day", col="delivery")
# clus.pc_groups <- clusters.pc$df
# clus.pc.df[mvs.rvs.shared, ] %>% drop_na()
# 
# group.pc <- clus.pc_groups %>% 
#   filter(cluster == 4|cluster == 5)
# 
# group.pc$gene_name <- sapply(group.pc$genes, function(x){gtf.dat[x, "gene_name"]})
# 
# apply(group.pc, 1, function(id){
#   gene_name <- id[3]
#   pdf(paste(results.dir, "DEG_Report_PC/", gene_name, ".pdf", sep = ""))
#   print(count.plot(dds, id[1]))
#   dev.off()
# })

#### Heatmap
heat.genes <- group.deg$genes
# heat.genes <- intersect(heat.genes,gtf.dat$gene_id[gtf.dat$gene_biotype %in% "protein_coding"])
heat.dat <-as.matrix(assay(vsd)[heat.genes, ])
# heat.dat   <- t(apply(heat.dat,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))}))
# heat.columns <- grep("naive|d6", colnames(heat.dat))
# heat.dat <- heat.dat[,heat.columns]
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
              row_names_gp            = gpar(fontsize = 8),
              column_names_gp         = gpar(fontsize = 12),
              # split                   = split.vec, #This would use vector I created above to split and cluster the genes
              # top_annotation          = ha1,
              # bottom_annotation       = anno_bottom,
              col                     = colorRamp2(c(min, 0, max), rev(colorRampPalette(brewer.pal(9, "RdBu"))(3))));

pdf(paste(results.dir, "lrt_all_DEGReport_heatmap.pdf", sep=""), height = 25, width = 15)
ht
decorate_heatmap_body("zscore", c(for(d in seq(0, n, 4)){
  grid.lines(x = c(d/n, d/n), y = c(0,1), gp= gpar(col = "black", lty="solid", lwd=2))},
  grid.lines(x = c(0,1), y = c(0,0), gp= gpar(col = "black", lty="solid", lwd=2)),
  grid.lines(x = c(0,1), y = c(1,1), gp= gpar(col = "black", lty="solid", lwd=2))))
dev.off()

### LRT Metacore Pathways
lrt_deg.pathways <- read.delim(file = paste(results.dir, "EKD11_LRT_DEGReport_metacore_pathways.txt", sep = ""),
                           sep = "\t",
                           header = TRUE)

lrt_deg.object_ids <- read.table(file = paste(results.dir, "EKD11_LRT_DEGReport_network_ids.txt", sep = ""),
                             sep = "\t", 
                             header = TRUE)


lrt_deg.path.genes <- lrt_deg.pathways$Network.Objects.from.Active.Data
names(lrt_deg.path.genes) <- lrt_deg.pathways$Maps

names(lrt_deg.path.genes)


lrt_deg.path.symbols <- list()
for (n in names(lrt_deg.path.genes)){
  gene.string <- lrt_deg.path.genes[[n]]
  gene.list <- str_split(gene.string, ", ")[[1]]
  gene.names <- sapply(gene.list, function(x){
    id <- lrt_deg.object_ids[lrt_deg.object_ids$Network.Object.Name %in% x, "Gene.Symbol"]
    id[[1]]
  })
  lrt_deg.path.symbols[[as.character(n)]] <- gene.names
}

lrt_deg.pathway.df <- data.frame(pathway = names(lrt_deg.path.symbols))
lrt_deg.pathway.df$num_genes_pathway <- lrt_deg.pathways$Total
lrt_deg.pathway.df$FDR <- lrt_deg.pathways$FDR
lrt_deg.pathway.df$p_value <- lrt_deg.pathways$pValue
lrt_deg.pathway.df$genes <- lrt_deg.path.symbols
lrt_deg.pathway.df$genes <- vapply(lrt_deg.pathway.df$genes, paste, collapse = ", ", character(1L))

lrt_deg.gene_list <- unique(unlist(lrt_deg.path.symbols))
write.table(lrt_deg.pathway.df, file = paste(results.dir, "lrt_deg_pathway_results.txt", sep = ""), sep = "\t", quote = F, row.names = F)

lrt_deg.pathway.sig <- lrt_deg.pathway.df[lrt_deg.pathway.df$FDR < 0.05, ]

lrt_deg.symbols.sig <- unlist(lrt_deg.path.symbols[lrt_deg.pathway.sig$pathway])
names(rev(sort(table(lrt_deg.symbols.sig))))



for (g in lrt_deg.gene_list){
  gene_id <- gtf.dat[gtf.dat$gene_name %in% g, "gene_id"]
  pdf(paste(results.dir, "LRT_metacore_profiles/", g, ".pdf", sep = ""))
  print(count.plot(dds, gene_id))
  dev.off()
}

lrt_metacore.plots <- list()
for (g in lrt_deg.gene_list){
  gene_id <- gtf.dat[gtf.dat$gene_name %in% g, "gene_id"]
  plot <- count.plot(dds, gene_id)
  lrt_metacore.plots[[as.character(g)]] <- plot
}
lrt_metacore.plots$ncol <- 4
pdf(paste(results.dir, "LRT_metacore_profile_facet.pdf", sep = ""), height = 60, width = 40)
do.call(plot_grid, lrt_metacore.plots)
dev.off()

lrt_metacore.sig.plots <- list()
lrt_meta_sig.genes <- names(rev(sort(table(lrt_deg.symbols.sig))))
for (g in lrt_meta_sig.genes){
  gene_id <- gtf.dat[gtf.dat$gene_name %in% g, "gene_id"]
  plot <- count.plot(dds, gene_id)
  lrt_metacore.sig.plots[[as.character(g)]] <- plot
}
lrt_metacore.sig.plots$ncol <- 4
pdf(paste(results.dir, "LRT_metacore_sig_profile_facet.pdf", sep = ""), height = 40, width = 30)
do.call(plot_grid, lrt_metacore.sig.plots)
dev.off()

lrt_metacore_sig.ids <- sapply(lrt_meta_sig.genes, function(x){gtf.dat[gtf.dat$gene_name %in% x, "gene_id"]})
heat.genes <- lrt_metacore_sig.ids
heat.dat <-as.matrix(assay(vsd)[heat.genes, ])
heat.columns <- grep("naive|d6", colnames(heat.dat))
heat.dat <- heat.dat[,heat.columns]
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
              row_names_gp            = gpar(fontsize = 8),
              column_names_gp         = gpar(fontsize = 12),
              # split                   = split.vec, #This would use vector I created above to split and cluster the genes
              # top_annotation          = ha1,
              # bottom_annotation       = anno_bottom,
              col                     = colorRamp2(c(min, 0, max), rev(colorRampPalette(brewer.pal(9, "RdBu"))(3))));

pdf(paste(results.dir, "lrt_metacore_sig_d6_heatmap.pdf", sep=""), height = 15, width = 15)
ht
decorate_heatmap_body("zscore", c(for(d in seq(0, n, 4)){
  grid.lines(x = c(d/n, d/n), y = c(0,1), gp= gpar(col = "black", lty="solid", lwd=2))},
  grid.lines(x = c(0,1), y = c(0,0), gp= gpar(col = "black", lty="solid", lwd=2)),
  grid.lines(x = c(0,1), y = c(1,1), gp= gpar(col = "black", lty="solid", lwd=2))))
dev.off()


#### Heatmap
heat.genes <- rownames(deg.lrt)
heat.genes <- intersect(heat.genes,gtf.dat$gene_id[gtf.dat$gene_biotype %in% "protein_coding"])
heat.dat <-as.matrix(assay(vsd.fr)[heat.genes, ])
# heat.dat   <- t(apply(heat.dat,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))}))
heat.columns <- grep("d6", colnames(heat.dat))
heat.dat <- heat.dat[,heat.columns]
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
              row_names_gp            = gpar(fontsize = 8),
              column_names_gp         = gpar(fontsize = 12),
              # split                   = split.vec, #This would use vector I created above to split and cluster the genes
              # top_annotation          = ha1,
              # bottom_annotation       = anno_bottom,
              col                     = colorRamp2(c(min, 0, max), rev(colorRampPalette(brewer.pal(9, "RdBu"))(3))));

pdf(paste(results.dir, "lrt_d6_no_control_heatmap.pdf", sep=""), height = 150, width = 15)
ht
dev.off()

manual.df <- read.table(paste(results.dir, "EKD11_LRT_manual_gene_list.txt", sep = ""), header = F, col.names = ("gene_name"))

manual.df$gene_id <- sapply(manual.df$gene_name, function(x){gtf.dat[gtf.dat$gene_name %in% x, "gene_id"]})

apply(manual.df, 1, function(x){
  gene_name <- x[1]
  pdf(paste(results.dir, "LRT_DEG_manual_profiles/", gene_name, ".pdf", sep = ""))
  count.plot(dds, x[2])
  dev.off()
})

manual_2.df <- read.table(paste(results.dir, "EKD11_LRT_manual_gene_list2.txt", sep = ""), header = F, col.names = ("gene_name"))
manual_2.df$gene_id <- sapply(manual_2.df$gene_name, function(x){gtf.dat[gtf.dat$gene_name %in% x, "gene_id"]})


apply(manual_2.df, 1, function(x){
  gene_name <- x[1]
  pdf(paste(results.dir, "LRT_DEG_manual2_profiles/", gene_name, ".pdf", sep = ""))
  count.plot(dds, x[2])
  dev.off()
})


## Design file - RMT VS SBP only
design.rvs <- design.fr
design.rvs <- design.rvs[-grep("^mt", rownames(design.rvs)),]
design.rvs$delivery <- factor(design.rvs$delivery)
design.rvs$day       <- factor(as.character(design.rvs$day), levels=as.character(sort(unique(design.rvs$day)))) # factors days in numberical order
design.rvs$replicate <- as.factor(design.rvs$replicate) # factors replicates (order not important)

## Importing counts
txi.rvs <- tximport(paste(nf.dir, "results/rsem/", design.rvs$lims.name, ".genes.results", sep=""), type="rsem")
txi.rvs$length[txi.rvs$length==0] <- 1


dat.rvs <- DESeqDataSetFromTximport(txi.rvs, design.rvs, ~ delivery + day + delivery:day) # For baseline counts - full rank in other file
rowData(dat.rvs) <- gtf.dat[rownames(dat.rvs),c("gene_id","gene_name","gene_source","gene_biotype","GRCm38")]
dds.rvs <- DESeq(dat.rvs)

rvs.lrt <- DESeq(dds.rvs, test="LRT", reduced = ~ delivery + day)

rvs.lrt.res <- results(rvs.lrt)
rvs.lrt.res <- as.data.frame(rvs.lrt.res)
colnames(rvs.lrt.res)[colnames(rvs.lrt.res) %in% "log2FoldChange"] <- "log2FC"
colnames(rvs.lrt.res)[colnames(rvs.lrt.res) %in% "padj"] <- "FDR"
rvs.lrt.res$DEG <- rvs.lrt.res$FDR <= 0.05 &
  rvs.lrt.res$baseMean>=30         &
  !is.na(rvs.lrt.res$log2FC)       &
  !is.na(rvs.lrt.res$FDR)

rvs.lrt.deg <- rvs.lrt.res[rvs.lrt.res$DEG, ]
rvs.lrt.deg$gene_name <- gtf.dat[rownames(rvs.lrt.deg), "gene_name"]


#Top 50 genes
top.rvs.res <- rvs.lrt.res[order(rvs.lrt.res$FDR), ]
top.rvs.res <- top.rvs.res[1:50, ]
top.rvs.res$gene_name  <- gtf.dat[rownames(top.rvs.res), "gene_name"]
top.rvs.res


Igsf8 <-  plotCounts(dds.fr, "ENSMUSG00000038034",
                     intgroup = c("day","delivery"), returnData = TRUE)
ggplot(Igsf8,
       aes(x = day, y = count, color = delivery, group = delivery)) + 
  geom_point() + stat_summary(fun=mean, geom="line") +
  scale_y_log10()


Dync1 <- plotCounts(dds.fr, "ENSMUSG00000029757",
                    intgroup = c("day","delivery"), returnData = TRUE)

ggplot(Dync1,
       aes(x = day, y = count, color = delivery, group = delivery)) + 
  geom_point() + stat_summary(fun=mean, geom="line") +
  scale_y_log10()

Hexa <-  plotCounts(dds.fr, "ENSMUSG00000025232",
                    intgroup = c("day","delivery"), returnData = TRUE)

ggplot(Hexa,
       aes(x = day, y = count, color = delivery, group = delivery)) + 
  geom_point() + stat_summary(fun=mean, geom="line") +
  scale_y_log10()

Rbm3 <- plotCounts(dds.fr, "ENSMUSG00000031167",
                   intgroup = c("day","delivery"), returnData = TRUE)

ggplot(Rbm3,
       aes(x = day, y = count, color = delivery, group = delivery)) + 
  geom_point() + stat_summary(fun=mean, geom="line") +
  scale_y_log10()

Nfkbia <- plotCounts(dds.fr, "ENSMUSG00000021025",
                     intgroup = c("day","delivery"), returnData = TRUE)

ggplot(Nfkbia,
       aes(x = day, y = count, color = delivery, group = delivery)) + 
  geom_point() + stat_summary(fun=mean, geom="line") +
  scale_y_log10()

cd86 <- plotCounts(dds.fr, "ENSMUSG00000022901",
                   intgroup = c("day","delivery"), returnData = TRUE)

ggplot(cd86,
       aes(x = day, y = count, color = delivery, group = delivery)) + 
  geom_point() + stat_summary(fun=mean, geom="line")



## Delivery - delivery comparison

mt.res <- list()

mt.res[["mt1.vs.mt0"]] <- results(dds, contrast=c("grp","mt.1","naive"))
mt.res[["mt2.vs.mt1"]] <- results(dds, contrast=c("grp", "mt.2", "mt.1"))
mt.res[["mt3.vs.mt2"]] <- results(dds, contrast=c("grp", "mt.3", "mt.2"))
mt.res[["mt4.vs.mt3"]] <- results(dds, contrast=c("grp", "mt.4", "mt.3"))
mt.res[["mt6.vs.mt4"]] <- results(dds, contrast=c("grp", "mt.6", "mt.4"))

mt.res <- lapply(mt.res,function(res){
  res <- as.data.frame(res)
  colnames(res)[colnames(res) %in% "log2FoldChange"] <- "log2FC"
  colnames(res)[colnames(res) %in% "padj"] <- "FDR"
  res$DEG <- res$FDR <= 0.05 &
    # abs(res$log2FC)>= 1 &  #only filters out logFC ±1 N.B CHANGE FOR RMT vs SBP COMPARISONS
    res$baseMean>=30         &
    !is.na(res$log2FC)       &
    !is.na(res$FDR)
  res
})

cbind(sapply(mt.res, function(x){sum(x$DEG)}))

rtmt.res <- list()

rtmt.res[["rtmt1.vs.rtmt0"]] <- results(dds, contrast=c("grp","rtmt.1","naive"))
rtmt.res[["rtmt2.vs.rtmt1"]] <- results(dds, contrast=c("grp", "rtmt.2", "rtmt.1"))
rtmt.res[["rtmt3.vs.rtmt2"]] <- results(dds, contrast=c("grp", "rtmt.3", "rtmt.2"))
rtmt.res[["rtmt4.vs.rtmt3"]] <- results(dds, contrast=c("grp", "rtmt.4", "rtmt.3"))
rtmt.res[["rtmt6.vs.rtmt4"]] <- results(dds, contrast=c("grp", "rtmt.6", "rtmt.4"))

rtmt.res <- lapply(rtmt.res,function(res){
  res <- as.data.frame(res)
  colnames(res)[colnames(res) %in% "log2FoldChange"] <- "log2FC"
  colnames(res)[colnames(res) %in% "padj"] <- "FDR"
  res$DEG <- res$FDR <= 0.05 &
    # abs(res$log2FC)>= 1 &  #only filters out logFC ±1 N.B CHANGE FOR RMT vs SBP COMPARISONS
    res$baseMean>=30         &
    !is.na(res$log2FC)       &
    !is.na(res$FDR)
  res
})

cbind(sapply(rtmt.res, function(x){sum(x$DEG)}))


sbp.res <- list()

sbp.res[["sbp1.vs.sbp0"]] <- results(dds, contrast=c("grp","sbp.1","naive"))
sbp.res[["sbp2.vs.sbp1"]] <- results(dds, contrast=c("grp", "sbp.2", "sbp.1"))
sbp.res[["sbp3.vs.sbp2"]] <- results(dds, contrast=c("grp", "sbp.3", "sbp.2"))
sbp.res[["sbp4.vs.sbp3"]] <- results(dds, contrast=c("grp", "sbp.4", "sbp.3"))
sbp.res[["sbp6.vs.sbp4"]] <- results(dds, contrast=c("grp", "sbp.6", "sbp.4"))

sbp.res <- lapply(sbp.res,function(res){
  res <- as.data.frame(res)
  colnames(res)[colnames(res) %in% "log2FoldChange"] <- "log2FC"
  colnames(res)[colnames(res) %in% "padj"] <- "FDR"
  res$DEG <- res$FDR <= 0.05 &
    # abs(res$log2FC)>= 1 &  #only filters out logFC ±1 N.B CHANGE FOR RMT vs SBP COMPARISONS
    res$baseMean>=30         &
    !is.na(res$log2FC)       &
    !is.na(res$FDR)
  res
})

cbind(sapply(sbp.res, function(x){sum(x$DEG)}))




mt.deg <- unique(unlist(sapply(mt.res, function(x){rownames(x[x$DEG,])})))
rtmt.deg <- unique(unlist(sapply(rtmt.res, function(x){rownames(x[x$DEG, ])})))
sbp.deg <- unique(unlist(sapply(sbp.res, function(x){rownames(x[x$DEG, ])})))

library("VennDiagram")
venn.cols <- c("#BCBD22", "#E377C2", "#17BECF")
venn <- venn.diagram(x = list(mt.deg, sbp.deg, rtmt.deg),
                     category.names = c("MT", "SBP", "RTMT"),
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
                     cat.pos = c(-27, 27, 135),
                     cat.dist = c(0.055, 0.055, 0.085),
                     cat.fontfamily = "sans",
                     cat.cex = 2.5,
                     rotation = 1);

# pdf(file = paste(results.dir, "venn/All_newcols.pdf", sep = ""))
grid.newpage()
grid.draw(venn)
dev.off()

common.deg <- Reduce(intersect, list(mt.deg, rtmt.deg, sbp.deg))

rtmt.mt.common <- intersect(mt.deg,rtmt.deg)
rtmt.mt.common <- setdiff(rtmt.mt.common, common.deg)
rtmt.mt.common

#refined count plot function
count.plot <- function(dds_object, gene_id){
  gene_name <- gtf.dat[gene_id, "gene_name"]
  cnt.cols <- c("mosquito" = "#BCBD22", "blood" = "#E377C2", "naive" = "Black", 'recent_blood' = "#17BECF")
  cnt.dat <- plotCounts(dds_object, gene_id,
                  intgroup = c("day","delivery"), returnData = TRUE)
  # print(cnt.dat)
  cnt.trans <- cnt.dat[1:60,]
  cnt.naive <- cnt.dat[61:64, ]
  naive.mean <- mean(cnt.naive$count)
  cnt.plot <- ggplot(cnt.trans,
                     aes(x = day, y = count, color = delivery, group = delivery)) + 
    ggtitle(gene_name) + 
    geom_point() + stat_summary(fun=mean, geom="line") + 
    geom_hline(yintercept = naive.mean, col = "black", linetype="dashed") +
    scale_x_discrete(limits=1:6) + 
    scale_colour_manual(values = cnt.cols)
  return(cnt.plot)
}

for (id in rtmt.mt.common){
  gene_name <- gtf.dat[id, "gene_name"]
  pdf(paste(results.dir, "rtmt_mt_deg_profiles_2/", gene_name, ".pdf", sep = ""))
  count.plot(dds, id)
  dev.off()
}
# 


protein.df <- read.table(file = "~/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/3-documents/EKD11_protein_list.txt",
                         sep = "\t",
                         header = TRUE)
protein.df$gene_id <- sapply(protein.df$protein, function(x){gtf.dat[gtf.dat$gene_name  %in% x, "gene_id"]})
