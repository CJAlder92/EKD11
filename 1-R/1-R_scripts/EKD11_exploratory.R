source('/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/1-R/1-R_scripts/packages.R')

## Directories
nf.dir      <- "/Users/alderc/1-projects/CAMP/1-AS_timecourse/2-EKD11/1-Pipeline/1-Nextflow/";
work.dir    <- "/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/";
r.dir       <- paste(work.dir, "1-R/1-R_scripts/",sep='');
tmp.dir     <- paste(work.dir,"tmp/",sep='');
data.dir    <- paste(work.dir,"4-data/",sep='');
results.dir <- paste(work.dir,"2-results/",sep='');
# topGo.dir   <- paste(results.dir, "5-pathway_analysis/topGO/", sep='');


for (dir in grep("\\.dir$",ls(),value=T)) {
  if (!file.exists(get(dir))) { dir.create(get(dir),recursive=TRUE, mode="0755"); }
}


## Genome
genome.dir <- "~/1-projects/9-Data/1-Reference_genomes/1-Mus_musculus/";
genome     <- "GRCm38";
gtf.file   <- paste(genome.dir,"Mus_musculus.GRCm38.86.gtf",sep='');
gtf.dat    <- import(gtf.file);
gtf.dat    <- as.data.frame(gtf.dat[gtf.dat$type %in% "gene",]);
gtf.dat$GRCm38 <- paste("chr",gtf.dat$seqnames,":",gtf.dat$start,"-",gtf.dat$end,sep='');
rownames(gtf.dat) <- gtf.dat$gene_id;


## Design file
design.file <- paste(nf.dir,"design.csv",sep='');
design      <- read.delim(design.file,header=TRUE,sep=",", stringsAsFactors = TRUE);
rownames(design) <- design$label;
#sets naive to be first level of delivery
design$delivery  <- relevel(design$delivery,"naive");
#orders days in numerical order
design$day       <- factor(as.character(design$day),levels=as.character(sort(unique(design$day))));
design$replicate <- as.factor(design$replicate)

# create a copy of column grp2 and removes .day from naive set so relevel will order all naive as first
design$grp2      <- factor(sub("naive.*","naive",as.character(design$grp)));
design$grp2      <- relevel(design$grp2,"naive");

design 

## Import RSEM counts and length information
txi <- tximport(paste(nf.dir,"results/rsem/",design$lims.name,".genes.results",sep=''), type="rsem")
txi$length[txi$length==0] <- 1


## Create DESeq2 object
dat <- DESeqDataSetFromTximport(txi, design, ~ grp) # To make Full rank model - Testing with controls removed to check for changes
rowData(dat) <- gtf.dat[rownames(dat),c("gene_id","gene_name","gene_source","gene_biotype","GRCm38")];


## Normalise / scale and stabilise variance
## As we've supplied tx length info from tximport it's necessary to call normalizationFactors() rather than sizeFactors() to get the normalisation factors

# was normlizationfactors used here?!
dds <- DESeq(dat);
rld <- rlog(dat, blind = TRUE);



#Prefiltering for Qusage
keep <- rowSums(counts(dds)) >= 10;
rld.filter <- dds[keep,];
rld.filter <- rlog(rld.filter);


?rlog

## Normalised counts.
# counts default is normalizationfactors
count.norm  <- counts(dds,normalized=TRUE);
count.dBase <- data.frame(count.norm,rowData(dds),stringsAsFactors=FALSE);
rownames(count.dBase) <- count.dBase$gene_id
write.table(count.dBase,file=gzfile(paste(results.dir,"genes.normalised_counts_with_controls.xls.gz",sep='')),col.names=TRUE,row.names=FALSE,sep="\t",quote=F);


## Rlog values
## Count data is converted to the log2 scale in a way which minimizes differences between samples for rows with small counts and which normalizes with respect to library size
rlog.dBase <- assay(rld);
rlog.dBase <- data.frame(rlog.dBase,rowData(dds),stringsAsFactors=FALSE);
write.table(rlog.dBase,file=gzfile(paste(results.dir,"genes.rlog_counts_with_controls.xls.gz",sep='')),col.names=TRUE,row.names=FALSE,sep="\t",quote=F);





## PCA plot
#part of DESeq
pca.dat  <- plotPCA(rld, intgroup = c("delivery","day","grp","label.short"),returnData=TRUE,ntop=1000);
perc.var <- round(100 * attr(pca.dat, "percentVar"));
pca.dat$day <- factor(as.character(pca.dat$day),as.character(sort(unique(pca.dat$day))));
pca.plot <- ggplot(pca.dat, aes(x=PC1, y=PC2, colour=delivery, shape=day, label=label.short)) + geom_point(size =3);
pca.plot <- pca.plot + xlab(paste0("PC1: ", perc.var[1], "% variance"));
pca.plot <- pca.plot + ylab(paste0("PC2: ", perc.var[2], "% variance"));
pca.plot <- pca.plot + geom_text(aes(label=label.short),hjust=0, vjust=0, size=3, nudge_x = 0.2, nudge_y = 0.4);
pca.cols <- c("mosquito" = "firebrick", "blood" = "darkblue", "naive" = "#E69F00", 'recent_blood' = 'forestgreen')
pca.plot <- pca.plot + scale_colour_manual(values = pca.cols)
pdf(paste(results.dir,"deseq2.pca_plot.pdf",sep=''),width=9,height=8);
print(pca.plot);
dev.off();


### All delivery pair-wise comparisons at each timepoint
res.list <- list();
design(dds) <- ~ grp;
dds <- DESeq(dds);
resultsNames(dds);  
res.list[["mt.1_vs_naive"]]      <- results(dds, contrast=c("grp","mt.1","naive"));
res.list[["mt.2_vs_naive"]]      <- results(dds, contrast=c("grp","mt.2","naive"));
res.list[["mt.3_vs_naive"]]      <- results(dds, contrast=c("grp","mt.3","naive"));
res.list[["mt.4_vs_naive"]]      <- results(dds, contrast=c("grp","mt.4","naive"));
res.list[["mt.6_vs_naive"]]      <- results(dds, contrast=c("grp","mt.6","naive"));

res.list[["rtmt.1_vs_naive"]]      <- results(dds, contrast=c("grp","rtmt.1","naive"));
res.list[["rtmt.2_vs_naive"]]      <- results(dds, contrast=c("grp","rtmt.2","naive"));
res.list[["rtmt.3_vs_naive"]]      <- results(dds, contrast=c("grp","rtmt.3","naive"));
res.list[["rtmt.4_vs_naive"]]      <- results(dds, contrast=c("grp","rtmt.4","naive"));
res.list[["rtmt.6_vs_naive"]]      <- results(dds, contrast=c("grp","rtmt.6","naive"));

res.list[["sbp.1_vs_naive"]]      <- results(dds, contrast=c("grp","sbp.1","naive"));
res.list[["sbp.2_vs_naive"]]      <- results(dds, contrast=c("grp","sbp.2","naive"));
res.list[["sbp.3_vs_naive"]]      <- results(dds, contrast=c("grp","sbp.3","naive"));
res.list[["sbp.4_vs_naive"]]      <- results(dds, contrast=c("grp","sbp.4","naive"));
res.list[["sbp.6_vs_naive"]]      <- results(dds, contrast=c("grp","sbp.6","naive"));

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



## Filter results
res.list <- lapply(res.list,function(res){
  res <- as.data.frame(res);
  colnames(res)[colnames(res) %in% "log2FoldChange"] <- "log2FC";
  colnames(res)[colnames(res) %in% "padj"] <- "FDR";
  res$DEG <- res$FDR <= 0.001 &
    abs(res$log2FC)>= 1 & #only filters out logFC ±1
    res$baseMean>=30         &
    !is.na(res$log2FC)       &
    !is.na(res$FDR);
  res;
})
deg.meta <- cbind(sapply(res.list,function(x){sum(x$DEG)}));

## Exporting results table

for (n in names(res.list)){
  tmp.df <- res.list[[n]]
  tmp.df <- tmp.df[tmp.df$DEG, ]
  tmp.df$gene_name <- gtf.dat[rownames(tmp.df), 'gene_name']
  write.table(tmp.df, paste(results.dir, n, '_deg.csv', sep = ''),
              row.names = TRUE,
              col.names = NA,
              sep = ',',
              quote = F)
  
}

## DEG gene list
deg.genes <- unique(unlist(sapply(res.list, function(x){
  deg <- rownames(x[x$DEG, ])
  return(deg)})))
deg.genes

## Individual expression profiles for DEG genes 
for (gid in deg.genes){
  sym <- gtf.dat[gid,"gene_name"]
  plot.file <- paste(results.dir,"2-gene_profiles/",sym,".png",sep='');
  if (!file.exists(plot.file)) {
    cnt.dat  <- plotCounts(dds, which(rownames(dds) %in% gid), intgroup = c("day","delivery"), returnData = TRUE);
    #cnt.dat$day <- as.numeric(as.character(cnt.dat$day))
    cnt.plot <- ggplot(cnt.dat, aes(x = day, y = count, color = delivery, fill=delivery, group = delivery)) + geom_point() + geom_smooth(se = FALSE, alpha=0.1, method = "loess", span=0.75);
    cnt.plot <- cnt.plot + ggtitle("ENSMUSG00000022965") + ylab("normalised count") + xlab("day");
    #cnt.plot <- cnt.plot + scale_y_log10();
    cnt.cols <- c("mosquito" = "firebrick", "blood" = "darkblue", "naive" = "#E69F00", 'recent_blood' = 'forestgreen')
    cnt.plot <- cnt.plot + scale_colour_manual(values = cnt.cols)
    png(file=paste(results.dir,"2-gene_profiles/",sym,".png",sep=''),height=500,width=700);
    print(cnt.plot);
    dev.off();
  }
}


## Checking sample distances for outliers

sampledist <- dist(t(assay(vsd)))
sampledist.matrix <- as.matrix(sampledist)
ht <- Heatmap(sampledist.matrix,
              name = 'sample_cor',
              cluster_rows= FALSE,
              cluster_columns = FALSE,
              clustering_distance_rows = 'pearson',
              clustering_distance_columns = 'pearson')
              # col = viridis::viridis(100))
pdf(paste(results.dir,"sample_cor_matrix.pdf",sep=''),width=15,height=15);
draw(ht)
decorate_heatmap_body('sample_cor', c(
  for (d in seq(0, ncol(sampledist.matrix), 4)){
    grid.lines(x=c(d/64, d/64), y=c(0,1), gp= gpar(col = 'black', lty='solid', lwd=2))
    grid.lines(y=c(d/64, d/64), x=c(0,1), gp= gpar(col = 'black', lty='solid', lwd=2))}))
dev.off()

# Cooks distance
cooks.mat <- log10(assays(dds)[['cooks']])
pdf(paste(results.dir,'cooks_distance.pdf', sep =''), width=15,height=10)
boxplot(cooks.mat, las = 2, range = 0,
        main = 'Sample distribution cooks distance',
        ylab = 'log10(Cooks distance)')
title(xlab = 'Sample', line = 4 )
dev.off()
cooks.melt <- melt(cooks.mat)
cooks.melt

# Dendogram
c <- cor(assay(rld), method = 'pearson')
c.dist <- as.dist(1-c)
dist.mat <- dist(t(assay(rld)), method = "euclidean")
head(dist.mat)
hr.cor <- hclust(c.dist, method = "complete", members=NULL)
hr.dist <- hclust(dist.mat, method = 'complete', members= NULL)
plot(hr.dist)

## Comparison of old dataset (Lin AS) with new (EKD11) to see why fewer DEG genes
Lin.path <- '~/1-projects/1-PIR_project/1-AS_Acute_Chronic_D3-60/3-results/6-deg/'
LinAS.file.list <- list.files(path = Lin.path,
                              pattern = '.tsv$')
LinAS.file.list

lin.res <- lapply(LinAS.file.list, function(f){
  df <- read.table(paste(Lin.path, f, sep = ''), sep = '\t', header = TRUE)
  return(df$gene_id)
})
names(lin.res) <- c('bvn.day3', 'bvn.day6', 'mvb.day3', 'mvb.day6', 'mvn.day3', 'mvn.day6')

res.names <- c('sbp.3_vs_naive', 'sbp.6_vs_naive', 'mt.3_vs_sbp.3', 'mt.6_vs_sbp.6', 'mt.3_vs_naive', 'mt.6_vs_naive')

# List withv all DEG identified in Lin AS experiment
lin.deg.list <- lapply(1:6, function(i){
  gene.list <- lin.res[[i]]
  df <- res.list[[res.names[i]]]
  df <- df[gene.list, ]
  return(df)
 })
names(lin.deg.list) <- res.names
lin.deg.list[[5]]

# Getting count data from Lin experiment
lin.counts <- read.table('~/1-projects/1-PIR_project/1-AS_Acute_Chronic_D3-60/3-results/genes.normalised_counts_test.acute.xls.gz',
                         header = TRUE,
                         row.names = 'gene_id')
lin.counts <- lin.counts[,1:24]
ekd11.counts <- count.dBase[, 1:64]
cnames.merged <- c(names(lin.counts), names(ekd11.counts))

merged.counts <- list()

for (d in c(3,6)){
  mvn.genes <- lin.res[[grep(pattern = paste('mvn.day', d, sep=''), x = names(lin.res))]]
  mvn.cols <- grep(pattern = paste('^mt.d', d,'|mosquito.d', d, '|naive', sep=''), x = cnames.merged, value = TRUE)
  mvn.lin <- lin.counts[mvn.genes, colnames(lin.counts) %in% mvn.cols]
  mvn.ekd11 <- ekd11.counts[mvn.genes, colnames(ekd11.counts) %in% mvn.cols]
  mvn.merged <- merge(mvn.lin, mvn.ekd11, by = 0)
  merged.counts[[paste('mvn.day', d, sep='')]] <- mvn.merged
  bvn.genes <- lin.res[[grep(pattern = paste('bvn.day', d, sep=''), x = names(lin.res))]]
  bvn.cols <- grep(pattern = paste('sbp.d', d,'|blood.d', d, '|naive', sep=''), x = cnames.merged, value = TRUE)
  bvn.lin <- lin.counts[bvn.genes, colnames(lin.counts) %in% bvn.cols]
  bvn.ekd11 <- ekd11.counts[bvn.genes, colnames(ekd11.counts) %in% bvn.cols]
  bvn.merged <- merge(bvn.lin, bvn.ekd11, by = 0)
  merged.counts[[paste('bvn.day', d, sep='')]] <- bvn.merged
  mvb.genes <- lin.res[[grep(pattern = paste('mvb.day', d, sep=''), x = names(lin.res))]]
  mvb.cols <- grep(pattern = paste('^mt.d', d,'|mosquito.d', d, '|sbp.d', d, '|blood.d', d, sep=''), x = cnames.merged, value = TRUE)
  mvb.lin <- lin.counts[mvn.genes, colnames(lin.counts) %in% mvb.cols]
  mvb.ekd11 <- ekd11.counts[mvb.genes, colnames(ekd11.counts) %in% mvb.cols]
  mvb.merged <- merge(mvb.lin, mvb.ekd11, by = 0)
  merged.counts[[paste('mvb.day', d, sep='')]] <- mvb.merged
}


# PCA for transmission groups
for (t in c('mosquito', 'blood', 'recent_blood')){
  rld.sub <- rld[ ,rld$delivery %in% t]
  pca.dat  <- plotPCA(rld.sub , intgroup = c("delivery","day","grp","label.short"),returnData=TRUE,ntop=1000);
  perc.var <- round(100 * attr(pca.dat, "percentVar"));
  pca.dat$day <- factor(as.character(pca.dat$day),as.character(sort(unique(pca.dat$day))));
  pca.plot <- ggplot(pca.dat, aes(x=PC1, y=PC2, colour=delivery, shape=day, label=label.short)) + geom_point(size =3);
  pca.plot <- pca.plot + xlab(paste0("PC1: ", perc.var[1], "% variance"));
  pca.plot <- pca.plot + ylab(paste0("PC2: ", perc.var[2], "% variance"));
  pca.plot <- pca.plot + geom_text(aes(label=label.short),hjust=0, vjust=0, size=3, nudge_x = 0.2, nudge_y = 0.4);
  pca.cols <- c("mosquito" = "firebrick", "blood" = "darkblue", "naive" = "#E69F00", 'recent_blood' = 'forestgreen')
  pca.plot <- pca.plot + scale_colour_manual(values = pca.cols)
  pdf(paste(results.dir, t, ".pca_plot.pdf",sep=''),width=9,height=8);
  print(pca.plot)
  dev.off()
}
 # PCA for day
for (d in c(1,2,3,4,6)){
  rld.sub <- rld[ ,rld$day %in% d]
  pca.dat  <- plotPCA(rld.sub , intgroup = c("delivery","day","grp","label.short"),returnData=TRUE,ntop=1000);
  perc.var <- round(100 * attr(pca.dat, "percentVar"));
  pca.dat$day <- factor(as.character(pca.dat$day),as.character(sort(unique(pca.dat$day))));
  pca.plot <- ggplot(pca.dat, aes(x=PC1, y=PC2, colour=delivery, shape=day, label=label.short)) + geom_point(size =3);
  pca.plot <- pca.plot + xlab(paste0("PC1: ", perc.var[1], "% variance"));
  pca.plot <- pca.plot + ylab(paste0("PC2: ", perc.var[2], "% variance"));
  pca.plot <- pca.plot + geom_text(aes(label=label.short),hjust=0, vjust=0, size=3, nudge_x = 0.2, nudge_y = 0.4);
  pca.cols <- c("mosquito" = "firebrick", "blood" = "darkblue", "naive" = "#E69F00", 'recent_blood' = 'forestgreen')
  pca.plot <- pca.plot + scale_colour_manual(values = pca.cols)
  pdf(paste(results.dir, "day",d, ".pca_plot.pdf",sep=''),width=9,height=8);
  print(pca.plot)
  dev.off()
}


## Analysis 
rt_vs_sbp <- res.list[grep('rtmt.._vs_sbp..', x = names(res.list), value = TRUE)]

rt_vs_sbp <- lapply(rt_vs_sbp, function(df){
  df <- df[df$DEG,]
  df$gene_name <- gtf.dat[rownames(df), 'gene_name']
  return(df)
})
rt_vs_sbp

### rtmt pathway 

for (d in c('mt', 'rtmt', 'sbp')){
  res <- list()
  res[["1.vs.2"]] <- results(dds, contrast=c("grp",paste(d,".1",sep=""), paste(d,".2",sep="")));
  res[["2.vs.3"]] <- results(dds, contrast=c("grp",paste(d,".2",sep=""), paste(d,".3",sep="")));
  res[["3.vs.4"]] <- results(dds, contrast=c("grp",paste(d,".3",sep=""), paste(d,".4",sep="")));
  res[["4.vs.6"]] <- results(dds, contrast=c("grp",paste(d,".4",sep=""), paste(d,".6",sep="")));
  
  
  res <- lapply(res,function(df){
    df <- as.data.frame(df);
    colnames(df)[colnames(df) %in% "log2FoldChange"] <- "log2FC";
    colnames(df)[colnames(df) %in% "padj"] <- "FDR";
    df$DEG <- df$FDR <= 0.01 &
      abs(df$log2FC)>= 1 & #only filters out logFC ±1
      df$baseMean>=30         &
      !is.na(df$log2FC)       &
      !is.na(df$FDR);
    df;
  })
    deg <- lapply(res, function(df){
      df <- df[df$DEG, ]
      df$gene_name <- gtf.dat[rownames(df), 'gene_name']
      return(df)
    })

    assign(paste(d, '.res',sep=''),res)
    assign(paste(d, '.deg',sep=''),deg)
  
}

##### COMPARING Transmission routes vs Naive #####

# Mosquito vs Naive results and DEG
mvn.res <- res.list[grep('^mt.._vs_naive', names(res.list))]
mvn.deg <- lapply(mvn.res, function(df){return(df[df$DEG,])})

# Recently MT vs Naive results and DEG
rtvn.res <- res.list[grep('^rtmt.._vs_naive', names(res.list))]
rtvn.deg <- lapply(rtvn.res, function(df){return(df[df$DEG,])})

# SBP vs Naive results and DEG 
svn.res <- res.list[grep('^sbp.._vs_naive', names(res.list))]
svn.deg <- lapply(svn.res, function(df){return(df[df$DEG,])})

# Genes involved in each transmission route
mt.genes <- unique(unlist(sapply(mvn.deg, row.names)))
rtmt.genes <- unique(unlist(sapply(rtvn.deg, row.names))) 
sbp.genes <- unique(unlist(sapply(svn.deg, row.names))) 

# Common genes involved in each transmission response
common.genes <- Reduce(intersect, list(mt.genes, rtmt.genes, sbp.genes))

# Genes shared by MT and RTMT, not used in SBP 
mt.rtmt.genes <- intersect(mt.genes, rtmt.genes)
mt.rtmt.genes <- setdiff(mt.rtmt.genes, common.genes)
# length(mt.rtmt.genes)

# Genes shared by MT and SBP, not used in RTMT 
mt.sbp.genes <- intersect(mt.genes, sbp.genes)
mt.sbp.genes <- setdiff(mt.sbp.genes, common.genes)
# length(mt.sbp.genes)

# Genes shared by RTMT and SBP, not used in MT
rtmt.sbp.genes <- intersect(rtmt.genes, sbp.genes)
rtmt.sbp.genes <- setdiff(rtmt.sbp.genes, common.genes)
# length(rtmt.sbp.genes)


gtf.dat[gtf.dat$gene_name %in% 'Ifng', 'gene_id']
#Gene Expression of individual genes
cnt.dat  <- plotCounts(dds, which(rownames(dds) %in% "ENSMUSG00000055170"), intgroup = c("day","delivery"), returnData = TRUE);
#cnt.dat$day <- as.numeric(as.character(cnt.dat$day))
cnt.plot <- ggplot(cnt.dat, aes(x = day, y = count, color = delivery, fill=delivery, group = delivery)) + geom_point() + geom_smooth(se = FALSE, alpha=0.1);
cnt.plot <- cnt.plot + ggtitle('Ifng') + ylab("normalised count") + xlab("day");
#cnt.plot <- cnt.plot + scale_y_log10();
cnt.cols <- c("mosquito" = "firebrick", "blood" = "darkblue", "naive" = "#E69F00", 'recent_blood' = 'forestgreen')
cnt.plot <- cnt.plot + scale_colour_manual(values = pca.cols)
ggplotly(cnt.plot)




# Heatmap of RTMT + MT genes genes

heat.genes <- mt.rtmt.genes

heat.dat   <- as.matrix(assay(rld)[heat.genes,]);
heat.dat   <- t(apply(heat.dat,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))}));
rownames(heat.dat) <- gtf.dat[heat.genes,"gene_name"];

ht <- Heatmap(heat.dat,
              name = "Z-Score",
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
              # split                   = split.vec, #This would use vector I created above to split and cluster the genes
              # top_annotation          = ha1,
              # bottom_annotation       = anno_bottom,
              col                     = viridis::viridis(100));

ht
