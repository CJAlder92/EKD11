source('/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/1-R/1-R_scripts/packages.R')

## Directories
nf.dir      <- "/Users/alderc/1-projects/CAMP/1-AS_timecourse/2-EKD11/1-Pipeline/1-Nextflow/";
work.dir    <- "/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/";
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
design      <- read.delim(design.file,header=TRUE,sep=",", stringsAsFactors = TRUE) %>% subset(delivery !='naive'); 

rownames(design) <- design$label;
#sets naive to be first level of delivery
design$delivery  <- relevel(design$delivery,"naive");
#orders days in numerical order
design$day       <- factor(as.character(design$day),levels=as.character(sort(unique(design$day))));
design$replicate <- as.factor(design$replicate)

# create a copy of column grp2 and removes .day from naive set so relevel will order all naive as first
design$grp2      <- factor(sub("naive.*","naive",as.character(design$grp)));
design$grp2      <- relevel(design$grp2,"naive");

is.factor(design$day)

## Import RSEM counts and length information
txi <- tximport(paste(nf.dir,"results/rsem/",design$lims.name,".genes.results",sep=''), type="rsem")
txi$length[txi$length==0] <- 1


## Create DESeq2 object
dat <- DESeqDataSetFromTximport(txi, design, ~ grp2)
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
write.table(count.dBase,file=gzfile(paste(results.dir,"genes.normalised_counts_without_controls.xls.gz",sep='')),col.names=TRUE,row.names=FALSE,sep="\t",quote=F);


## Rlog values
## Count data is converted to the log2 scale in a way which minimizes differences between samples for rows with small counts and which normalizes with respect to library size
rlog.dBase <- assay(rld);
rlog.dBase <- data.frame(rlog.dBase,rowData(dds),stringsAsFactors=FALSE);
write.table(rlog.dBase,file=gzfile(paste(results.dir,"genes.rlog_counts_without_controls.xls.gz",sep='')),col.names=TRUE,row.names=FALSE,sep="\t",quote=F);





## PCA plot
#part of DESeq
pca.dat  <- plotPCA(rld, intgroup = c("delivery","day","grp","label.short"),returnData=TRUE,ntop=1000);
perc.var <- round(100 * attr(pca.dat, "percentVar"));
pca.dat$day <- factor(as.character(pca.dat$day),as.character(sort(unique(pca.dat$day))));
pca.plot <- ggplot(pca.dat, aes(x=PC1, y=PC2, colour=delivery, shape=day, label=label.short)) + geom_point(size =3);
pca.plot <- pca.plot + xlab(paste0("PC1: ", perc.var[1], "% variance"));
pca.plot <- pca.plot + ylab(paste0("PC2: ", perc.var[2], "% variance"));
pca.plot <- pca.plot + geom_text(aes(label=label.short),hjust=0, vjust=0, size=3, nudge_x = 0.2, nudge_y = 0.4);
pca.cols <- c("mosquito" = "firebrick", "blood" = "darkblue", "naive" = "#E69F00", 'recent_blood' = 'springgreen')
pca.plot <- pca.plot + scale_colour_manual(values = pca.cols)
pdf(paste(results.dir,"deseq2.pca_plot_without_controls.pdf",sep=''),width=9,height=8);
print(pca.plot);
dev.off();


### All delivery pair-wise comparisons at each timepoint
res.list <- list();
design(dds) <- ~ grp2;
dds <- DESeq(dds);
resultsNames(dds);
res.list[['mt.1_vs_sbp.1']]      <- results(dds, contrast=c("grp2","mt.1","sbp.1"));
res.list[['mt.2_vs_sbp.2']]      <- results(dds, contrast=c("grp2","mt.2","sbp.2"));
res.list[['mt.3_vs_sbp.3']]      <- results(dds, contrast=c("grp2","mt.3","sbp.3"));
res.list[['mt.4_vs_sbp.4']]      <- results(dds, contrast=c("grp2","mt.4","sbp.4"));
res.list[['mt.6_vs_sbp.6']]      <- results(dds, contrast=c("grp2","mt.6","sbp.6"));

res.list[['mt.1_vs_rtmt.1']]      <- results(dds, contrast=c("grp2","mt.1","rtmt.1"));
res.list[['mt.2_vs_rtmt.2']]      <- results(dds, contrast=c("grp2","mt.2","rtmt.2"));
res.list[['mt.3_vs_rtmt.3']]      <- results(dds, contrast=c("grp2","mt.3","rtmt.3"));
res.list[['mt.4_vs_rtmt.4']]      <- results(dds, contrast=c("grp2","mt.4","rtmt.4"));
res.list[['mt.6_vs_rtmt.6']]      <- results(dds, contrast=c("grp2","mt.6","rtmt.6"));

res.list[['rtmt.1_vs_sbp.1']]      <- results(dds, contrast=c("grp2","rtmt.1","sbp.1"));
res.list[['rtmt.2_vs_sbp.2']]      <- results(dds, contrast=c("grp2","rtmt.2","sbp.2"));
res.list[['rtmt.3_vs_sbp.3']]      <- results(dds, contrast=c("grp2","rtmt.3","sbp.3"));
res.list[['rtmt.4_vs_sbp.4']]      <- results(dds, contrast=c("grp2","rtmt.4","sbp.4"));
res.list[['rtmt.6_vs_sbp.6']]      <- results(dds, contrast=c("grp2","rtmt.6","sbp.6"));



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
cbind(sapply(res.list,function(x){sum(x$DEG)}));

deg.genes <- unique(unlist(sapply(res.list, function(x){
  deg <- rownames(x[x$DEG, ])
  return(deg)})))
deg.genes


## Exporting results table

for (n in names(res.list)){
  tmp.df <- res.list[[n]]
  tmp.df <- tmp.df[tmp.df$DEG, ]
  tmp.df$gene_name <- gtf.dat[rownames(tmp.df), 'gene_name']
  write.table(tmp.df, paste(results.dir, n, '_deg_without_controls_grp.csv', sep = ''),
              row.names = TRUE,
              col.names = NA,
              sep = ',',
              quote = F)
  
}

for (gid in deg.genes){
  sym <- gtf.dat[gid,"gene_name"]
  plot.file <- paste(results.dir,"2-gene_profiles/",sym,".png",sep='');
  if (!file.exists(plot.file)) {
    cnt.dat  <- plotCounts(dds, which(rownames(dds) %in% gid), intgroup = c("day","delivery"), returnData = TRUE);
    #cnt.dat$day <- as.numeric(as.character(cnt.dat$day))
    cnt.plot <- ggplot(cnt.dat, aes(x = day, y = count, color = delivery, fill=delivery, group = delivery)) + geom_point() + geom_smooth(se = FALSE, alpha=0.1, method = "loess", span=0.75);
    cnt.plot <- cnt.plot + ggtitle(sym) + ylab("normalised count") + xlab("day");
    #cnt.plot <- cnt.plot + scale_y_log10();
    cnt.cols <- c("mosquito" = "firebrick", "blood" = "darkblue", "naive" = "#E69F00")
    cnt.plot <- cnt.plot + scale_colour_manual(values = pca.cols)
    png(file=paste(results.dir,"2-gene_profiles/",sym,".png",sep=''),height=500,width=700);
    print(cnt.plot);
    dev.off();
  }
}

mvb3.df <- mvb3.tmp[mvb3.genes, c('baseMean', 'log2FC')]
tpm <- as.data.frame(txi$abundance)
names(tpm) <- design$label

tpm.group <- function(df, treatment, day){
  genes <- rownames(df)
  cols <- grep(pattern = paste(treatment, '.d', day, sep = ''), x = names(tpm))
  dat <- tpm[genes, cols]
  return(rowMeans(dat))
}

mvb3.df$tpm.M <- tpm.group(mvb3.df, 'mosquito', 3)
mvb3.df$tpm.B <- tpm.group(mvb3.df, 'blood', 3)
mvb3.df$gene.sym <- gtf.dat[rownames(mvb3.df), 'gene_name']

mvb.list <- lapply(names(res.list[9:12]), function(n){
  df <- res.list[[n]]
  x <- as.numeric(str_extract_all(n, "[0-9]+$"))
  df <- df[df$DEG, c('baseMean', 'log2FC')]
  df$tpm.M <- tpm.group(df, 'mosquito', x)
  df$tpm.B <- tpm.group(df, 'blood', x)
  df$gene.sym <- gtf.dat[rownames(df), 'gene_name']
  df <- df[order(df$log2FC), ]
  return(df)
})
names(mvb.list) <- names(res.list[9:12])

grep(x = names(res.list[10]), pattern = '\\..$', value = TRUE)