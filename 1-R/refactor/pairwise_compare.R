### DESeq Analysis
### Runs pairwise comparisons of conditions, which can then be used to determine differentially expressed genes
### Author: Chris Alder

## Load packages used throughout analysis
rm(list = ls())
source('/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/1-R/refactor/packages.R')

##############
### INPUTS ###
##############

## Directories - Set according to your own specific directory tree
nf.dir      <- "/Users/alderc/1-projects/CAMP/1-AS_timecourse/2-EKD11/1-Pipeline/1-Nextflow/"
work.dir    <- "/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/"
r.dir       <- paste(work.dir, "1-R/refactor/",sep='') #change once refactor is completed
tmp.dir     <- paste(work.dir,"tmp/",sep='')
data.dir    <- paste(work.dir,"4-data/",sep='')
results.dir <- paste(work.dir,"2-results_2020/",sep='') # For 2020 Analysis

# Checks for directories before creating (to prevent overwriting)
for (dir in grep("\\.dir$",ls(),value=T)) {
  if (!file.exists(get(dir))) { dir.create(get(dir),recursive=TRUE, mode="0755"); }
}

## Files
design.file <- paste(nf.dir,"design_2020.csv",sep="")
count.file  <- paste(results.dir, "EKD11_counts.csv", sep="")
ensembl.file <- paste(results.dir, "ensembl_decode.csv", sep="")

design.baseline <- "naive"

##############
#### MAIN ####
##############

## Ensembl file
ensembl.decode <- read.delim(ensembl.file, header=TRUE, sep=",")
rownames(ensembl.decode) <- ensembl.decode$gene_id

## Design file
design           <- read.delim(design.file, header=TRUE, sep=",", stringsAsFactors = TRUE)
rownames(design) <- design$label
design$delivery  <- relevel(design$delivery, design.baseline) # sets design.baseline as first factor (baseline)
design$grp       <- relevel(design$grp, design.baseline) # same as above but for analysis with controls 
design$day       <- factor(as.character(design$day), levels=as.character(sort(unique(design$day)))) # factors days in numberical order
design$replicate <- as.factor(design$replicate) # factors replicates (order not important)

## Count file
count.mat <- read.delim(count.file, header=TRUE, sep=",")
all(colnames(count.mat) == rownames(design)) # Check if the order of the count matrix matches the design dataframe

## Creating dds from count matrix
dat <- DESeqDataSetFromMatrix(countData = count.mat, colData = design, design = ~ grp )
dds <- DESeq(dat)

## Pairwise comparisons

# Results list
res.list <- list()
design(dds) <- ~ grp
# design(dds) <- formula(~ delivery + day + delivery:day)
dds <- DESeq(dds)
resultsNames(dds)

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

# Filter results
res.list <- lapply(res.list,function(res){
  res <- as.data.frame(res)
  colnames(res)[colnames(res) %in% "log2FoldChange"] <- "log2FC"
  colnames(res)[colnames(res) %in% "padj"] <- "FDR"
  res$DEG <- res$FDR <= 0.01 &
    abs(res$log2FC)>= 0.5 & #only filters out logFC ±1 N.B CHANGE FOR RMT vs SBP COMPARISONS
    res$baseMean>=30         &
    !is.na(res$log2FC)       &
    !is.na(res$FDR)
  res
})

res.list <- lapply(res.list,function(res){
  res <- as.data.frame(res)
  colnames(res)[colnames(res) %in% "log2FoldChange"] <- "log2FC"
  colnames(res)[colnames(res) %in% "padj"] <- "FDR"
  res$DEG2 <- res$FDR <= 0.01 &
    # abs(res$log2FC)>= 0.5 & #only filters out logFC ±1 N.B CHANGE FOR RMT vs SBP COMPARISONS
    res$baseMean>=30         &
    !is.na(res$log2FC)       &
    !is.na(res$FDR)
  res
})
names(deg.list)
deg.list[[15]]
data.frame(DEG = cbind(sapply(res.list,function(x){sum(x$DEG)})),
           DEG2 = cbind(sapply(res.list,function(x){sum(x$DEG2)})))

# test <- results(dds, lfcThreshold=0.5, altHypothesis="greater", contrast=c("grp","rtmt.6","sbp.6"))

deg.list <- lapply(res.list, function(df){
  df <- df[df$DEG2, ]
  names <- colnames(df)
  df$gene_name <- ensembl.decode[rownames(df), "gene_name"]
  df[ ,c("gene_name", names)]
})

# Write DEG tables

for (n in names(deg.list)){
  df <- deg.list[[n]] %>% rownames_to_column("gene_id")
  write.table(df,
              file = paste(results.dir, "1-deg/", n, ".csv", sep=""),
              col.names = TRUE,
              row.names = FALSE,
              sep = ",",
              quote = FALSE)
}

length(unique(unlist(sapply(deg.list, row.names))))
save(deg.list, file = paste(r.dir, "objects/deg_list.rds", sep = ""))

##############
#### END #####
##############
