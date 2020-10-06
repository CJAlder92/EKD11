### DESeq Analysis - SBP and RMT Only
### Counts and normalisation
### Runs DESeq2 program to create count matrix from fastq files for downstream analysis
### Author: Chris Alder

## Clear workspace and load packages used throughout analysis
source('/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/1-R/1-R_scripts/packages.R')

##############
### INPUTS ###
##############

## Directories - Set according to your own specific directory tree
nf.dir      <- "/Users/alderc/1-projects/CAMP/1-AS_timecourse/2-EKD11/1-Pipeline/1-Nextflow/"
work.dir    <- "/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/"
r.dir       <- paste(work.dir, "1-R/refactor/",sep='') #change once refactor is completed
tmp.dir     <- paste(work.dir,"tmp/",sep='')
data.dir    <- paste(work.dir,"4-data/",sep='')
results.dir <- paste(work.dir,"2-results/",sep='')
genome.dir  <- "/Users/alderc/1-projects/9-Data/1-Reference_genomes/1-Mus_musculus/"

# Checks for directories before creating (to prevent overwriting)
for (dir in grep("\\.dir$",ls(),value=T)) {
  if (!file.exists(get(dir))) { dir.create(get(dir),recursive=TRUE, mode="0755"); }
}

## Files
gtf.file <- paste(genome.dir, "Mus_musculus.GRCm38.86.gtf", sep="")
design.file <- paste(nf.dir,"design.csv",sep='')
design.baseline <- "naive" # This is the treatment group you want to use for baseline measurements and test others against - Usually control

# ## Output names
# count.output <- "genes.normalised_counts_with_controls.xls.gz" 
# log.output   <-  "genes.vst_counts_with_controls.xls.gz"
# matrix.output <- "EKD11_counts.csv"

##############
#### MAIN ####
##############

## Loading GTF file
gtf.dat    <- import(gtf.file)
gtf.dat    <- as.data.frame(gtf.dat[gtf.dat$type %in% "gene",])
gtf.dat$GRCm38 <- paste("chr",gtf.dat$seqnames,":",gtf.dat$start,"-",gtf.dat$end,sep='')
rownames(gtf.dat) <- gtf.dat$gene_id

## Design file
design.file      <- paste(nf.dir, "design.csv", sep='')
design           <- read.delim(design.file, header=TRUE, sep=",", stringsAsFactors = TRUE)
design           <- design[]
rownames(design) <- design$label
design$delivery  <- relevel(design$delivery, design.baseline) # sets design.baseline as first factor (baseline)
design$grp       <- relevel(design$grp, design.baseline) # same as above but for analysis with controls 
design$day       <- factor(as.character(design$day), levels=as.character(sort(unique(design$day)))) # factors days in numberical order
design$replicate <- as.factor(design$replicate) # factors replicates (order not important)


## Importing counts
txi <- tximport(paste(nf.dir, "results/rsem/", design$lims.name, ".genes.results", sep=""), type="rsem")
txi$length[txi$length==0] <- 1

## DESeq2
dat <- DESeqDataSetFromTximport(txi, design, ~ grp) # For baseline counts - full rank in other file
# dat <- DESeqDataSetFromTximport(txi, design, ~ delivery + day + delivery:day ) # test to see if new method works
rowData(dat) <- gtf.dat[rownames(dat),c("gene_id","gene_name","gene_source","gene_biotype","GRCm38")]
dds <- DESeq(dat)