### DESeq Analysis
### Produces gene profile graphs using ggplot of genes inputted via gene list (ENSEMBL IDs)
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
results.dir <- paste(work.dir,"2-results/",sep='')
genome.dir  <- "/Users/alderc/1-projects/9-Data/1-Reference_genomes/1-Mus_musculus/"
# Output directory
profiles.dir<- paste(results.dir, "/2-gene_profiles/test/", sep = "")

# Checks for directories before creating (to prevent overwriting)
for (dir in grep("\\.dir$",ls(),value=T)) {
  if (!file.exists(get(dir))) { dir.create(get(dir),recursive=TRUE, mode="0755"); }
}


## Files or list
design.file <- paste(nf.dir,"design.csv",sep="")
count.file  <- paste(results.dir, "EKD11_counts.csv", sep="")
ensembl.file <- paste(results.dir, "ensembl_decode.csv", sep="")

# Example gene list, insert own or file here
gene.list <- c("Ifngr2", "Tlr2", "Nrp2", "Tnfrsf4", "Nfkbia", "Cd44", "Cd86", "Tnfaip3", "Mcemp1", "Csf1", "Bst1")  # string or file path 

##############
#### MAIN ####
##############

for (gid in gene.list){
  sym <- gtf.dat[gtf.dat$gene_name %in% gid,"gene_id"]
  plot.file <- paste(results.dir,"2-gene_profiles/rmt_sbp/",gid,".png",sep='');
  if (!file.exists(plot.file)) {
    cnt.dat  <- plotCounts(dds, which(rownames(dds) %in% sym), intgroup = c("day","delivery"), returnData = TRUE);
    #cnt.dat$day <- as.numeric(as.character(cnt.dat$day))
    cnt.plot <- ggplot(cnt.dat, aes(x = day, y = count, color = delivery, fill=delivery, group = delivery)) + geom_point() + geom_smooth(se = FALSE, alpha=0.1, method = "loess", span=0.75);
    cnt.plot <- cnt.plot + ggtitle(gid) + ylab("normalised count") + xlab("day");
    #cnt.plot <- cnt.plot + scale_y_log10();
    cnt.cols <- c("mosquito" = "grey50", "blood" = "darkblue", "naive" = "#E69F00", 'recent_blood' = 'forestgreen')
    cnt.plot <- cnt.plot + scale_colour_manual(values = cnt.cols)
    png(file=paste(profiles.dir ,gid,".png",sep=''),height=500,width=700);
    print(cnt.plot);
    dev.off();
  }
}

##############
#### END #####
##############