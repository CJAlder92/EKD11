rm(list = ls())

library(tidyverse)
library(mgcv)
library(reshape2)
library(DESeq2)

nf.dir      <- "/Users/alderc/1-projects/CAMP/1-AS_timecourse/2-EKD11/1-Pipeline/1-Nextflow/"
work.dir    <- "/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/"
r.dir       <- paste(work.dir, "1-R/",sep='') #change once refactor is completed
tmp.dir     <- paste(work.dir,"tmp/",sep='')
data.dir    <- paste(work.dir,"4-data/",sep='')
results.dir <- paste(work.dir,"2-results_2020/",sep='')
genome.dir  <- "/Users/alderc/1-projects/9-Data/1-Reference_genomes/1-Mus_musculus/"

load(paste(r.dir, "Projects/deseq.Rdata", sep = ""))
load(paste(r.dir, "Projects/results.Rdata", sep = ""))

######################################################################
######################################################################



count.db <- counts(dds, normalized = TRUE)
# vsd.db <- assay(vsd) #log counts


#Mosqutio DEG gene list
mvn <- grep(pattern = '^mt.\\d+_vs_naive', names(res.list), value= TRUE)
mos.deg <- deg.list[mvn]
mvn.gene_list <- unique(unlist(sapply(mos.deg, rownames)))

rt.mt.df<- read.table(paste(results.dir, "rt_mt_df.csv", sep = ""), sep = ",", header = TRUE)
rt.mt.ids <- rt.mt.df$gene_id
# For testing 


.db <- count.db[, grepl('^mt\\.d\\d+\\.r\\d|naive', colnames(count.db))]
mos.melt <- melt(as.matrix(mos.db))
colnames(mos.melt) <- c('gene_id', 'sample', 'count')
mos.melt <- mos.melt[mos.melt$gene_id %in% rt.mt.ids, ]
mos.melt$day <- sapply(mos.melt$sample, function(x){
  ifelse(grepl('naive', x), 
         return(0),
         return(as.numeric(str_match(x, '\\.d(.+)\\.')[,2])))
}
)




mos.curves <- gam(count ~ s(day,
                            by = gene_id,
                            bs = 'tp',
                            k = 6),
                  sp = rep(0.01, length(unique(mos.melt$gene_id))),
                  data = mos.melt,
                  family = "nb")

save(mos.curves, file = paste(r.dir, "mos.gam.rds", sep = ""))

#recent mosqutio DEG gene list
rvn <- grep(pattern = 'rtmt.\\d+_vs_naive', names(res.list), value= TRUE)
rtmt.deg <- deg.list[rvn]
rvn.gene_list <- unique(unlist(sapply(rtmt.deg, rownames)))

# rt.mt.df<- read.table(paste(results.dir, "rt_mt_df.csv", sep = ""), sep = ",", header = TRUE)
# rt.mt.ids <- rt.mt.df$gene_id
# For testing 


rtmt.db <- count.db[, grepl('mt\\.d\\d+\\.r\\d|naive', colnames(count.db))]
rtmt.melt <- melt(as.matrix(rtmt.db))
colnames(rtmt.melt) <- c('gene_id', 'sample', 'count')
rtmt.melt <- rtmt.melt[rtmt.melt$gene_id %in% rt.mt.ids, ]
rtmt.melt$day <- sapply(rtmt.melt$sample, function(x){
  ifelse(grepl('naive', x), 
         return(0),
         return(as.numeric(str_match(x, '\\.d(.+)\\.')[,2])))
}
)

rtmt.curves <- gam(count ~ s(day,
                            by = gene_id,
                            bs = 'tp',
                            k = 6),
                  sp = rep(0.01, length(unique(rtmt.melt$gene_id))),
                  data = rtmt.melt,
                  family = "nb")

###################################################

#recent mosqutio DEG gene list
svn <- grep(pattern = 'sbp.\\d+_vs_naive', names(res.list), value= TRUE)
sbp.deg <- deg.list[svn]
svn.gene_list <- unique(unlist(sapply(sbp.deg, rownames)))

# rt.mt.df<- read.table(paste(results.dir, "rt_mt_df.csv", sep = ""), sep = ",", header = TRUE)
# rt.mt.ids <- rt.mt.df$gene_id
# For testing 


sbp.db <- count.db[, grepl('sbp\\.d\\d+\\.r\\d|naive', colnames(count.db))]
sbp.melt <- melt(as.matrix(sbp.db))
colnames(sbp.melt) <- c('gene_id', 'sample', 'count')
sbp.melt <- sbp.melt[sbp.melt$gene_id %in% rt.mt.ids, ]
sbp.melt$day <- sapply(sbp.melt$sample, function(x){
  ifelse(grepl('naive', x), 
         return(0),
         return(as.numeric(str_match(x, '\\.d(.+)\\.')[,2])))
}
)

sbp.curves <- gam(count ~ s(day,
                             by = gene_id,
                             bs = 'tp',
                             k = 6),
                   sp = rep(0.01, length(unique(sbp.melt$gene_id))),
                   data = sbp.melt,
                   family = "nb")

save(sbp.curves, file = paste(r.dir, "3-RData/sbp.gam.rds", sep = ""))


