### 
rm(list = ls())

# Libraries
library(DESeq2)
library(ComplexHeatmap)
library(tidyverse)
library(qusage)

work.dir    <- '/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/';
nf.dir      <- "/Users/alderc/1-projects/CAMP/1-AS_timecourse/2-EKD11/1-Pipeline/1-Nextflow/"
work.dir    <- "/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/"
r.dir       <- paste(work.dir, "1-R/",sep='') #change once refactor is completed
tmp.dir     <- paste(work.dir,"tmp/",sep='')
data.dir    <- paste(work.dir,"4-data/",sep='')
results.dir <- paste(work.dir,"2-results_2020/",sep='')
genome.dir  <- "/Users/alderc/1-projects/9-Data/1-Reference_genomes/1-Mus_musculus/"

 
module.anno <- read.csv(file=paste(module.dir, 'Annotation of Blood modules.csv', sep=''),
                        header=TRUE,
                        sep=',');
module.dat <- read.csv(file=paste(module.dir, 'Blood modules.csv', sep=''),
                       header=TRUE,
                       sep=',');

rlog.dat <- read.table(file= paste(results.dir, 'genes.rlog_counts.xls', sep=''),
                       header=TRUE,
                       sep='\t');

#List of Module IDs
module.num <- module.anno$X
module.anno

module.geneset <- lapply(module.num, function(mod){
  dat <- module.dat[module.dat$Module == mod, 'X']; 
  dat <- dat[!is.na(dat)]
  mod <- dat
})
names(module.geneset) <- module.num


#Check to see if Merge was correct
for (mod in names(module.geneset)){
  anno <- as.numeric(module.anno[module.anno$X == mod, 'Number.of.genes']);
  if (length(module.geneset[[mod]]) == anno){
    print(paste("Geneset", mod, 'contains correct number of genes:', anno, sep=' '))
  } else {
    print('Module gene number does not match')
  }
}

#Get Columns we need 
names(rlog.dat);
acute.cols <- grep(pattern = 'd3\\.|d6\\.|d9\\.', x = names(rlog.dat), value = TRUE)
rlog.db <- rlog.dat[ , acute.cols];
rownames(rlog.db) <- rlog.dat$gene_id;

#Import DESeq results list
load(paste(results.dir, 'as.results.list.RData', sep = ''))
res.list[grep(pattern = 'mosquito\\.._vs_naive',
     x = names(res.list),
     value = TRUE)]


#List of all DEG genes from DESeq
deg.list <- lapply(names(res.list), function(n){
  if(grepl(pattern = 'mosquito\\.._vs_naive|blood\\.._vs_naive', x = n)){
    df <- res.list[[n]]
    genes <- rownames(df[df$DEG, ])
    return(genes)
  }    
})

deg.list <- unique(unlist(deg.list)


# For Mos vs Naive
mvn.cols <- grepl(pattern = 'naive|mosquito', x= names(rlog.db));
mvn.rlog <- rlog.db[, mvn.cols];

results.list <- list()
for (d in c(3,6,9)){
  cols <- grepl(pattern = paste('mosquito\\.d', d, '\\.|naive', sep='')  , x = names(mvn.rlog));
  dat <- mvn.rlog[, cols]
  print(names(dat))
  labels <- sapply(names(dat), function(x){
    ifelse(grepl(pattern = 'mosquito', x = x), 'mt','naive')
  })
  contrast = 'mt-naive';
  results.list[[paste('mvn.d', d, sep='')]] <- qusage(dat, labels, contrast, module.geneset.deg, filter.genes = TRUE);
}

# For Blood vs Naive
bvn.cols <- grepl(pattern = 'naive|blood', x= names(rlog.db));
bvn.rlog <- rlog.db[, bvn.cols];

for (d in c(3,6,9)){
  cols <- grepl(pattern = paste('blood\\.d', d, '\\.|naive', sep='')  , x = names(bvn.rlog));
  dat <- bvn.rlog[, cols]
  print(names(dat))
  labels <- sapply(names(dat), function(x){
    ifelse(grepl(pattern = 'blood', x = x), 'sbp','naive')
  })
  contrast = 'sbp-naive';
  results.list[[paste('bvn.d', d, sep='')]] <- qusage(dat, labels, contrast, module.geneset, filter.genes = TRUE);
}


mvn_all = grepl(pattern = 'mvn', x = names(results.list));
mvn_all = results.list[mvn_all]

bvn_all = grepl(pattern = 'bvn', x = names(results.list));
bvn_all = results.list[bvn_all]


# Creating Matrix for all the Mos Qusage tables
mvn.dat <- list()
n <- 1
for (m in mvn_all){
  name <- names(mvn_all)[n];
  df <- qsTable(m, number = 41)
  df <- df[, c("pathway.name", "log.fold.change")]
  colnames(df) <- c('Module', paste(name, 'logFC', sep='_'))
  mvn.dat[[name]]  <- as.data.frame(df)
  n <- n + 1 
}
mvn.heat <- Reduce(function(x,y) merge(x = x, y = y, by = "Module", no.dups=T), 
                   mvn.dat)

# mvn.heat$Module <- module.anno$Biological.process[match(mvn.heat$Module, module.anno$X)]


mvn.heat$anno <- paste(mvn.heat$Module, module.anno$Biological.process[match(mvn.heat$Module, module.anno$X)])
# mvn.heat$anno <- module.anno$Biological.process[match(mvn.heat$Module, module.anno$X)]

mvn.heat

rownames(mvn.heat) <- mvn.heat$Module
mvn.heat <- mvn.heat[,-1]
mvn.heat <- mvn.heat[module.num, ]
rownames(mvn.heat) <- mvn.heat$anno
colnames(mvn.heat) <- c('Day 3', 'Day 6', 'Day 9', 'anno')
# mvn.heat <- mvn.heat[,-4] # Not using for the bubble plot
mvn.heat.fc <- 2^mvn.heat[,1:3]
mvn.heat.fc$anno <- mvn.heat$anno
mvn.heat.fc

bvn.dat <- list()
n <- 1
for (m in bvn_all){
  name <- names(bvn_all)[n];
  df <- qsTable(m, number = 41)
  df <- df[, c("pathway.name", "log.fold.change")]
  colnames(df) <- c('Module', paste(name, 'logFC', sep='_'))
  bvn.dat[[name]]  <- as.data.frame(df)
  n <- n + 1 
}
bvn.heat <- Reduce(function(x,y) merge(x = x, y = y, by = "Module", no.dups=T), 
                   bvn.dat)

# bvn.heat$Module <- module.anno$Biological.process[match(bvn.heat$Module, module.anno$X)]

bvn.heat$anno <- paste(bvn.heat$Module, module.anno$Biological.process[match(bvn.heat$Module, module.anno$X)])

rownames(bvn.heat) <- bvn.heat$Module
bvn.heat <- bvn.heat[,-1]
bvn.heat <- bvn.heat[module.num, ]
rownames(bvn.heat) <- bvn.heat$anno
colnames(bvn.heat) <- c('Day 3', 'Day 6', 'Day 9', 'anno')
# bvn.heat <- as.matrix(bvn.heat[, -5])
bvn.heat.fc <- 2^bvn.heat[,1:3]
bvn.heat.fc$anno <- bvn.heat$anno




library(ggplot2)
library(reshape2)
library(RColorBrewer);

# d <- ggplot(diamonds, aes(x = cut, y = clarity))
# 
# d + geom_count(aes(size = stat(prop), group = carat, colour = carat)) +
#   scale_colour_gradient2() +
#   scale_size_area(max_size = 10) 


mvn.db <- mvn.heat.fc[ , !names(mvn.heat.fc) %in% 'Acute'];
mvn.db$anno <- factor(mvn.db$anno, levels = mvn.db$anno)
names(mvn.db) <- str_replace(pattern = 'Day ', replacement =  '', string = names(mvn.db))
mvn.db <- melt(mvn.db, id.vars = 'anno')

mvn.db$prop <- apply(mvn.db, 1, function(prop){
  mod <- unlist(strsplit(as.character(prop[1]), split= ' '))[1]
  day <- (prop[2])
  geneset <- module.geneset.deg[[mod]]
  num.genes <- length(geneset)
  res.name <- grep(pattern = paste('mosquito\\.', day, '_vs_naive', sep='')  , x = names(res.list), value = TRUE);
  res <- res.list[[res.name]]
  res <- res[geneset, ]
  res <- res[res$DEG == T,]
  prop <- nrow(res) / num.genes
  
})
mvn.db
mvn.db[order(mvn.db$prop), ]

x <- ggplot(mvn.db, aes(x = variable, y = reorder(anno, desc(anno)))) + ggtitle('Mosquito Transmission') + xlab('Day') + ylab('Module')

x <- x + geom_count(aes(size = prop, group = value , colour = value), alpha=0.75) +
  labs(size = "Proportion of DEG in Module", colour = 'Enrichment score') +
  scale_color_gradientn(colours = rev(brewer.pal(11, "RdBu")), limits = c(0,2.5)) +
  scale_size_area(max_size = 9) +
  # coord_fixed(ratio=0.25) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white'),
        # panel.border = element_rect(fill = "transparent"),
        plot.background = element_blank(),
        legend.background = element_rect(fill = 'transparent'))

bvn.db <- bvn.heat.fc[ , !names(bvn.heat.fc) %in% 'Acute'];
bvn.db$anno <- factor(bvn.db$anno, levels = bvn.db$anno)
names(bvn.db) <- str_replace(pattern = 'Day ', replacement =  '', string = names(bvn.db))
bvn.db <- melt(bvn.db, id.vars = 'anno')

bvn.db$prop <- apply(bvn.db, 1, function(prop){
  mod <- unlist(strsplit(as.character(prop[1]), split= ' '))[1]
  day <- prop[2]
  geneset <- module.geneset.deg[[mod]]
  num.genes <- length(geneset)
  res.name <- grep(pattern = paste('blood\\.', day, '_vs_naive', sep='')  , x = names(res.list), value = TRUE);
  res <- res.list[[res.name]]
  res <- res[geneset, ]
  res <- res[res$DEG == T,]
  prop <- nrow(res) / num.genes
  
})
bvn.db
bvn.db[order(bvn.db$prop), ]



x.blood <- ggplot(bvn.db, aes(x = variable, y = reorder(anno, desc(anno)))) + ggtitle('SBP Infection') + xlab('Day') + ylab('Module')

x.blood <- x.blood + geom_count(aes(size = prop, group = value , colour = value), alpha=0.75) +
  labs(size = "Proportion of DEG in Module", colour = 'Enrichment Score') +
  scale_color_gradientn(colours = rev(brewer.pal(11, "RdBu")), limits = c(0,2.5)) +
  scale_size_area(max_size = 9) + 
  # coord_fixed(ratio=0.4) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white'),
        # panel.border = element_rect(fill = "transparent"),
        plot.background = element_blank(),
        legend.background = element_rect(fill = 'transparent'))
# ggsave(filename =paste(results.dir,"/figures/blood_modular_bubble.png",sep=''), height = 10, width = 8)

x.plot <- egg::ggarrange(x, x.blood, nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")
ggsave(filename =paste(results.dir,"/figures/combo_modular_bubble_new2.tiff",sep=''), height = 14, width = 16, bg = 'transparent', dpi = 300, units = 'in', device ='tiff', limitsize = TRUE)

