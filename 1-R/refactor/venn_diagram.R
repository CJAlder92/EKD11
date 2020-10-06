### DESeq Analysis
### Creates Venn Diagrams using DE genes found within each treatment group (compared with Naive)
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
results.dir <- paste(work.dir,"2-results_2020/",sep='')
genome.dir  <- "/Users/alderc/1-projects/9-Data/1-Reference_genomes/1-Mus_musculus/"

## Files - Old data
# deg.object <- paste(r.dir, "objects/deg_list.rds", sep = "")
# 
# 
# load(deg.object)

# From Paper 2020 data
deg.list <- lapply(res.list, function(x){
  x <- x[x$DEG, ]
})

##############
#### MAIN ####
##############

## Load deg list R objects
names(deg.list)

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

venn.cols <- c("grey50", "deeppink2", "dodgerblue2")
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

pdf(file = paste(results.dir, "venn/All.pdf", sep = ""))
grid.newpage()
grid.draw(venn)
dev.off()

# Individual days

for (i in 1:5){
  m <- rownames(mvn.deg[[i]])
  s <- rownames(svn.deg[[i]])
  r <- rownames(rvn.deg[[i]])
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
                       
                       
                       
                       # Set names
                       cat.fontface = "bold",
                       cat.default.pos = "outer",
                       cat.pos = c(-27, 27, 135),
                       cat.dist = c(0.055, 0.055, 0.085),
                       cat.fontfamily = "sans",
                       rotation = 1);
  pdf(file = paste(results.dir, "venn/Day_", x, ".pdf", sep = ""))
  grid.newpage()
  grid.draw(venn)
  dev.off()
}

## Up and Down regulated (overview)
# Up 
m.up <- unique(unlist(sapply(mvn.deg, function(df){
  genes <- rownames(df[df$log2FC > 0, ])
})))
s.up <- unique(unlist(sapply(svn.deg, function(df){
  genes <- rownames(df[df$log2FC > 0, ])
})))
r.up <- unique(unlist(sapply(rvn.deg, function(df){
  genes <- rownames(df[df$log2FC > 0, ])
})))

venn.cols <- c("firebrick", "darkblue", "forestgreen")
venn <- venn.diagram(x = list(m.up, s.up, r.up),
                     main = paste("Upregulated genes overview", sep = " "),
                     main.fontface = "bold",
                     main.fontfamily = "sans",
                     category.names = c( "M_UP", "S_UP" ,"R_UP"),
                     filename = NULL,
                     
                     #Circles
                     lwd = 2,
                     lty = "blank",
                     fill = venn.cols,
                     
                     # Numbers
                     fontface = "bold",
                     fontfamily = "sans",
                     
                     # Set names
                     cat.fontface = "bold",
                     cat.default.pos = "outer",
                     cat.pos = c(-27, 27, 135),
                     cat.dist = c(0.055, 0.055, 0.085),
                     cat.fontfamily = "sans",
                     rotation = 1)


pdf(file = paste(results.dir, "3-venn/Up_all.pdf", sep = ""))
grid.newpage()
grid.draw(venn)
dev.off()

#Down 
m.down <- unique(unlist(sapply(mvn.deg, function(df){
  genes <- rownames(df[df$log2FC < 0, ])
})))
s.down <- unique(unlist(sapply(svn.deg, function(df){
  genes <- rownames(df[df$log2FC < 0, ])
})))
r.down <- unique(unlist(sapply(rvn.deg, function(df){
  genes <- rownames(df[df$log2FC < 0, ])
})))

venn.cols <- c("firebrick", "darkblue", "forestgreen")
venn <- venn.diagram(x = list(m.down, s.down, r.down),
                     main = paste("Downregulated genes overview", sep = " "),
                     main.fontface = "bold",
                     main.fontfamily = "sans",
                     category.names = c( "M_DOWN", "S_DOWN" ,"R_DOWN"),
                     filename = NULL,
                     
                     #Circles
                     lwd = 2,
                     lty = "blank",
                     fill = venn.cols,
                     
                     # Numbers
                     fontface = "bold",
                     fontfamily = "sans",
                     
                     # Set names
                     cat.fontface = "bold",
                     cat.default.pos = "outer",
                     cat.pos = c(-27, 27, 135),
                     cat.dist = c(0.055, 0.055, 0.085),
                     cat.fontfamily = "sans",
                     rotation = 1)


pdf(file = paste(results.dir, "3-venn/Down_all.pdf", sep = ""))
grid.newpage()
grid.draw(venn)
dev.off()


## Up and Down Reg (Put together?)
for (i in 1:5){
  m <- mvn.deg[[i]]
  s <- svn.deg[[i]]
  r <- rvn.deg[[i]]
  m.up <- rownames(m[m$log2FC > 0, ])
  m.down <- rownames(m[m$log2FC < 0, ])
  s.up <- rownames(s[s$log2FC > 0, ])
  s.down <- rownames(s[s$log2FC < 0, ])
  r.up <- rownames(r[r$log2FC > 0, ])
  r.down <- rownames(r[r$log2FC < 0, ])
  ifelse(i == 5, x <- 6, x <- i)
  venn <- venn.diagram(x = list(m.up, s.up, r.up),
                       main = paste("Day", x, "Upregulated genes", sep = " "),
                       main.fontface = "bold",
                       main.fontfamily = "sans",
                       category.names = c( "M_UP", "S_UP" ,"R_UP"),
                       filename = NULL,
                       
                       #Circles
                       lwd = 2,
                       lty = "blank",
                       fill = venn.cols,
                       
                       # Numbers
                       fontface = "bold",
                       fontfamily = "sans",
                       
                       # Set names
                       cat.fontface = "bold",
                       cat.default.pos = "outer",
                       cat.pos = c(-27, 27, 135),
                       cat.dist = c(0.055, 0.055, 0.085),
                       cat.fontfamily = "sans",
                       rotation = 1)
  

  pdf(file = paste(results.dir, "venn/Day_", x, "_UP.pdf", sep = ""))
  grid.newpage()
  grid.draw(venn)
  dev.off()
}

## Unique genes to RMT and MT

rt_mt.unique <- intersect(mvn.genes, rvn.genes)
rt_mt.unique <- setdiff(rt_mt.unique, common.genes) # Remove genes shared by all, and with SBP

day4.rmt <- rvn.deg[[4]] 
day4.rmt <- day4.rmt[rownames(day4.rmt) %in% rt_mt.unique, ]
day6.rmt <- rvn.deg[[5]]
day6.rmt <- day6.rmt[rownames(day6.rmt) %in% rt_mt.unique, ]


day6.rmt.up <- day6.rmt[day6.rmt$log2FC > 0, c("gene_name", "log2FC")]
day6.rmt.down <- day6.rmt[day6.rmt$log2FC < 0, c("gene_name", "log2FC") ]


day4.mt <- mvn.deg[[4]] 
day4.mt <- day4.mt[rownames(day4.mt) %in% rt_mt.unique, ]
day6.mt <- mvn.deg[[5]]
day6.mt <- day6.mt[rownames(day6.mt) %in% rt_mt.unique, ]


day6.mt.up <- day6.mt[day6.mt$log2FC > 0, c("gene_name", "log2FC")]
day6.mt.down <- day6.mt[day6.mt$log2FC < 0, c("gene_name", "log2FC") ]

rt_mt.df <- data.frame(matrix(nrow = length(rt_mt.unique), ncol = 10))
rownames(rt_mt.df) <- rt_mt.unique
treat <- c("MT.", "RMT.")
day <- c(1,2,3,4,6)
colnames(rt_mt.df) <- as.list(outer(treat, day, paste0))
for(g in rt_mt.unique){
  for (i in 1:5){
    m <- mvn.deg[[i]]
    r <- rvn.deg[[i]]
    y.m <- i * 2 - 1
    y.r <- i * 2
    if(g %in% rownames(m)){
      ifelse(m[g, "log2FC"] > 0, x.m <- "1", x.m <- "-1")
    } else { x.m <- 0 }
    rt_mt.df[g, y.m] <- x.m
    if(g %in% rownames(r)){
      ifelse(r[g, "log2FC"] > 0, x.r <- "1", x.r <- "-1")
    } else { x.r <- 0 }
    rt_mt.df[g, y.r] <- x.r
    }
}

View(rt_mt.df)
rt_mt.vec <- apply(rt_mt.df,1,function(x){paste(x,collapse='.')})
rt_mt.vec
rt_mt.levels <- rev(sort(unique(rt_mt.vec)))
rt_mt.names  <- LETTERS[1:length(rt_mt.levels)];
names(rt_mt.names) <- rt_mt.levels;
rt_mt.vec <- as.character(rt_mt.names[rt_mt.vec]);

#Matrix
rt_mt.mat <- as.matrix(rt_mt.df)


## Heatmap 
rownames(rt_mt.mat) <- gtf.dat[rownames(rt_mt.mat), "gene_name"]
ht <- Heatmap(rt_mt.mat,
              name = "DE",
              show_row_names          = TRUE,
              cluster_columns         = FALSE,
              rect_gp                 = gpar(col = "white", lwd = 1),
              row_dend_width          = unit(30, "mm"),
              column_dend_height      = unit(30, "mm"),
              row_names_gp            = gpar(fontsize = 6),
              column_names_gp         = gpar(fontsize = 10),
              split                   = rt_mt.vec,
              col                     = c("darkblue","lightgrey","firebrick"))
ht
?Heatmap

groups <- c(0,2,4,6,8,10)
pdf(paste(results.dir,"rt_mt_unique_HM.pdf",sep=''),width=12,height=70);
draw(ht)
for (n in 1:length(rt_mt.levels)){
  decorate_heatmap_body('DE', {
    for (d in 1:10){
      ifelse(d %in% groups, x <- 'solid', x <- 'dashed')
      grid.lines(x=c(d/10, d/10), y=c(0,1), gp= gpar(col = 'black', lty=x, lwd=2))}
  }, slice=n)};
dev.off()


# Write dataframe to csv
rt_mt.df$group <- rt_mt.vec
rt_mt.df$name <- gtf.dat[rownames(rt_mt.df), "gene_name"]

rt_mt.df <- rt_mt.df %>% 
  rownames_to_column('id') %>% 
  select(id, name, everything()) %>% 
  arrange(group)

write.table(x = rt_mt.df,
            file = paste(results.dir, "rt_mt_unique_df.csv", sep =""),
            sep = ",",
            quote = FALSE,
            row.names =  FALSE
            )

# Day 6
m.6 <- mvn.deg[[5]]
s.6 <- svn.deg[[5]]
r.6 <- rvn.deg[[5]]


day6.common.genes <- Reduce(intersect, list(rownames(m.6), rownames(s.6), rownames(r.6)))
day6.common.genes

r.m.6 <- intersect(rownames(m.6), rownames(r.6)) 
r.m.6 <- setdiff(r.m.6, day6.common.genes)
setdiff(r.m.6, rownames(s.6))
length(r.m.6)

r.m.6.df <- data.frame(gene_id = r.m.6,
                       gene_name = gtf.dat[r.m.6, "gene_name"])

write.table(r.m.6.df, file = paste(results.dir, "rmt.mt.d6.tsv", sep = ""), sep = "\t", row.names = F, quote = F)


# Day 4
day4.all <- unique(Reduce(union, list(rownames(svn.deg[[4]]), rownames(mvn.deg[[4]]), rownames(rvn.deg[[4]]))))

day4.df <- data.frame(gene_id = day4.all,
                      gene_name = gtf.dat[day4.all, "gene_name"])
write.table(day4.df, file = paste(results.dir, "day4_all.tsv", sep = ""), sep = "\t", row.names = F, quote = F)





# SBP unique
svn.unique <- setdiff(svn.genes, unique(union(mvn.genes, rvn.genes)))

svn.unique.df <- data.frame(gene_id = svn.unique,
                            gene_name = gtf.dat[svn.unique, "gene_name"])

write.table(svn.unique.df, file = paste(results.dir, "SBP_unique.tsv", sep = ""), sep = "\t", row.names = F, quote = F)

svn.unique == svn.new
##### Heatmap

## Z-score
heat.dat   <- as.matrix(assay(vsd)[svn.unique, ]);
# heat.dat   <- as.matrix(assay(vsd)[r.m.6, ]);
sample.remove <- grep("d1|d2|d3|d4", colnames(heat.dat), value = TRUE)
heat.dat <- heat.dat[ , !colnames(heat.dat) %in% sample.remove]
heat.dat   <- t(apply(heat.dat,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))}));

# ## TPM
# heat.dat   <- as.matrix(tpm.db[pathway.rt_mt.d1_2.genes.ensembl, ]);
# heat.dat   <- log(heat.dat + 1)
# sample.remove <- grep("d4|d6", colnames(heat.dat), value = TRUE)
# heat.dat <- heat.dat[ , !colnames(heat.dat) %in% sample.remove]


rownames(heat.dat) <- gtf.dat[rownames(heat.dat), "gene_name"];

heat.dat



#order so days are together 
heat.list <- colnames(heat.dat);
heat.list <- as.data.frame(heat.list);
heat.list$day <- ifelse(grepl("naive", heat.list$heat.list),
                        0,
                        str_match(string=heat.list$heat.list, pattern= "\\.d(.+)\\.")[,2]);
heat.list$delivery <- ifelse(grepl("naive", heat.list$heat.list),
                             "naive",
                             str_match(string=heat.list$heat.list, pattern="^(.*)\\..*\\.")[,2]);
heat.list$delivery <- ordered(factor(heat.list$delivery, levels= c('naive', 'mt', 'sbp', 'rtmt')))
heat.list$rep <- as.integer(str_match(string=heat.list$heat.list, pattern=".*r(.)$")[,2]);
heat.list <- heat.list[order(as.numeric(heat.list$day), as.numeric(as.factor(heat.list$delivery)), heat.list$rep),];
heat.anno <- heat.list
heat.list <- as.vector(heat.list$heat.list);


heat.dat <- heat.dat[,heat.list] # Use the same col order as the for the zscore heatmap

ha1 = HeatmapAnnotation(
  df  = data.frame(
    day = as.factor(heat.anno$day)),
  delivery = as.factor(heat.anno$delivery),
  col = list(
    delivery = c(
      "sbp" = "deeppink2",
      "mt"    = "grey50",
      "naive"    = "black",
      "rtmt" = "dodgerblue"),
    day = c(
      "0" = "#FFFFCC",
      "1" = "#D9F0A3",
      "2" = "#ADDD8E",
      "3" = "#78C679",
      "4" = "#31A354",
      "6" = "#006837")
    
    
  )
)

ht <- Heatmap(heat.dat,
              name = "z-score",
              show_row_names          = TRUE,
              show_column_names       = TRUE,
              column_names_side       = "bottom",
              show_column_dend        = FALSE,
              show_row_dend           = FALSE, 
              cluster_rows            = TRUE,
              cluster_columns         = TRUE,
              rect_gp                 = gpar(col = "white", lwd = 1),
              row_dend_width          = unit(30, "mm"),
              column_dend_height      = unit(30, "mm"),
              show_heatmap_legend     = TRUE,
              clustering_distance_rows = 'pearson',
              column_names_max_height = unit(8, "cm"),
              row_names_gp            = gpar(fontsize = 10),
              column_names_gp         = gpar(fontsize = 10),
              # split                   = split.vec, #This would use vector I created above to split and cluster the genes
              top_annotation          = ha1,
              # bottom_annotation       = anno_bottom,
              col                     = colorRamp2(c(min(heat.dat), 0, max(heat.dat)), rev(colorRampPalette(brewer.pal(9, "RdBu"))(3))));
# col                     = viridis::viridis(100));

ht
colnames(heat.dat)

names(res.list)

test <- lapply(res.list[1:10], function(df){
  fil.res <- df[svn.unique, ]
  fil.res
})


test.deg <- lapply(test, function(df){
  df <- df[df$DEG, ]
  df
})

svn.matching <- unique(unlist(sapply(test.deg, row.names)))

svn.new <- setdiff(svn.unique, svn.matching)

lapply(mvn.deg, function(df){
  df <- df[svn.new, ]
  df
})


# Dendogram
library(dendextend)
c <- cor(assay(vsd), method = 'pearson')
c.dist <- as.dist(1-c)
# dist.mat <- dist(t(assay(vsd)), method = "euclidean")
hr.cor <- hclust(c.dist, method = "complete", members=NULL)

dend <- (as.dendrogram(hr.cor))

delivery <- vsd$delivery
cols <- c("black", "deeppink2", "grey50", "dodgerblue") #order in levels of delivery
col_grp <- cols[delivery]
col_grp
col_grp <- col_grp[order.dendrogram(dend)]
par(mfrow = c(1,1))
pdf(paste(results.dir, "dendrogram_colour.pdf", sep = ""), width = 14, height = 7)
dend <- dend  %>% 
  set("labels_colors", col_grp) %>% #change label colors to GROUP
  plot(main = "Dendrogram of sample correlation")
legend("topleft", legend = levels(delivery), fill = cols, cex = 0.75)
dev.off()


Heatmap(as.matrix(hr.cor))
