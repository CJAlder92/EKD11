source('/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/1-R/refactor/packages.R')

nf.dir      <- "/Users/alderc/1-projects/CAMP/1-AS_timecourse/2-EKD11/1-Pipeline/1-Nextflow/"
work.dir    <- "/Users/alderc/1-projects/1-PIR_project/3-EKD11_AS_Early_TP/"
r.dir       <- paste(work.dir, "1-R/",sep='') #change once refactor is completed
tmp.dir     <- paste(work.dir,"tmp/",sep='')
data.dir    <- paste(work.dir,"4-data/",sep='')
results.dir <- paste(work.dir,"2-results_2020/",sep='')
genome.dir  <- "/Users/alderc/1-projects/9-Data/1-Reference_genomes/1-Mus_musculus/"



#### MT vs Naive Unique Data
mt.test.dat <- data.frame(delivery = c(rep("MT", 5), rep("RMT", 5), rep("SBP", 5)),
                          day = rep(c(1,2,3,4,6),3),
                          count = c(10,12,20,40,60,8,11,20,35,58,10,8,9,11,11))

cnt.plot <- ggplot(mt.test.dat, aes(x = day, y = count, color = delivery, fill=delivery, group = delivery)) + geom_smooth(se = FALSE, alpha=0.1, method = "loess", span=0.75);
cnt.plot <- cnt.plot + geom_hline(yintercept = 8, col = "black", linetype="dashed")
cnt.plot <- cnt.plot + ggtitle("Gene X") + ylab("normalised count") + xlab("day")
cnt.plot <- cnt.plot + scale_x_discrete(limits=1:6)
cnt.cols <- c("MT" = "#BCBD22", "SBP" = "#E377C2", "RMT" = "#17BECF")
cnt.plot <- cnt.plot + scale_colour_manual(values = cnt.cols)
png(file=paste(results.dir, "MT_Example_geneprofile.pdf",sep=''),height=500,width=700);
print(cnt.plot);
dev.off();


sbp.test.dat <- data.frame(delivery = c(rep("MT", 5), rep("RMT", 5), rep("SBP", 5)),
                          day = rep(c(1,2,3,4,6),3),
                          count = c(8,8,9,10,12,9,8,10,12,12,10,12,20,40,60))

cnt.plot <- ggplot(sbp.test.dat, aes(x = day, y = count, color = delivery, fill=delivery, group = delivery)) + geom_smooth(se = FALSE, alpha=0.1, method = "loess", span=0.75);
cnt.plot <- cnt.plot + geom_hline(yintercept = 8, col = "black", linetype="dashed")
cnt.plot <- cnt.plot + ggtitle("Gene Y") + ylab("normalised count") + xlab("day")
cnt.plot <- cnt.plot + scale_x_discrete(limits=1:6)
cnt.cols <- c("MT" = "#BCBD22", "SBP" = "#E377C2", "RMT" = "#17BECF")
cnt.plot <- cnt.plot + scale_colour_manual(values = cnt.cols)
png(file=paste(results.dir, "SBP_Example_geneprofile.pdf",sep=''),height=500,width=700);
print(cnt.plot);
dev.off();


#### Examples
arpc4 <- "ENSMUSG00000079426"

res.arpc4 <- lapply(res.list[1:15], function(df){
  df <- df[arpc4, ]
})
arpc4.df <- data.frame(log2FC = cbind(sapply(res.arpc4, function(x){x$log2FC})),
                       FDR = cbind(sapply(res.arpc4, function(x){x$FDR})),
                       DEG = cbind(sapply(res.arpc4, function(x){x$DEG})))
write.table(arpc4.df, file = paste(results.dir, "arpc4_df.csv", sep = ""), sep = ",")


cbx7 <- "ENSMUSG00000053411"

res.cbx7 <- lapply(res.list[1:15], function(df){
  df <- df[cbx7, ]
})
cbx7.df <- data.frame(log2FC = cbind(sapply(res.cbx7, function(x){x$log2FC})),
                       FDR = cbind(sapply(res.cbx7, function(x){x$FDR})),
                       DEG = cbind(sapply(res.cbx7, function(x){x$DEG})))

write.table(cbx7.df, file = paste(results.dir, "cbx7_df.csv", sep = ""), sep = ",")



##### What to do?
res.list <- lapply(res.list,function(res){
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

deg2.list <- lapply(res.list, function(df){
  df <- df[df$DEG, ]
})

## Creating DEG table for each condition (compared against naive)
mvn.deg <- deg2.list[grep('^mt.._vs_naive', names(deg2.list))]
rvn.deg <- deg2.list[grep('^rtmt.._vs_naive', names(deg2.list))]
svn.deg <- deg2.list[grep('^sbp.._vs_naive', names(deg2.list))]


## Create DEG gene lists for each condition
mvn.genes <- unique(unlist(sapply(mvn.deg, row.names)))
rvn.genes <- unique(unlist(sapply(rvn.deg, row.names))) 
svn.genes <- unique(unlist(sapply(svn.deg, row.names)))

mvn.pc <- intersect(mvn.genes, gtf.dat$gene_id[gtf.dat$gene_biotype %in% "protein_coding"])
rvn.pc <- intersect(rvn.genes, gtf.dat$gene_id[gtf.dat$gene_biotype %in% "protein_coding"])
svn.pc <- intersect(svn.genes, gtf.dat$gene_id[gtf.dat$gene_biotype %in% "protein_coding"])


## Venn Diagram 
library("VennDiagram")
venn.cols <- c("#BCBD22", "#E377C2", "#17BECF")
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
                     rotation = 1);

pdf(file = paste(results.dir, "venn/All_PC_no_FC.pdf", sep = ""))
grid.newpage()
grid.draw(venn)
dev.off()

compare_sbp.test.dat <- data.frame(delivery = c(rep("MT", 5), rep("RMT", 5), rep("SBP", 5)),
                          day = rep(c(1,2,3,4,6),3),
                          count = c(10,12,20,40,60,8,11,20,35,58,10,8,9,20,20))

cnt.plot <- ggplot(compare_sbp.test.dat, aes(x = day, y = count, color = delivery, fill=delivery, group = delivery)) + geom_smooth(se = FALSE, alpha=0.1, method = "loess", span=0.75);
cnt.plot <- cnt.plot + geom_hline(yintercept = 8, col = "black", linetype="dashed")
cnt.plot <- cnt.plot + ggtitle("Gene X") + ylab("normalised count") + xlab("day")
cnt.plot <- cnt.plot + scale_x_discrete(limits=1:6)
cnt.cols <- c("MT" = "#BCBD22", "SBP" = "#E377C2", "RMT" = "#17BECF")
cnt.plot <- cnt.plot + scale_colour_manual(values = cnt.cols)
png(file=paste(results.dir, "CompareSBP_Example_geneprofile.pdf",sep=''),height=500,width=700);
print(cnt.plot);
dev.off();





