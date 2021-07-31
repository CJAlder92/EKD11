### Functions
count.plot <- function(dds_object, gene_id){
  gene_name <- gtf.dat[gene_id, "gene_name"]
  cnt.cols <- c("mosquito" = "#BCBD22", "blood" = "#E377C2", "naive" = "Black", 'recent_blood' = "#17BECF")
  cnt.dat <- plotCounts(dds_object, gene_id,
                        intgroup = c("day","delivery"), returnData = TRUE)
  # print(cnt.dat)
  cnt.trans <- cnt.dat[1:40,]
  cnt.naive <- cnt.dat[41:44, ]
  naive.mean <- mean(cnt.naive$count)
  cnt.plot <- ggplot(cnt.trans,
                     aes(x = day, y = count, color = delivery, group = delivery)) + 
    ggtitle(gene_name) + 
    geom_point() + stat_summary(fun=mean, geom="line") + 
    geom_hline(yintercept = naive.mean, col = "black", linetype="dashed") +
    scale_x_discrete(limits=factor(1:6)) + 
    scale_colour_manual(values = cnt.cols)
  return(cnt.plot)
}

#### EDIT for sample removed dataset
tpm.plot <- function(gene_id, tpm.data = tpm.db){
  gene_name <- gtf.dat[gene_id, "gene_name"]
  cnt.cols <- c("mosquito" = "#BCBD22", "blood" = "#E377C2", "naive" = "Black", 'recent_blood' = "#17BECF")
  tpm.df <- tpm.data[gene_id, ]
  tpm.df <- melt(as.matrix(tpm.df))
  tpm.df<- cbind(tpm.df ,design[as.character(tpm.df$Var2), c("delivery", "day")])
  tpm.naive <- tpm.df[tpm.df$delivery == "naive", ]
  naive.mean <- mean(tpm.naive$value)
  tpm.inf <- tpm.df[!tpm.df$delivery == "naive", ]
  tpm.plot <- ggplot(tpm.inf,
                     aes(x = day, y = value, group = delivery, fill = delivery)) + 
    ggtitle(gene_name) + 
    geom_bar(position =  position_dodge2(preserve = "single", padding = 0), stat = "summary", fun= "mean", width = 0.75, alpha = 0.7) +
    geom_point(position = position_dodge(width=0.75), aes(fill = delivery)) +
    geom_hline(yintercept = naive.mean, col = "black", linetype="dashed") +
    # scale_x_discrete(limits=factor(1,2,3,4,6)) +
    scale_fill_manual("Delivery", values = cnt.cols, labels = c("SBP", "RMT")) +
    # scale_color_manual(values = cnt.cols) +
    scale_x_discrete("Blood cycles post-infection", labels= c(1,2,3,4,6)) +
    scale_y_continuous("TPM") +
    theme(panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          plot.background = element_blank(),
          legend.position = "none",
          axis.line = element_line(color="black", size = 2),
          text = element_text(size = 30, face = "bold"))
  plot(tpm.plot)
}

heatmap.draw <- function(gene_list, heat.dat = assay(vsd), row_names = FALSE, cluster_samples = FALSE){
  sample.order <- design[colnames(heat.dat), ] %>% 
    arrange(delivery, day) %>% 
    pull(label) %>% as.vector()
  
  heat.genes <- gene_list
  heat.dat <- heat.dat[heat.genes,sample.order ]
  heat.dat   <- t(apply(heat.dat,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))})) %>% data.frame() %>% drop_na() %>% as.matrix()
  rownames(heat.dat) <- gtf.dat[rownames(heat.dat),"gene_name"]
  
  
  heat.dat <- heat.dat[,sample.order]
  max <- max(heat.dat)
  min <- min(heat.dat)
  col.scale <- colorRamp2(seq(min, max, length.out = 100), rev(colorRampPalette(brewer.pal(7, "RdYlBu"))(100)))
  n <- ncol(heat.dat)
  solid.lines <- c(0,4,24,44) # Solid lines for decorate heatmap
  ht <- Heatmap(heat.dat,
                name = "zscore",
                show_row_names          = row_names,
                show_column_names       = TRUE,
                column_names_side       = "bottom",
                show_column_dend        = cluster_samples,
                cluster_columns         = cluster_samples,
                cluster_rows            = TRUE, 
                rect_gp                 = gpar(col = "white", lwd = 1),
                row_dend_width          = unit(30, "mm"),
                column_dend_height      = unit(30, "mm"),
                show_heatmap_legend     = TRUE,
                clustering_distance_rows = 'euclidean',
                # clustering_method_rows = "complete",
                column_names_max_height = unit(8, "cm"),
                row_names_gp            = gpar(fontsize = 8),
                column_names_gp         = gpar(fontsize = 12),
                col                     = col.scale);
  # col                     = colorRamp2(c(min, 0, max), rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(3))));
  
  draw(ht)
  
  decorate_heatmap_body("zscore", c(for(d in seq(0, n, 1)){
    ifelse(d %in% solid.lines, l.type <- 'solid', l.type <- 'dashed')
    grid.lines(x = c(d/n, d/n), y = c(0,1), gp= gpar(col = "black", lty='solid', lwd=2))},
    grid.lines(x = c(0,1), y = c(0,0), gp= gpar(col = "black", lty="solid", lwd=2)),
    grid.lines(x = c(0,1), y = c(1,1), gp= gpar(col = "black", lty="solid", lwd=2))))
}

split.heatmap <- function(gene_list, heat.dat = assay(vsd), z.both = TRUE, row_names = FALSE, cluster_samples = FALSE){
  #order so days are together 
  sample.order <- design[colnames(heat.dat), ] %>% 
    arrange(delivery, day) %>% 
    pull(label) %>% as.vector()
  
  heat.genes <- intersect(gene_list, rownames(heat.dat))
  heat.dat <- heat.dat[heat.genes,sample.order ]
  if(z.both == TRUE){
    heat.dat   <- t(apply(heat.dat,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))})) %>% data.frame() %>% drop_na() %>% as.matrix()}
  rownames(heat.dat) <- gtf.dat[rownames(heat.dat),"gene_name"]
  
  heat.rmt <- heat.dat[, -grep("sbp", colnames(heat.dat))]
  heat.sbp <- heat.dat[, -grep("rtmt", colnames(heat.dat))]
  # heat.sbp <- heat.sbp - rowMeans(heat.sbp)
  if(z.both == FALSE){
    heat.sbp   <- t(apply(heat.sbp,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))})) %>% data.frame() %>% drop_na() %>% as.matrix()
    heat.rmt   <- t(apply(heat.rmt,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))})) %>% data.frame() %>% drop_na() %>% as.matrix()}
  
  
  common.max <- max(c(heat.rmt, heat.sbp))
  common.min <- min(c(heat.rmt, heat.sbp))
  
  col.scale <- colorRamp2(seq(common.min, common.max, length.out = 100), rev(colorRampPalette(brewer.pal(7, "RdYlBu"))(100)))
  # solid.lines <- c(0,4,24,44) # Solid lines for decorate heatmap
  ht.rmt <- Heatmap(heat.rmt,
                    name = "zscore",
                    show_row_names          = FALSE,
                    show_column_names       = TRUE,
                    column_names_side       = "bottom",
                    cluster_columns         = cluster_samples,
                    show_column_dend        = cluster_samples,
                    cluster_rows            = TRUE, 
                    rect_gp                 = gpar(col = "white", lwd = 1),
                    row_dend_width          = unit(30, "mm"),
                    column_dend_height      = unit(30, "mm"),
                    show_heatmap_legend     = TRUE,
                    clustering_distance_rows = 'euclidean',
                    # clustering_method_rows = "complete",
                    column_names_max_height = unit(8, "cm"),
                    row_names_gp            = gpar(fontsize = 8),
                    column_names_gp         = gpar(fontsize = 12),
                    # split = kclus$cluster,
                    col                     = col.scale);
  # col                     = colorRamp2(c(min, 0, max), rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(3))));
  
  
  
  
  
  ht.sbp <- Heatmap(heat.sbp,
                    name = "sbp",
                    show_row_names          = row_names,
                    show_column_names       = TRUE,
                    column_names_side       = "bottom",
                    cluster_columns         = cluster_samples,
                    show_column_dend        = cluster_samples,
                    cluster_rows            = FALSE, 
                    rect_gp                 = gpar(col = "white", lwd = 1),
                    row_dend_width          = unit(30, "mm"),
                    column_dend_height      = unit(30, "mm"),
                    show_heatmap_legend     = FALSE,
                    clustering_distance_rows = 'euclidean',
                    # clustering_method_rows = "complete",
                    column_names_max_height = unit(8, "cm"),
                    row_names_gp            = gpar(fontsize = 8),
                    column_names_gp         = gpar(fontsize = 12),
                    # split = kclus$cluster,
                    col                     = col.scale);
  # col                     = colorRamp2(c(min, 0, max), rev(colorRampPalette(brewer.pal(7, "RdYlBu"))(3))));
  
  
  
  ht_list = ht.rmt + ht.sbp
  ht_list = draw(ht_list, main_heatmap = "zscore") # main_heatmap - dictates row order of other heatmaps based on declared heatmap
  
}


split.heatmap2 <- function(gene_list, heat.dat = assay(vsd), row_names = FALSE, cluster_samples = FALSE){
  #order so days are together 
  sample.order <- design[colnames(heat.dat), ] %>% 
    arrange(delivery, day) %>% 
    pull(label) %>% as.vector()
  
  heat.genes <- intersect(gene_list, rownames(heat.dat))
  heat.dat <- heat.dat[heat.genes,sample.order ]
  heat.dat   <- t(apply(heat.dat,1,function(x){((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))})) %>% data.frame() %>% drop_na() %>% as.matrix()
  rownames(heat.dat) <- gtf.dat[rownames(heat.dat),"gene_name"]
  
  heat.naive <- heat.dat[, grep("naive", colnames(heat.dat))]
  heat.rmt <- heat.dat[, grep("rtmt", colnames(heat.dat))]
  heat.sbp <- heat.dat[, grep("sbp", colnames(heat.dat))]
  # heat.sbp <- heat.sbp - rowMeans(heat.sbp)
  
  
  common.max <- max(c(heat.naive, heat.rmt, heat.sbp))
  common.min <- min(c(heat.naive, heat.rmt, heat.sbp))
  
  col.scale <- colorRamp2(seq(common.min, common.max, length.out = 100), rev(colorRampPalette(brewer.pal(7, "RdYlBu"))(100)))
  
  naive.anno <- HeatmapAnnotation(Treatment = anno_block(gp = gpar(fill = "black"),
                                                         labels = c("Naive"), 
                                                         labels_gp = gpar(col = "white", fontsize = 10)),
                                  Day = anno_text(c(0,0,0,0), rot = 0, location = 0.5))
  ht.naive <- Heatmap(heat.naive,
                      name = "naive",
                      show_row_names          = FALSE,
                      show_column_names       = FALSE,
                      column_names_side       = "bottom",
                      cluster_columns         = cluster_samples,
                      show_column_dend        = cluster_samples,
                      cluster_rows            = TRUE, 
                      show_row_dend           = TRUE, 
                      rect_gp                 = gpar(col = "white", lwd = 1),
                      row_dend_width          = unit(30, "mm"),
                      column_dend_height      = unit(30, "mm"),
                      show_heatmap_legend     = FALSE,
                      clustering_distance_rows = 'euclidean',
                      # clustering_method_rows = "complete",
                      column_names_max_height = unit(8, "cm"),
                      row_names_gp            = gpar(fontsize = 8),
                      column_names_gp         = gpar(fontsize = 12),
                      # split = kclus$cluster,
                      col                     = col.scale,
                      bottom_annotation = naive.anno);
  
  
  rmt.anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "#17BECF"),
                                                 labels = c("RMT"), 
                                                 labels_gp = gpar(col = "white", fontsize = 10)),
                                day = anno_text(str_match(string = colnames(heat.rmt), pattern = "d(.)" )[,2],
                                                location = 0.5, rot = 0))
  
  ht.rmt <- Heatmap(heat.rmt,
                    column_title = "Blood cycles post-infection", 
                    column_title_side = "bottom",
                    name = "zscore",
                    show_row_names          = FALSE,
                    show_column_names       = FALSE,
                    column_names_side       = "bottom",
                    cluster_columns         = cluster_samples,
                    show_column_dend        = cluster_samples,
                    cluster_rows            = TRUE,
                    show_row_dend           = FALSE, 
                    rect_gp                 = gpar(col = "white", lwd = 1),
                    row_dend_width          = unit(30, "mm"),
                    column_dend_height      = unit(30, "mm"),
                    show_heatmap_legend     = TRUE,
                    clustering_distance_rows = 'euclidean',
                    column_names_max_height = unit(8, "cm"),
                    row_names_gp            = gpar(fontsize = 8),
                    column_names_gp         = gpar(fontsize = 12),
                    # split = kclus$cluster,
                    col                     = col.scale,
                    bottom_annotation = rmt.anno);
  # col                     = colorRamp2(c(min, 0, max), rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(3))));
  
  
  sbp.anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "#E377C2"),
                                                 labels = c("SBP"), 
                                                 labels_gp = gpar(col = "white", fontsize = 10)),
                                day = anno_text(str_match(string = colnames(heat.sbp), pattern = "d(.)" )[,2],
                                                location = 0.5, rot = 0))
  
  
  ht.sbp <- Heatmap(heat.sbp,
                    name = "sbp",
                    show_row_names          = row_names,
                    show_column_names       = FALSE,
                    column_names_side       = "bottom",
                    cluster_columns         = cluster_samples,
                    show_column_dend        = cluster_samples,
                    cluster_rows            = FALSE, 
                    rect_gp                 = gpar(col = "white", lwd = 1),
                    row_dend_width          = unit(30, "mm"),
                    column_dend_height      = unit(30, "mm"),
                    show_heatmap_legend     = FALSE,
                    clustering_distance_rows = 'euclidean',
                    # clustering_method_rows = "complete",
                    column_names_max_height = unit(8, "cm"),
                    row_names_gp            = gpar(fontsize = 8),
                    column_names_gp         = gpar(fontsize = 12),
                    # split = kclus$cluster,
                    col                     = col.scale,
                    bottom_annotation = sbp.anno);
  # col                     = colorRamp2(c(min, 0, max), rev(colorRampPalette(brewer.pal(7, "RdYlBu"))(3))));
  
  
  
  ht_list = ht.naive + ht.rmt + ht.sbp
  ht_list = draw(ht_list, main_heatmap = "zscore") # main_heatmap - dictates row order of other heatmaps based on declared heatmap
  
}
