#! /usr/bin/env Rscript
################################################################################
############################### LOAD LIBRARIES #################################
library(DESeq2)
library(preprocessCore)
library(reshape2)
library(tidyverse)
library(forcats)
library(hrbrthemes)
library(viridis)
library(ggpubr)
library(SummarizedExperiment)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)

############################## DEFINE FUNCTIONS ################################
# mutate dataframe for plotting
mutDF = function(res, absFoldChange, fdr){
  res = res[order(res$padj),]
  
  results = as.data.frame(dplyr::mutate(
    as.data.frame(res),
    sig=case_when(res$log2FoldChange <= -absFoldChange & 
                    res$padj < fdr ~ "Down",
                  res$log2FoldChange >= absFoldChange & 
                    res$padj < fdr ~ "Up",
                  TRUE ~ "Not Sig."),
    row.names=rownames(res)))
  return(results)
}

# get significant up & down regions
getUpDownReg = function(res, rld, group, myDesign, absFoldChange, n){
  aCols <- c(1:ifelse(myDesign == "Renin (primary kidney cell) vs Non-renin control",
                      13, 9))     # first n rows are renin
  bCols <- c(ifelse(myDesign == "Renin (primary kidney cell) vs Non-renin control",
                    14, 10):dim(rld)[2])   # ENCODE
  colNums =  c(aCols, bCols)
  # make the lists
  down <- subset(res, log2FoldChange <= -absFoldChange )
  down <- down[ order( down$log2FoldChange ), ]
  down_reg <- rownames(head(down, n=n))
  
  up <- subset(res, log2FoldChange >= absFoldChange )
  up <- up[ order( -up$log2FoldChange ), ]
  up_reg <- rownames(head(up, n=n))
  
  # this gives us the rows we want
  rows <- match(up_reg, row.names(rld))
  mat_up <- assay(rld)[rows,colNums]
  mat_up <- mat_up - rowMeans(mat_up)
  
  rows <- match(down_reg, row.names(rld))
  mat_down <- assay(rld)[rows,colNums]
  mat_down <- mat_down - rowMeans(mat_down)
  
  return(list(mat_down, mat_up))
}

# melt df
meltDF = function(df, renin_samples){
  df_melt <- melt(df)
  df_melt <- df_melt %>% 
    rename(
      consensus_peak_ID = Var1,
      sample_name = Var2,
      Log2FoldChange=value
    )
  df_melt <- df_melt %>%
    # convert consensus_peak_ID to factor and reverse order of levels
    mutate(consensus_peak_ID=factor(consensus_peak_ID,
                                    levels=rev(sort(unique(consensus_peak_ID))))) %>%
    # create a new variable from count
    mutate(FCfactor=cut(df_melt$Log2FoldChange,
                        breaks = seq(-6, 6, by = 1),
                        right = FALSE)) %>%
    # change level order
    mutate(FCfactor=factor(as.character(FCfactor),
                           levels=rev(levels(FCfactor)))) %>%
    # add group
    mutate(df_melt, Group=ifelse(df_melt$sample_name %in% 
                                   renin_samples, "Renin", "Non-renin control"))
  
  return(df_melt)
}

# PLOTS
# MA 
plotDESeqMA = function(res, group, myDesign){
  know_gene = c(
    "Ren1", "Akr1b7"
  )
  # label = res[order(-abs(res$log2FoldChange)),][1:10, ]
  # label_up = res[1:10, ]
  # label_down = res[res$log2FoldChange < -2 ,][order(res$padj),][1:10, ]
  # labels = rbind(label_up, label_down)
  label = res %>% drop_na()
  label = label[label$gene %in%  know_gene,]
  label = label[abs(label$log2FoldChange) >= absFoldChange,]
  label = label[abs(label$padj) <= fdr ,]
  labels = as.data.frame(label)
  print(labels)
  
  # label = rbind(label_up, label_down) ggplot(res_df, aes(x=log2(baseMean), y=log2FoldChange, colour=significant))
  plot <- ggplot2::ggplot(res, ggplot2::aes(x=log2(baseMean), y=log2FoldChange)) +
    ggplot2::geom_point(ggplot2::aes(col = sig), size = 1) +
    ggplot2::scale_color_manual(values = c("#1f77b4", "grey", "#e1812c")) +
    ggplot2::ggtitle(paste("DESeq2 analysis\n", group, myDesign, sep = " ")) +
    ggrepel::geom_text_repel(max.overlaps = Inf,
                             size = 6,
                             data=label,
                             direction = "both",
                             force = 5,
                             ggplot2::aes(label=labels$label)) +
    theme_classic()
  # ggsave(paste(group, myDesign,  "DESeq2 analysis",".svg"),
  #        plot=plot, width = 95, height = 65, units = "mm" )
  plot
}
# plotDESeqMA = function(res, group, myDesign, normalization, fdr){
#   know_gene = c(
#     "Ren1", "Akr1b7"
#   )
#   # label = res[order(-abs(res$log2FoldChange)),][1:10, ]
#   # label_up = res[1:10, ]
#   # label_down = res[res$log2FoldChange < -2 ,][order(res$padj),][1:10, ]
#   # labels = rbind(label_up, label_down)
#   label = res %>% drop_na()
#   label = label[label$gene %in%  know_gene,]
#   label = label[abs(label$log2FoldChange) >= absFoldChange,]
#   label = label[abs(label$padj) <= fdr ,]
#   labels = as.data.frame(label)
#   print(labels)
#   df <- data.frame(res$baseMean, res$log2FoldChange, res$padj <= fdr)
#   # plot <- DESeq2::plotMA(df, ylim = c(-6, 6), colSig = "#6495edff")
#   plot <- ggplot(res_df, aes(x=log2(baseMean), y=log2FoldChange, colour=significant)) +
#     geom_point() +
#     ggrepel::geom_text_repel(max.overlaps = Inf,
#                            size = 4,
#                            data=label,
#                            direction = "both",
#                            force = 5,
#                            ggplot2::aes(label=labels$label)) +
#     theme_classic()
# }

# scree plot
plotScree = function(vsd){
  rv <- rowVars(assay(vsd))
  ## select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
  ## perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(vsd)[select,]))
  ## the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  ##plot the "percentVar"
  scree_plot=data.frame(percentVar)
  scree_plot[,2]<- c(1:dim(vsd)[2])
  colnames(scree_plot)<-c("variance","component_number")
  options(repr.p.width = 30, repr.p.height = 50)
  
  ggplot(scree_plot[1:10,], mapping=aes(x=component_number, y=variance)) +
    geom_bar(stat="identity", fill = "#FF6666") +
    labs(title="Scree plot")+
    xlab("component_number")+
    scale_x_continuous(breaks=c(1:10), labels=c(1:10),limits=c(0,10)) +
    theme_classic()
}

# PCA plot 
plotMyPCA = function(vsd, myDesign){
  data <- plotPCA(vsd, intgroup=contrast, returnData = TRUE)
  percentVar <- round(100 * attr(data, "percentVar"))
  
  p1 <- ggplot(data, aes(PC1, PC2,  color=group)) +   
    geom_point(size=2) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    labs(title=paste(myDesign,  "PCA", sep = " ")) +
    theme_classic() +
    theme(legend.position= "top")
  
  
  p2 <- plotScree(vsd)
  
  grid.arrange(p1, p2, nrow = 1, widths = c(2,1))
}

# volcano plot
plotVolcano = function(res, group, myDesign){
  know_gene = c(
    "Ren1", "Akr1b7"
  )
  # label = res[order(-abs(res$log2FoldChange)),][1:10, ]
  # label_up = res[1:10, ]
  # label_down = res[res$log2FoldChange < -2 ,][order(res$padj),][1:10, ]
  # labels = rbind(label_up, label_down)
  label = res %>% drop_na()
  label = label[label$gene %in%  know_gene,]
  label = label[abs(label$log2FoldChange) >= absFoldChange,]
  label = label[abs(label$padj) <= fdr ,]
  labels = as.data.frame(label)
  print(labels)
  
  # label = rbind(label_up, label_down)
  plot <- ggplot2::ggplot(res, ggplot2::aes(log2FoldChange, -log10(pvalue))) +
    ggplot2::geom_point(ggplot2::aes(col = sig), size = 1) +
    ggplot2::scale_color_manual(values = c("#1f77b4", "grey", "#e1812c")) +
    ggplot2::ggtitle(paste("DESeq2 analysis\n", group, myDesign, sep = " ")) +
    ggrepel::geom_text_repel(max.overlaps = Inf,
                             size = 6,
                             data=label,
                             direction = "both",
                             force = 5,
                             ggplot2::aes(label=labels$label)) +
    theme_classic()
  # ggsave(paste(group, myDesign,  "DESeq2 analysis",".svg"),
  #        plot=plot, width = 95, height = 65, units = "mm" )
  plot
}

# heatmap clustered
plotClusteredHeatmap = function(res, rld, group, myDesign, contrast, absFoldChange, n){
  output = getUpDownReg(res, rld, group, myDesign, absFoldChange, n)
  down = output[[1]]
  up = output[[2]]
  
  df <- as.data.frame(colData(rld)[contrast])
  pheatmap(down, fontsize=8, fontsize_col = 4, annotation_col=df,
           main = paste(group, myDesign, "top 20 down regions"),
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))
                  (length(seq(-6, 6, by = 1))), 
           breaks = seq(-6, 6, by = 1))
  
  pheatmap(up, fontsize=8, fontsize_col = 4, annotation_col=df, 
           main = paste(group, myDesign,"top 20 up regions"),
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))
                  (length(seq(-6, 6, by = 1))), 
           breaks = seq(-6, 6, by = 1))
}

# heatmap 
plotHeatmap = function(df, group, region, myDesign){
  textcol <- "grey40"
  ggplot(df, aes(consensus_peak_ID, sample_name)) + 
    geom_tile(mapping = aes(fill = Log2FoldChange, width=0.9, height=0.9),
              size=0.2) +
    scale_fill_gradient2(high = "darkred", 
                         mid = "white", 
                         low = "darkblue", 
                         midpoint = 0, 
                         n.breaks = 6, 
                         limits=c(-8, 8), 
                         name = NULL) +
    facet_grid(Group ~ ., space="free", scales="free_y", switch="y") +
    labs(title=paste("Top 20", region, "regions", sep = " "))+
    scale_x_discrete(limits = rev, guide = guide_axis(check.overlap = TRUE))+
    scale_y_discrete(guide = guide_axis(check.overlap = TRUE))+
    ylab("sample group")+
    xlab("consensus peak ID")+
    theme(legend.position="right",legend.direction="vertical",
          legend.title=element_text(colour=textcol),
          legend.margin=margin(grid::unit(0,"cm")),
          legend.text=element_text(colour=textcol,size=7,face="bold"),
          legend.key.height=grid::unit(0.8,"cm"),
          legend.key.width=grid::unit(0.2,"cm"),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          strip.text.x = element_text(size = 2),
          axis.text.x=element_text(angle = 90, 
                                   vjust = 0.2, 
                                   hjust=0.95, 
                                   size = 5,
                                  colour=textcol),
          axis.ticks=element_line(size=0.1),
          plot.background=element_blank(),
          panel.border=element_blank(),
          plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
          plot.title=element_text(hjust=0,size=12))
}

plotDESeqHeatmaps = function(res, rld, group, myDesign, renin_samples, contrast, absFoldChange, n){
  output = getUpDownReg(res, rld, group, myDesign, absFoldChange, n)
  down = output[[1]]
  up = output[[2]]
  
  down_melt = meltDF(down, renin_samples)
  up_melt = meltDF(up, renin_samples)
  
  p1 = plotHeatmap(down_melt, group, "down", myDesign)
  p2 = plotHeatmap(up_melt, group, "up", myDesign)
  grid.arrange(p1, p2, nrow = 1, 
               top = text_grob(paste(group, myDesign, sep = " "), vjust=2))
}
df = up_melt
region = "up"
################################# LOAD IMAGE ###################################
setwd("/project/shefflab/processed/gomez_atac/results_pipeline/differential")
primary_dir = "primary"
tumoral_dir = "tumoral"
# load("renin_primary_encode_susztak_atac.RData")
load("renin_encode_atac_0723.RData")

################################ DEFINE FILE PATH ##############################
# myDesign <- "Renin (primary kidney cell) vs Non-renin control"
# myDesign <- "Renin (As4.1 cell line) vs Non-renin control"

# renin_normal_report = paste(ifelse(myDesign == "Renin (primary kidney cell) vs Non-renin control", 
#                                    primary_dir, tumoral_dir), 
#                             "/deseq/deseq_normal_report_no_renin_low.csv", sep="" )
#                             
# renin_high_report = paste(ifelse(myDesign == "Renin (primary kidney cell) vs Non-renin control", 
#                                  primary_dir, tumoral_dir), 
#                           "/deseq/deseq_high_report_no_renin_low.csv", sep="" )
# 
# renin_normal_report_shrink = paste(ifelse(myDesign == "Renin (primary kidney cell) vs Non-renin control", 
#                                           primary_dir, tumoral_dir), 
#                                    "/deseq/deseq_normal_report_no_renin_low_shrink.csv", sep="" )
# 
# renin_high_report_shrink = paste(ifelse(myDesign == "Renin (primary kidney cell) vs Non-renin control", 
#                                         primary_dir, tumoral_dir), 
#                                  "/deseq/deseq_high_report_no_renin_low_shrink.csv", sep="" )

renin_normal_p_report = "gitk_univ/deseq/deseq_normal_primary_report.csv"
renin_normal_t_report = "gitk_univ/deseq/deseq_normal_tumoral_report.csv"
renin_high_p_report = "gitk_univ/deseq/deseq_high_primary_report.csv"

recruited_report = "gitk_univ/deseq/recruited_report.csv"

renin_normal_p_report_shrink = "gitk_univ/deseq/deseq_normal_primary_report_shrink.csv"
renin_normal_t_report_shrink = "gitk_univ/deseq/deseq_normal_tumoral_report_shrink.csv"
renin_high_p_report_shrink = "gitk_univ/deseq/deseq_high_primary_report_shrink.csv"

recruited_report_shrink = "gitk_univ/deseq/recruited_report_shrink.csv"
#################################### MAIN ######################################
# deseq2 analysis
contrast = "Group"

dds <- DESeqDataSetFromMatrix(countData = fc,
                              colData = samples,
                              design = ~Group, 
                              rowRanges = reducedConsensus)

normalizationFactors(dds) <- nfactors
dds <- DESeq(dds)

# result in dataframe format
# res_ren_normal <- results(dds, c(contrast, "Renin_Normal", "Non-Renin_Control"))
# res_ren_high <- results(dds, c(contrast, "Renin_High", "Non-Renin_Control"))

# result in GRanges format
res_ren_normal_p <- results(dds, c(contrast, "Renin_Normal", "Non-Renin_Control"), 
                          format = "GRanges")
res_ren_normal_t <- results(dds, c(contrast, "Renin_Normal", "Non-Renin_Control"), 
                          format = "GRanges")
res_ren_high_p <- results(dds, c(contrast, "Renin_High", "Non-Renin_Control"), 
                        format = "GRanges")

res_recruited <- results(dds, c(contrast, "Recruited", "Native"), 
                                        format = "GRanges")

# save reports
# renin_normal_primary vs ENCODE
write.table( x = data.frame(res_ren_normal_p), 
             file = renin_normal_p_report, 
             sep=",", 
             col.names=TRUE, 
             row.names=FALSE, 
             quote=FALSE )
# renin_normal_tumoral vs ENCODE
write.table( x = data.frame(res_ren_normal_t), 
             file = renin_normal_t_report, 
             sep=",", 
             col.names=TRUE, 
             row.names=FALSE, 
             quote=FALSE )
# renin_high_primary vs ENCODE
write.table( x = data.frame(res_ren_high_p), 
             file = renin_high_p_report, 
             sep=",", 
             col.names=TRUE, 
             row.names=FALSE, 
             quote=FALSE )
#recruited
write.table( x = data.frame(res_recruited), 
             file = recruited_report, 
             sep=",", 
             col.names=TRUE, 
             row.names=FALSE, 
             quote=FALSE )

# shrink log2FoldChange
resLFC_n_p <- lfcShrink(dds, 
                      contrast = c(contrast, "Renin_Normal", "Non-Renin_Control"), 
                      type="ashr", 
                      format = "GRanges")
resLFC_n_t <- lfcShrink(dds, 
                      contrast = c(contrast, "Renin_Normal", "Non-Renin_Control"), 
                      type="ashr", 
                      format = "GRanges")
resLFC_h_p <- lfcShrink(dds, 
                      contrast = c(contrast, "Renin_High", "Non-Renin_Control"), 
                      type="ashr", 
                      format = "GRanges")
#recruited
resLFC_recruited <- lfcShrink(dds, 
                        contrast = c(contrast, "Recruited", "Native"), 
                        type="ashr", 
                        format = "GRanges")

# save files 
# renin_normal_priamry vs ENCODE
write.table( x = data.frame(resLFC_n_p), 
             file = renin_normal_p_report_shrink, 
             sep=",", 
             col.names=TRUE, 
             row.names=FALSE, 
             quote=FALSE )
# renin_normal_tumoral vs ENCODE
write.table( x = data.frame(resLFC_n_t), 
             file = renin_normal_t_report_shrink, 
             sep=",", 
             col.names=TRUE, 
             row.names=FALSE, 
             quote=FALSE )
# renin_high_primary vs ENCODE
write.table( x = data.frame(resLFC_h_p), 
             file = renin_high_p_report_shrink, 
             sep=",", 
             col.names=TRUE, 
             row.names=FALSE, 
             quote=FALSE )

# recruited
write.table( x = data.frame(resLFC_recruited), 
             file = recruited_report_shrink, 
             sep=",", 
             col.names=TRUE, 
             row.names=FALSE, 
             quote=FALSE )

resLFC_n_t = readGeneric(renin_normal_t_report_shrink, chr = 1, start = 2, end = 3, 
                         keep.all.metadata = TRUE, header = TRUE, sep = ",")
resLFC_n_p = readGeneric(renin_normal_p_report_shrink, chr = 1, start = 2, end = 3, 
                         keep.all.metadata = TRUE, header = TRUE, sep = ",")
resLFC_h_p = readGeneric(renin_high_p_report_shrink, chr = 1, start = 2, end = 3, 
                         keep.all.metadata = TRUE, header = TRUE, sep = ",")
# plots
n = 20
absFoldChange = 2
fdr = 0.01

# MA
# plotDESeqMA(res_ren_normal, "Normal", myDesign, "Quantiles Normalized", fdr)
# plotDESeqMA(res_ren_high, "High", myDesign, "Quantiles Normalized", fdr)
# plotDESeqMA(res_ren_t, "High", myDesign, "Quantiles Normalized", fdr)
plotDESeqMA(res_ren_normal, "Normal", myDesign)
plotDESeqMA(res_ren_high, "High", myDesign)
plotDESeqMA(res_ren_t, "High", myDesign)

# PCA
vsd = vst(dds)
plotMyPCA(vsd, myDesign)

# volcano plot
res_ren_normal = mutDF(resLFC_n_p, absFoldChange, fdr)
res_ren_high = mutDF(resLFC_h_p, absFoldChange, fdr)
res_ren_t = mutDF(resLFC_n_t, absFoldChange, fdr)

# tumoral
as_up_gene <- read.csv("/home/bx2ur/Downloads/ph_up_genes.txt", sep = "\t",skip = 1, header=FALSE)
as_up_gene = separate_rows(as_up_gene, V2, sep = ",\\s")
as_up_gene = as_up_gene %>%  separate(V2, c("name", "position"), sep = "\\s")
colnames(as_up_gene) <- c('gene',"name",'position')
as_up_gene = as.data.frame(as_up_gene)
# res_ren_normal$name = row.names(res_ren_normal)
res_ren_normal$name = paste(res_ren_normal$seqnames, res_ren_normal$start, res_ren_normal$end, sep = "_")
res_ren_normal <- res_ren_normal %>% left_join(as_up_gene, by = "name")
res_ren_normal$label <- paste(res_ren_normal$gene, res_ren_normal$position, sep=" ")

# res_ren_high$name = row.names(res_ren_high)
res_ren_high$name = paste(res_ren_high$seqnames, res_ren_high$start, res_ren_high$end, sep = "_")
res_ren_high <- res_ren_high %>% left_join(as_up_gene, by = "name")
res_ren_high$label <- paste(res_ren_high$gene, res_ren_high$position, sep=" ")

#res_ren_t
res_ren_t$name = paste(res_ren_t$seqnames, res_ren_t$start, res_ren_t$end, sep = "_")
res_ren_t <- res_ren_t %>% left_join(as_up_gene, by = "name")
res_ren_t$label <- paste(res_ren_t$gene, res_ren_t$position, sep=" ")

res_ren_normal[res_ren_normal$sig == "up",]
plotVolcano (res_ren_normal, "Normal", myDesign)
plotVolcano (res_ren_high, "High", myDesign)
plotVolcano (res_ren_t, "High", myDesign)

test = subset(res_ren_normal, padj < fdr )
test = subset(res_ren_normal, abs(log2FoldChange) >= absFoldChange )
test = test %>% drop_na()
pheatmap(test, fontsize=8, fontsize_col = 4, annotation_col=df,
         main = paste(group, myDesign, "top 20 down regions"),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))
         (length(seq(-6, 6, by = 1))), 
         breaks = seq(-6, 6, by = 1))
# heatmap
plotClusteredHeatmap(res_ren_normal, vsd, "Normal", myDesign,
                     contrast, absFoldChange, n)
plotClusteredHeatmap(res_ren_high, vsd, "High", myDesign,
                     contrast, absFoldChange, n)

plotDESeqHeatmaps(res_ren_normal, vsd, "Normal", myDesign,
                  renin_samples, contrast, absFoldChange, n)
plotDESeqHeatmaps(res_ren_high, vsd, "High", myDesign,
                  renin_samples, contrast, absFoldChange, n)
group= "Normal"
#################################### SAVE ######################################
# save.image("renin_primary_encode_susztak_atac.RData")
save.image("renin_encode_atac_feature_count.RData")
