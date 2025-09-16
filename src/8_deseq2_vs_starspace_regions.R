#! /usr/bin/env Rscript
################################################################################
############################### LOAD LIBRARIES #################################
library("genomation")
library("GenomicDistributions")
library("GenomicDistributionsData")
library("UpSetR")
library("gplots")
library("tidyverse")
library(venn)
library(tornadoplot)

################################ DEFINE FILE PATH ##############################
setwd("/project/shefflab/processed/gomez_atac/results_pipeline/differential")
setwd("//scratch/bx2ur/data")
primary_dir = "primary"
tumoral_dir = "tumoral"

## GITK univ
primary_normal_file = "gitk_univ/deseq/deseq_normal_primary_report_shrink.csv"
as_normal_file = "gitk_univ/deseq/deseq_normal_tumoral_report_shrink.csv"
primary_high_file = "gitk_univ/deseq/deseq_high_primary_report_shrink.csv"

################################## LOAD FILES ##################################
primary_normal <- readGeneric(primary_normal_file, chr = 1, start = 2, end = 3,
                              keep.all.metadata = TRUE, zero.based = FALSE,
                              remove.unusual = FALSE, header = TRUE, 
                              skip = 0, sep = ",")
primary_high <- readGeneric(primary_high_file, chr = 1, start = 2, end = 3,
                            keep.all.metadata = TRUE, zero.based = FALSE,
                            remove.unusual = FALSE, header = TRUE, 
                            skip = 0, sep = ",")

as_normal <- readGeneric(as_normal_file, chr = 1, start = 2, end = 3,
                         keep.all.metadata = TRUE, zero.based = FALSE,
                         remove.unusual = FALSE, header = TRUE, 
                         skip = 0, sep = ",")

starspace_regions = read.csv('/scratch/bx2ur/code/bedembed/outputs/testall/cluster_regions.csv',sep=',')

#################################### MAIN ######################################
#deseq2 regions
absFoldChange = 2
fdr = 0.01

primary_normal <- subset(primary_normal, pvalue < fdr )
primary_normal_up <- subset(primary_normal, log2FoldChange >=absFoldChange )
primary_normal_down <- subset(primary_normal, log2FoldChange <= -absFoldChange )

primary_high <- subset(primary_high, padj < fdr )
primary_high_up <- subset(primary_high, log2FoldChange >= absFoldChange )
primary_high_down <- subset(primary_high, log2FoldChange <= -absFoldChange )

as_normal <- subset(as_normal, padj < fdr )
as_normal_up <- subset(as_normal, log2FoldChange >= absFoldChange )
as_normal_down <- subset(as_normal, log2FoldChange <= -absFoldChange )

# define up down regions

lst_up <- list(paste(as.data.frame(as_normal_up)$seqnames, as.data.frame(as_normal_up)$start, as.data.frame(as_normal_up)$end, sep="_"),
               paste(as.data.frame(primary_normal_up)$seqnames, as.data.frame(primary_normal_up)$start, as.data.frame(primary_normal_up)$end, sep="_"),
               paste(as.data.frame(primary_high_up)$seqnames, as.data.frame(primary_high_up)$start, as.data.frame(primary_high_up)$end, sep="_")
)

lst_down <- list(paste(as.data.frame(as_normal_down)$seqnames, as.data.frame(as_normal_down)$start, as.data.frame(as_normal_down)$end, sep="_"),
                 paste(as.data.frame(primary_normal_down)$seqnames, as.data.frame(primary_normal_down)$start, as.data.frame(primary_normal_down)$end, sep="_"),
                 paste(as.data.frame(primary_high_down)$seqnames, as.data.frame(primary_high_down)$start, as.data.frame(primary_high_down)$end, sep="_")
)


itemList_up <- venn(lst_up, show.plot=FALSE)
lengths(attributes(itemList_up)$intersections)

inc = c(attributes(itemList_up)$intersections$`A:B`, 
        attributes(itemList_up)$intersections$`A:C`,
        attributes(itemList_up)$intersections$`B:C`,
        attributes(itemList_up)$intersections$`A:B:C`)

itemList_down <- venn(lst_down, show.plot=FALSE)
lengths(attributes(itemList_down)$intersections)

dec = c(attributes(itemList_down)$intersections$`A:B`, 
        attributes(itemList_down)$intersections$`A:C`,
        attributes(itemList_down)$intersections$`B:C`,
        attributes(itemList_down)$intersections$`A:B:C`)

# glst_up <- list(primary_high_up, primary_normal_up, as_normal_up)
# up_regions = do.call(c, as(glst_up, "GRangesList"))
# 
# glst_down <- list(primary_high_down, primary_normal_down, as_normal_down)
# down_regions = do.call(c, as(glst_down, "GRangesList"))
# 
# glst_nd <- list(primary_high_nd, primary_normal_nd, as_normal_nd)
# nd_regions = do.call(c, as(glst_nd, "GRangesList"))

inc = inc[!(inc %in% dec)]
dec = dec[!(dec %in% inc)]
#nd = unique(paste(as.data.frame(nd_regions)$seqnames, as.data.frame(nd_regions)$start, as.data.frame(nd_regions)$end, sep="_"))


# starspace regions
ss_up = subset(starspace_regions, title == 'Renin')
ss_down = subset(starspace_regions, title == 'Non-renin')
ss_other = subset(starspace_regions, title == 'other')

r = paste(ss_up$id, ss_up$start, ss_up$end, sep="_")
nr = paste(ss_down$id, ss_down$start, ss_down$end, sep="_")
o = paste(ss_other$id, ss_other$start, ss_other$end, sep="_")

i = r[!(r %in% inc)]
i = i[!(i %in% dec)]

d = nr[!(nr %in% dec)]
d = d[!(d %in% inc)]

nd = o[!(o %in% dec)]
nd = nd[!(nd %in% inc)]

nd = unique(c(nd,i,d))

# upset
listInput_all = list(Increase = inc,
                     Decrease =dec,
                     Non_differential=nd,
                     Renin = paste(ss_up$id, ss_up$start, ss_up$end, sep="_"),
                     Non_renin = paste(ss_down$id, ss_down$start, ss_down$end, sep="_"),
                     Other = paste(ss_other$id, ss_other$start, ss_other$end, sep="_"))

upset(fromList(listInput_all), mainbar.y.label = "overlapped regions", 
      sets = c("Increase", "Decrease", "Non_differential", "Renin", "Non_renin", "Other"),
      keep.order = TRUE, order.by = "freq", text.scale = c(1.75,1.75,1.75,1.75,2,1.75))

sets <- venn(listInput_all, show.plot=FALSE)
lengths(attributes(sets)$intersections)

meta = "~/code/renin_atac/metadata/renin_encode_atac.csv"
# meta data
samples <- read.table(meta, header=TRUE, sep = ",")
samples <- samples[(samples$Group != "Renin_Low"), ] # remove Renin_Low samples

# list of peak files
samples$bigwig= str_replace(samples$bamReads, "_sort_dedup.bam", "_smooth.bw")
files <- list(samples$bigwig)[[1]][1:19]
sample_name <- list(samples$SampleID)[[1]][1:19]


# Declare files
# These bigwigs are derived from PMID: 28212747
sources <- c(
  "/project/shefflab/processed/gomez_atac/results_pipeline/variable/results_pipeline/As4.1_2/aligned_mm10/As4.1_2_smooth.bw",
  "/project/shefflab/processed/gomez_atac/results_pipeline/variable/results_pipeline/Ren1c--_Ren1c-YFP_2/aligned_mm10/Ren1c--_Ren1c-YFP_2_smooth.bw",
  "/project/shefflab/processed/gomez_atac/results_pipeline/variable/results_pipeline/Ren1c-YFP_2/aligned_mm10/Ren1c-YFP_2_smooth.bw",
  "/project/shefflab/processed/gomez_atac/results_pipeline/encode_ATAC/results_pipeline/ENCSR976LWP_1/aligned_mm10/ENCSR976LWP_1_smooth.bw",
  "/project/shefflab/processed/gomez_atac/results_pipeline/encode_ATAC/results_pipeline/ENCSR136XSY_1/aligned_mm10/ENCSR136XSY_1_smooth.bw",
  "/project/shefflab/processed/gomez_atac/results_pipeline/encode_ATAC/results_pipeline/ENCSR064IHX_1/aligned_mm10/ENCSR064IHX_1_smooth.bw"
)
bigwigs <- setNames(sources, c("TN", "PH", "PN", "NR_1", "NR_2", "NR_3"))
# all samples
bigwigs <- setNames(files, sample_name)


inc_gr = data.frame(r)
names(inc_gr)[1] <- "peaks"
inc_gr = inc_gr %>%
  extract(peaks, c("seqname", "start", "end"), regex = "(.*)_([^_]+)_([^_]+)$")
inc_gr <- GenomicRanges::makeGRangesFromDataFrame(inc_gr)

dec_gr = data.frame(nr)
names(dec_gr)[1] <- "peaks"
dec_gr = dec_gr %>%
  extract(peaks, c("seqname", "start", "end"), regex = "(.*)_([^_]+)_([^_]+)$")
dec_gr <- GenomicRanges::makeGRangesFromDataFrame(dec_gr)

# These are a preprocessed set of atac peaks
nd_o = data.frame(attributes(sets)$intersections$`Non_differential:Other`)
names(nd_o)[1] <- "peaks"
nd_o = nd_o %>%
  extract(peaks, c("seqname", "start", "end"), regex = "(.*)_([^_]+)_([^_]+)$")
nd_o = nd_o [order(nd_o$seqname,nd_o$start ),]
write.table(nd_o, "/project/shefflab/processed/gomez_atac/nd_o.bed", row.names=FALSE, quote=FALSE, sep = "\t", col.names = FALSE) 
nd_o <- GenomicRanges::makeGRangesFromDataFrame(nd_o)

nd_r = data.frame(attributes(sets)$intersections$`Non_differential:Renin`)
names(nd_r)[1] <- "peaks"
nd_r = nd_r %>%
  extract(peaks, c("seqname", "start", "end"), regex = "(.*)_([^_]+)_([^_]+)$")
nd_r = nd_r [order(nd_r$seqname,nd_r$start ),]
write.table(nd_r, "/project/shefflab/processed/gomez_atac/nd_r.bed", row.names=FALSE, quote=FALSE, sep = "\t", col.names = FALSE) 
nd_r <- GenomicRanges::makeGRangesFromDataFrame(nd_r)

nd_nr = data.frame(attributes(sets)$intersections$`Non_differential:Non_renin`)
names(nd_nr)[1] <- "peaks"
nd_nr = nd_nr %>%
  extract(peaks, c("seqname", "start", "end"), regex = "(.*)_([^_]+)_([^_]+)$")
nd_nr = nd_nr [order(nd_nr$seqname,nd_nr$start ),]
write.table(nd_nr, "/project/shefflab/processed/gomez_atac/nd_nr.bed", row.names=FALSE, quote=FALSE, sep = "\t", col.names = FALSE) 
nd_nr <- GenomicRanges::makeGRangesFromDataFrame(nd_nr)

inc_r = data.frame(attributes(sets)$intersections$`Increase:Renin`)
names(inc_r)[1] <- "peaks"
inc_r = inc_r %>%
  extract(peaks, c("seqname", "start", "end"), regex = "(.*)_([^_]+)_([^_]+)$")
inc_r = inc_r [order(inc_r$seqname,inc_r$start ),]
write.table(inc_r, "/project/shefflab/processed/gomez_atac/inc_r.bed", row.names=FALSE, quote=FALSE, sep = "\t", col.names = FALSE) 
inc_r <- GenomicRanges::makeGRangesFromDataFrame(inc_r)

dec_o = data.frame(attributes(sets)$intersections$`Decrease:Other`)
names(dec_o)[1] <- "peaks"
dec_o = dec_o %>%
  extract(peaks, c("seqname", "start", "end"), regex = "(.*)_([^_]+)_([^_]+)$")
dec_o = dec_o [order(dec_o$seqname,dec_o$start ),]
write.table(dec_o, "/project/shefflab/processed/gomez_atac/dec_o.bed", row.names=FALSE, quote=FALSE, sep = "\t", col.names = FALSE)
dec_o <- GenomicRanges::makeGRangesFromDataFrame(dec_o)

inc_o = data.frame(attributes(sets)$intersections$`Increase:Other`)
names(inc_o)[1] <- "peaks"
inc_o = inc_o %>%
  extract(peaks, c("seqname", "start", "end"), regex = "(.*)_([^_]+)_([^_]+)$")
inc_o = inc_o [order(inc_o$seqname,inc_o$start ),]
write.table(inc_o, "/project/shefflab/processed/gomez_atac/inc_o.bed", row.names=FALSE, quote=FALSE, sep = "\t", col.names = FALSE)
inc_o <- GenomicRanges::makeGRangesFromDataFrame(inc_o)

dec_r = data.frame(attributes(sets)$intersections$`Decrease:Renin`)
names(dec_r)[1] <- "peaks"
dec_r =dec_r %>%
  extract(peaks, c("seqname", "start", "end"), regex = "(.*)_([^_]+)_([^_]+)$")
dec_r = dec_r [order(dec_r$seqname,dec_r$start ),]
write.table(dec_r, "/project/shefflab/processed/gomez_atac/dec_r.bed", row.names=FALSE, quote=FALSE, sep = "\t", col.names = FALSE)
dec_r <- GenomicRanges::makeGRangesFromDataFrame(dec_r)

dec_nr = data.frame(attributes(sets)$intersections$`Decrease:Non_renin`)
names(dec_nr)[1] <- "peaks"
dec_nr =dec_nr %>%
  extract(peaks, c("seqname", "start", "end"), regex = "(.*)_([^_]+)_([^_]+)$")
dec_nr = dec_nr [order(dec_nr$seqname,dec_nr$start ),]
write.table(dec_nr, "/project/shefflab/processed/gomez_atac/dec_nr.bed", row.names=FALSE, quote=FALSE, sep = "\t", col.names = FALSE)
dec_nr <- GenomicRanges::makeGRangesFromDataFrame(dec_nr)

inc_nr = data.frame(attributes(sets)$intersections$`Increase:Non_renin`)
names(inc_nr)[1] <- "peaks"
inc_nr =inc_nr %>%
  extract(peaks, c("seqname", "start", "end"), regex = "(.*)_([^_]+)_([^_]+)$")
inc_nr = inc_nr [order(inc_nr$seqname,inc_nr$start ),]
write.table(inc_nr, "/project/shefflab/processed/gomez_atac/inc_nr.bed", row.names=FALSE, quote=FALSE, sep = "\t", col.names = FALSE)
inc_nr <- GenomicRanges::makeGRangesFromDataFrame(inc_nr)

###########load files##########
setwd("/scratch/bx2ur/data")
dec_nr = readGeneric("dec_nr.bed", chr = 1, start = 2, end = 3, strand = NULL,
                     meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                     remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
dec_o = readGeneric("dec_o.bed", chr = 1, start = 2, end = 3, strand = NULL,
                    meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                    remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
dec_r = readGeneric("dec_r.bed", chr = 1, start = 2, end = 3, strand = NULL,
                    meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                    remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")

inc_nr = readGeneric("inc_nr.bed", chr = 1, start = 2, end = 3, strand = NULL,
                     meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                     remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
inc_o = readGeneric("inc_o.bed", chr = 1, start = 2, end = 3, strand = NULL,
                    meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                    remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
inc_r = readGeneric("inc_r.bed", chr = 1, start = 2, end = 3, strand = NULL,
                    meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                    remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")

nd_nr = readGeneric("nd_nr.bed", chr = 1, start = 2, end = 3, strand = NULL,
                    meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                    remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
nd_o = readGeneric("nd_o.bed", chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                   remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
nd_r = readGeneric("nd_r.bed", chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                   remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")

# feats <- GRangesList("nd_o" = nd_o, "nd_r" = nd_r, "nd_nr" = nd_nr,
#                      "inc_r" = inc_r,  "inc_o" = inc_o, "dec_o" = dec_o, "dec_nr" = dec_nr)
# 
# feats <- GRangesList("nd_o" = nd_o, "inc_r" = inc_r,  "dec_o" = dec_o, "dec_nr" = dec_nr)
feats <- GRangesList("inc_o" = inc_o, "inc_r" = inc_r,  "dec_o" = dec_o, "dec_nr" = dec_nr)
feats <- GRangesList("nd_o" =sample(nd_o, size=25560, replace = TRUE), 
                     "nd_r" =  sample(nd_r, size=11146, replace = TRUE),  
                     "nd_nr" =  sample(nd_nr, size=5639, replace = TRUE))

#MOTIF#
RARA = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/RARA_filtered.bed', 
                   chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                   remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
MAFK = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/MAFK_filtered.bed', 
                   chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                   remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
SRF = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/SRF_filtered.bed', 
                  chr = 1, start = 2, end = 3, strand = NULL,
                  meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                  remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
PPARA = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/PPARA_filtered.bed', 
                    chr = 1, start = 2, end = 3, strand = NULL,
                    meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                    remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
MAFF = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/MAFF_filtered.bed', 
                   chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                   remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
NR1H4 = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/NR1H4_filtered.bed', 
                    chr = 1, start = 2, end = 3, strand = NULL,
                    meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                    remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
BATF = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/BATF_filtered.bed', 
                   chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                   remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
ESR1 = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/ESR1_filtered.bed', 
                   chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                   remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
MEF2C = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/MEF2C_filtered.bed', 
                    chr = 1, start = 2, end = 3, strand = NULL,
                    meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                    remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
BACH1 = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/BACH1_filtered.bed', 
                    chr = 1, start = 2, end = 3, strand = NULL,
                    meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                    remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
NFE2 = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/NFE2_filtered.bed', 
                   chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                   remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
MEF2D = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/MEF2D_filtered.bed', 
                    chr = 1, start = 2, end = 3, strand = NULL,
                    meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                    remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
MEIS1 = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/MEIS1_filtered.bed', 
                    chr = 1, start = 2, end = 3, strand = NULL,
                    meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                    remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
FOSL1 = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/FOSL1_filtered.bed', 
                    chr = 1, start = 2, end = 3, strand = NULL,
                    meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                    remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
FOS = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/FOS_filtered.bed', 
                  chr = 1, start = 2, end = 3, strand = NULL,
                  meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                  remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
RARB = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/RARB_filtered.bed', 
                   chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                   remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
JUND = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/JUND_filtered.bed', 
                   chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                   remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
FOSL2 = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/FOSL2_filtered.bed', 
                    chr = 1, start = 2, end = 3, strand = NULL,
                    meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                    remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
BACH2 = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/BACH2_filtered.bed', 
                    chr = 1, start = 2, end = 3, strand = NULL,
                    meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                    remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
JUN = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/JUN_filtered.bed', 
                  chr = 1, start = 2, end = 3, strand = NULL,
                  meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                  remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
ATF1 = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/ATF1_filtered.bed', 
                   chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                   remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
ATF3 = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/ATF3_filtered.bed', 
                   chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                   remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
NR4A2 = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/NR4A2_filtered.bed', 
                    chr = 1, start = 2, end = 3, strand = NULL,
                    meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                    remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
FOSB = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/FOSB_filtered.bed', 
                   chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                   remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
JUNB = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/JUNB_filtered.bed', 
                   chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                   remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
CDX1 = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/CDX1_filtered.bed', 
                   chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                   remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")
HMGA1 = readGeneric('/scratch/bx2ur/data/motif_bed/RNAseq_filtered_motifs/HMGA1_filtered.bed', 
                    chr = 1, start = 2, end = 3, strand = NULL,
                    meta.cols = NULL, keep.all.metadata = FALSE, zero.based = FALSE,
                    remove.unusual = FALSE, header = FALSE, skip = 0, sep = "\t")

feats <- GRangesList('RARA' = RARA , 'MAFK' = MAFK , 'SRF' = SRF , 'PPARA' = PPARA , 
                     'MAFF' = MAFF , 'NR1H4' = NR1H4 , 'BATF' = BATF , 'ESR1' = ESR1 , 
                     'MEF2C' = MEF2C , 'BACH1' = BACH1 , 'NFE2' = NFE2 , 'MEF2D' = MEF2D , 
                     'MEIS1' = MEIS1 , 'FOSL1' = FOSL1 , 'FOS' = FOS , 'RARB' = RARB , 
                     'JUND' = JUND , 'FOSL2' = FOSL2 , 'BACH2' = BACH2 , 'JUN' = JUN , 
                     'ATF1' = ATF1 , 'ATF3' = ATF3 , 'NR4A2' = NR4A2 , 'FOSB' = FOSB , 
                     'JUNB' = JUNB , 'CDX1' = CDX1 , 'HMGA1' = HMGA1)

feats <- GRangesList('ESR1' =  ESR1, 'MEF2D' =  MEF2D, 'MEIS1' =  MEIS1, 'ATF3' =  ATF3)
# Recall that the files are named, which the next function uses for labels
# These bigwigs are derived from PMID: 28212747
sources <- c(
  "TN_average.bw",
  "PH_average.bw",
  "PN_average.bw",
  "NR_average.bw"
)

sources <- c(
  "/scratch/bx2ur/data/PN/WT_1_smooth.bw",
  "/scratch/bx2ur/data/PN/WT_2_smooth.bw",
  "/scratch/bx2ur/data/PN/Ren1c-YFP_1_smooth.bw",
  "/scratch/bx2ur/data/PN/Ren1c-YFP_2_smooth.bw",
  "/scratch/bx2ur/data/PH/Ren1c-YFP_Low_Na+Cap_1_smooth.bw",
  "/scratch/bx2ur/data/PH/Ren1c-YFP_Low_Na+Cap_2_smooth.bw",
  "/scratch/bx2ur/data/PH/Ren1c--_Ren1c-YFP_1_smooth.bw",
  "/scratch/bx2ur/data/PH/Ren1c--_Ren1c-YFP_2_smooth.bw",
  "/scratch/bx2ur/data/PH/Recruited_1_smooth.bw",
  "/scratch/bx2ur/data/PH/Recruited_2_smooth.bw",
  "/scratch/bx2ur/data/PH/Recruited_GFP_smooth.bw",
  "/scratch/bx2ur/data/PH/ASKO_1_smooth.bw",
  "/scratch/bx2ur/data/PH/ASKO_2_smooth.bw",
  "/scratch/bx2ur/data/TN/As4.1_2_smooth.bw",
  "/scratch/bx2ur/data/TN/As4.1_H89_Vehicle_1_smooth.bw",
  "/scratch/bx2ur/data/TN/As4.1_H89_Vehicle_2_smooth.bw",
  "/scratch/bx2ur/data/TN/As4.1_Jq1_-_1_smooth.bw",
  "/scratch/bx2ur/data/TN/As4.1_Jq1_-_2_smooth.bw",
  "/scratch/bx2ur/data/TN/As4.1_A485_Vehicle_2_smooth.bw",
  "/scratch/bx2ur/data/NR/P0_mouse_kidney_bulk_ATAC-seq_sample_1_smooth.bw",
  "/scratch/bx2ur/data/NR/P0_mouse_kidney_bulk_ATAC-seq_sample_2_smooth.bw",
  "/scratch/bx2ur/data/NR/3-week_mouse_kidney_bulk_ATAC-seq_sample_1_smooth.bw",
  "/scratch/bx2ur/data/NR/3-week_mouse_kidney_bulk_ATAC-seq_sample_2_smooth.bw",
  "/scratch/bx2ur/data/NR/8-week_mouse_kidney_bulk_ATAC-seq_sample_1_smooth.bw",
  "/scratch/bx2ur/data/NR/8-week_mouse_kidney_bulk_ATAC-seq_sample_2_smooth.bw",
  "/scratch/bx2ur/data/NR/ENCSR023QZX_1_smooth.bw",
  "/scratch/bx2ur/data/NR/ENCSR023QZX_2_smooth.bw",
  "/scratch/bx2ur/data/NR/ENCSR389CLN_1_smooth.bw",
  "/scratch/bx2ur/data/NR/ENCSR389CLN_2_smooth.bw",
  "/scratch/bx2ur/data/NR/ENCSR732OTZ_1_smooth.bw",
  "/scratch/bx2ur/data/NR/ENCSR732OTZ_2_smooth.bw",
  "/scratch/bx2ur/data/NR/ENCSR758IRM_1_smooth.bw",
  "/scratch/bx2ur/data/NR/ENCSR758IRM_2_smooth.bw"
)

bigwigs <- setNames(sources, c("WT_1", "WT_2", "Ren1c-YFP_1", "Ren1c-YFP_2",
                               "Ren1c-YFP_Low_Na+Cap_1", "Ren1c-YFP_Low_Na+Cap_2", "Ren1c--_Ren1c-YFP_1", "Ren1c--_Ren1c-YFP_2", 
                               "Recruited_1", "Recruited_2", "Recruited_GFP", "ASKO_1", "ASKO_2", 
                               "As4.1_2", "As4.1_H89_Vehicle_1", "As4.1_H89_Vehicle_2", "As4.1_Jq1_-_1", "As4.1_Jq1_-_2", "As4.1_A485_Vehicle_2", 
                               "P0_1", "P0_2", "3W_1", "3W_2", "8W_1", "8W_2", "ENCSR023QZX_1", "ENCSR023QZX_2", 
                               "ENCSR389CLN_1", "ENCSR389CLN_2", "ENCSR732OTZ_1", "ENCSR732OTZ_2", "ENCSR758IRM_1", "ENCSR758IRM_2"))
bigwigs <- setNames(sources, c("TN", "PH", "PN", "NR"))
bigwigs <- BigWigFileList(bigwigs)

# Give the function locations and data423464 60 26 13//// 
p = tornado_plot(features = feats, data = bigwigs, width=1000)
p = build_tornado(features = feats, data = bigwigs, width=8000)
df  <- flatten_features(p)

# Function to normalize mean values based on linear regression
normalize_mean <- function(df, feature_set, sample_name, window_size = 500) {
  # Filter dataframe by feature_set and sample_name
  filter_condition <- df$feature_set == feature_set & df$sample_name == sample_name
  subset_df <- df[filter_condition, ]
  
  bg_df <- subset_df[(subset_df$position < -3000 | subset_df$position > 3000), ]
  # Fit linear regression model within the window
  lm_model <- lm(mean ~ position, data = bg_df)
  
  subset_df$mean <- subset_df$mean / predict(lm_model, newdata = data.frame(position = subset_df$position))
  
  # window_size = window_size/25
  # 
  # # Iterate through each position and apply local normalization
  # for (i in 1:nrow(subset_df)) {
  #   # store background here
  #   min_signal = 1
  #   bg = c()
  #   # min slope
  #   s =  0.00000000001
  #   # Define the window around the current position
  #   window_start <- max(1, i - window_size)
  #   window_end <- min(nrow(subset_df), i + window_size)
  #   
  #   # Subset the dataframe within the window
  #   window_df <- subset_df[window_start:window_end, ]
  #   
  #   # Fit linear regression model within the window
  #   lm_model <- lm(mean ~ position, data = window_df)
  #   
  #   # Check if the slope is close to zero
  #   # if (!( window_start >= 80 &&  window_start <= 240)){
  #     if (abs(coef(lm_model)[2]) < s) {
  #       s = abs(coef(lm_model)[2])
  #       # Normalize mean by subtracting the predicted value at the current position
  #       min_signal <- mean(window_df$mean)
  #       
  #       # if (signal < min_signal)
  #       #   min_signal = signal
  #       bg = c(bg, min_signal)
  #     } 
  #   
  # }
  # bg = mean(bg)
  # cat(feature_set, sample_name, ":", bg, "\n")
  # # Replace the original mean values in the subsetted dataframe with normalized values
  # subset_df$mean <- subset_df$mean / bg
  
  # Update the original dataframe with the modified subset
  df[filter_condition, ] <- subset_df
  
  return(df)
}

# Apply the normalization function for each combination of feature_set and sample_name
unique_feature_sets <- unique(df$feature_set)
unique_sample_names <- unique(df$sample_name)

for (sample_name in unique_sample_names) {
  for (feature_set in unique_feature_sets) {
    df <- normalize_mean(df, feature_set, sample_name)
  }
}

# Estimate standard error of the mean
df$se <- sqrt(df$sd^2 / df$n)

df_TN = df[df$sample_name=="TN",] 
df_PN = df[df$sample_name=="PN",]
df_PH = df[df$sample_name=="PH",]

TN_PN = df_TN[,c("mean","position","feature_set", "sample_name", "se")]
TN_PN$mean = log(df_TN$mean / df_PN$mean )
TN_PN$se = log(df_TN$se / df_PN$se)
TN_PN$position = df_TN$position
TN_PN$feature_set = df_TN$feature_set
TN_PN$sample_name = "TN/PN"

PH_PN = df_PH[,c("mean","position","feature_set", "sample_name", "se")]
PH_PN$mean = log(df_PH$mean / df_PN$mean )
PH_PN$se = log(df_PH$se / df_PN$se)
PH_PN$position = df_PH$position
PH_PN$feature_set = df_PH$feature_set
PH_PN$sample_name = "PH/PN"

TN_PH = df_TN[,c("mean","position","feature_set", "sample_name", "se")]
TN_PH$mean = log(df_TN$mean / df_PH$mean )
TN_PH$se = log(df_TN$se / df_PH$se)
TN_PH$position = df_TN$position
TN_PH$feature_set = df_TN$feature_set
TN_PH$sample_name = "TN/PH"

df_plot = rbind(TN_PN,TN_PH, PH_PN)

# by TF class
#bzip
bZip = c("MAFK", "MAFF", "BATF", "BACH1", "NFE2", "FOSL1", "FOS", "JUND", "FOSL2", 
         "BACH2", "JUN", "ATF1", "ATF3", "FOSB", "JUNB")
df_bzip = df_plot[df_plot$feature_set %in% bZip,] 
df_bzip$feature_class = "bZip"
#C4 zinc fingers
C4 = c("RARA", "PPARA", "NR1H4", "ESR1", "RARB", "NR4A2")
df_c4 = df_plot[df_plot$feature_set %in% C4,] 
df_c4$feature_class = "C4 zinc fingers"
#Homeo
Homeo = c("MEIS1", "CDX1")
df_homeo = df_plot[df_plot$feature_set %in% Homeo,] 
df_homeo$feature_class = "Homeo"
#MADS box
MADS = c("MEF2C", "MEF2D", "SRF")
df_mads = df_plot[df_plot$feature_set %in% MADS,] 
df_mads$feature_class = "MADS box"
# TATA-binding
TATA = c("HMGA1")
df_tata = df_plot[df_plot$feature_set %in% TATA,] 
df_tata$feature_class = "TATA-binding"

df_plot = rbind(df_bzip, df_c4, df_homeo, df_mads, df_tata)


# Plotting the result
require(ggplot2)

ggplot(df, aes(position, fill = sample_name)) +
  # geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.2) +
  geom_line(aes(y = mean, colour = sample_name), linewidth = 0.2) +
  # geom_smooth(aes(y = mean, colour = sample_name), linewidth = 0.2, se = FALSE, method = "loess", formula = y ~ x, span = 0.3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~ feature_set) +
  # ylim(-1, 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5)

ggplot(df_plot[df_plot$feature_class=="TATA-binding",], aes(position, fill = feature_set)) +
  # geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.2) +
  # geom_line(aes(y = mean, colour = feature_set), linewidth = 0.2) +
  geom_smooth(aes(y = mean, colour = feature_set), linewidth = 0.4, se = FALSE, method = "loess", formula = y ~ x, span = 0.3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~ sample_name) +
  ylim(-1, 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4)

# Calculating alternative metrics
require(matrixStats)
measure <- list(median = matrixStats::colMedians,
                mad = matrixStats::colMads,
                n = nrow)
df <- flatten_features(tor, measure = measure)
