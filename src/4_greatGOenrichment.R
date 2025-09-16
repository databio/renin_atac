require("rGREAT")

setwd("/project/shefflab/processed/gomez_atac/results_pipeline/differential")
bed <- read.csv("primary_normal_high_up.bed", sep = "\t")
job = submitGreatJob(bed, species = "mm10")
                        
tb = getEnrichmentTables(job)
names(tb)
go_bio = subset(tb$`GO Biological Process`, tb$`GO Biological Process`$Binom_Fold_Enrichment >= 2 )
plot_this = go_bio [ order( go_bio$Binom_Adjp_BH ), ][1:10,]

x<-plot_this$name
y<- -log10(plot_this$Binom_Adjp_BH)
x <- factor( x, levels = x )
data<-data.frame(x,y)
library(ggplot2)
# Basic barplot
p<-ggplot(data=data, aes(x=x, y=y)) +
  geom_bar(stat="identity", fill="darkseagreen3", width=0.8)
  
# Horizontal bar plot
p + 
  coord_flip()+
  ylab("-log10(FDR)")+
  xlab("")+
  theme_bw() +
  theme(aspect.ratio = 1/1,
        axis.text = element_text(size = 12,colour = "black" ),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


ggplot (data=go_bio, aes(x=name, y=log10(Binom_Adjp_BH))) +        # y takes on negative values
  geom_bar (position = position_dodge(), stat = "identity", fill="salmon", width=0.8) + 
  coord_flip () + 
  scale_x_discrete(name = "", position = "top") +     # x axis (before coord_flip) on opposite side
  scale_y_continuous(name = "-log10(FDR)",
                     breaks = seq(0, -120, by = -30),  # y axis values (before coord_flip) 
                     labels = seq(0, 120, by =  30))+
  theme(aspect.ratio = 1/2,
        axis.text = element_text(size = 12, colour = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 




 
