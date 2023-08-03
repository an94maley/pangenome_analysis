# install THD package from github if required
# devtools::install_github("rasigadelab/thd")
remove(list = ls())
library(tidyverse)
library(ggplot2)
library(ggpattern)
library(ggbreak)
library(pheatmap)
library(openxlsx)
library(thd)
library(ape)
library(phytools)
library(cowplot)
setwd("C:/Users/22674/Desktop/pangenome_analysis")


####Figure 2####
data_pie1 <- read.xlsx("core_pan_genes.xlsx",
                       sheet = "pie_plot",
                       rows = c(1,2,3,4),
                       cols = c(1,2))
data_pie1$percent <- data_pie1$freq/sum(data_pie1$freq)*100
data_pie2 <- read.xlsx("core_pan_genes.xlsx",
                       sheet = "pie_plot",
                       rows = c(1,2,3),
                       cols = c(4,5))
data_pie2$percent <- data_pie2$freq/sum(data_pie2$freq)*100
data_pie2$class <- factor(data_pie2$class, levels = c("Others", "Core genes"))
pie_plot1 <- ggplot(data_pie1, aes(x="", y=percent, fill=class)) +
  geom_bar(width=0.8, stat="identity") +
  scale_fill_manual(values = c("#00EEEE","#32CD32","#FFFF00")) +
  coord_polar(theta="y") +
  geom_text(aes(label=paste0(round(percent,2),"%")),
            position=position_stack(vjust=0.5),size=6) +
  theme_classic() +
  theme(panel.background=element_rect(fill="transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=13),
        legend.spacing.y = unit(2, "cm")
  )
ggsave(pie_plot1, filename = "pie_plot1.pdf", width = 7, height = 7)
pie_plot2 <- ggplot(data_pie2, aes(x="", y=percent)) +
  geom_bar_pattern(aes(pattern_fill=class), 
                   stat="identity",
                   width=0.8,
                   pattern_angle=30,
                   colour="black",
                   fill="#FFFF00",
                   pattern="stripe",
                   pattern_color=c("#FF0000", "black"),
                   pattern_spacing = 0.03,
                   pattern_density = 0.5,
                   pattern_linetype = 1) +
  coord_polar("y", start = 0) +
  scale_pattern_fill_manual(values=c("black", "#FF0000")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        legend.title=element_blank(),
        legend.text = element_text(size=13)
        )
ggsave(pie_plot2, filename = "pie_plot2.pdf", width = 7, height = 7)
# pie_plot1 and pie_plot2 were combined into Figure 2 using Adobe Illustrator software


####Figure 3####
data_genecluster <- read.xlsx("core_pan_genes.xlsx", sheet = "bar_plot")
data_genecluster$genome <- factor(data_genecluster$genome)
data_genecluster$class <- factor(data_genecluster$class, levels = c("Cloud genes", "Shell genes", "Soft-core genes", "Core genes"))
plot_genecluster <- ggplot(data_genecluster, aes(x=genome, y=freq, fill=class)) +
  geom_bar(width=0.8, stat="identity") +
  scale_fill_manual(values = c("#00EEEE","#32CD32","#FFFF00","#FF0000")) +
  scale_y_continuous(breaks=c(seq(0,600,100), seq(1000,1600,200), seq(17500,17600,50)),
                     limits = c(0,17600),
                     expand = c(0,0)) +
  scale_y_break(c(1600, 17450), scales = 0.3, space = 0.25) +
  scale_y_break(c(650, 1000), scales = 0.3, space = 0.25) +
  labs(x="No. of genomes in clusters", y="No. of gene clusters") +
  theme_classic() +
  theme(panel.background=element_rect(fill="transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size=14, color="black"),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=13, color="black"),
        axis.title.x = element_text(size=16, color="black"),
        axis.title.y = element_text(size=17, color="black"),
        legend.position = "top",
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.spacing.x = unit(2, "mm"),
        legend.text = element_text(size=14)
        )
ggsave(plot_genecluster, filename = "Figure 3.pdf", width = 15, height = 9)


####Figure 6####
heat_data <- read.xlsx("average_nucleotide_identity.xlsx", colNames = T, rowNames = T)
data <- read.xlsx("basic_characteristic.xlsx",
                  colNames = T, rowNames = T)
annotation_row <- data %>% select(Lineage)
annotation_col <- data %>% select(DR)
# create heatmap
pdf(file="Figure 6.pdf", height=8, width=11)
pheatmap(heat_data,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         annotation_row = annotation_row,
         annotation_col = annotation_col,
         fontsize_row = 5.6, fontsize_col = 5.4, angle_col = 45,
         legend_breaks = c(0.998,0.9985,0.999,0.9995,1),
         border_color = "NA",
         legend_labels = c("99.80%","99.85%","99.90%","99.95%","100%")
         )
dev.off()


####Figure 7####
kegg_data <- read.xlsx("KEGG_COG_histogram.xlsx",
                        colNames = T, sheet = "kegg_detail")
kegg_data$type <- factor(kegg_data$type, levels=c("Core genome","Accessory genome","Unique genes"))
order <- c("Amino acid metabolism", "Biosynthesis of other secondary metabolites", "Carbohydrate metabolism", "Energy metabolism", "Glycan biosynthesis and metabolism", "Lipid metabolism", "Metabolism of cofactors and vitamins", "Metabolism of other amino acids", "Metabolism of terpenoids and polyketides", "Nucleotide metabolism", "Overview", "Xenobiotics biodegradation and metabolism", "Cancers", "Cardiovascular diseases", "Drug resistance", "Endocrine and metabolic diseases", "Immune diseases", "Infectious diseases", "Neurodegenerative diseases", "Substance dependence", "Folding sorting and degradation", "Replication and repair", "Transcription", "Translation", "Circulatory system", "Development", "Digestive system", "Endocrine system", "Environmental adaptation", "Excretory system", "Immune system", "Nervous system", "Sensory system", "Membrane transport", "Signal transduction", "Signaling molecules and interaction", "Cell growth and death", "Cell motility", "Cellular commiunity", "Transport and catabolism")
kegg_data$secondary <- factor(kegg_data$secondary, levels = order)
kegg_data_major <- read.xlsx("KEGG_COG_histogram.xlsx",
                             colNames = T, sheet = "kegg_sum")
kegg_data_major$type <- factor(kegg_data_major$type, levels = c("Core genome","Accessory genome","Unique genes"))
kegg_plot1 <- ggplot(kegg_data, aes(x=secondary,y=value, fill=type)) +
                geom_bar(stat = "identity", width = 0.6) +
                scale_y_continuous(breaks = seq(0,40,5), limits = c(0, 58), expand = c(0,0)) +
                scale_fill_manual(values = c("#FF0000", "#FFD700", "#00EEEE")) +
                labs(x="Secondary KEGG pathway", y="Percentage (%)") +
                theme_bw() +
                theme(panel.background=element_rect(fill="transparent"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.text.y = element_text(size=13, color="black"),
                      axis.text.x = element_text(angle=30, hjust=1, size=12, color="black"),
                      axis.title = element_text(size=15, color="black", face = "bold"),
                      legend.position = "left", 
                      legend.background = element_blank(), 
                      legend.title = element_blank(),
                      legend.spacing.y = unit(5, "mm"),
                      legend.text = element_text(size=13)
                      )
kegg_plot2 <- ggplot(kegg_data_major, aes(x=major, y=value, fill=type)) +
                geom_bar(stat = "identity", width = 0.75, position = "dodge") +
                scale_y_continuous(breaks = seq(0,90,10), limits = c(0,85), expand = c(0,0)) +
                scale_fill_manual(values = c("#FF0000", "#FFD700", "#00EEEE")) +
                labs(x="Major KEGG pathway", y="Percentage (%)") +
                coord_flip() +
                theme_bw() +
                theme(panel.background=element_rect(fill="transparent"),
                      panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      axis.text.y = element_text(size=13, color="black"),
                      axis.text.x = element_text(hjust=0.5, size=13, color="black"),
                      axis.title = element_text(size=13, color="black", face = "bold"),
                      legend.position = "none")
p2_grob <- ggplotGrob(kegg_plot2)
kegg_plot <- kegg_plot1 + annotation_custom(grob = p2_grob,xmin = 5,xmax = 40,ymin = 14,ymax = 56)
ggsave(kegg_plot, filename="Figure 7.pdf", height = 8, width = 15)


####Figure S3####
data_each <- read.xlsx("core_pan_genes.xlsx",
                       sheet = "individual")
data_each$prop <- data_each$prop*100
order <- sort(unique(data_each$isolates), decreasing = T)
data_each$isolates <- factor(data_each$isolates, levels = order)
data_each$type <- factor(data_each$type, levels = c("Soft-core genes", "Shell genes", "Cloud genes"))
data_each1 <- data_each[c(seq(1,108,1)),]
data_each2 <- data_each[c(seq(109,216,1)),]

plot1 <- ggplot(data_each1, aes(x=isolates, y=prop, fill=type)) +
  geom_bar(width=0.65, stat="identity") +
  scale_fill_manual(values = c("#FFFF00","#32CD32","#00EEEE")) +
  scale_y_continuous(breaks = seq(0,100,10), expand = c(0,0)) +
  labs(y="Percentage (%)") +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.key.width = unit(0.5,units = "cm"),
        legend.key.height = unit(0.5,units = "cm"),
        legend.text = element_text(size=13),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14, color="black", face = "bold"),
        axis.text = element_text(size=13, color="black")
  )
plot2 <- ggplot(data_each2, aes(x=isolates, y=prop, fill=type)) +
  geom_bar(width=0.65, stat="identity") +
  scale_fill_manual(values = c("#FFFF00","#32CD32","#00EEEE")) +
  scale_y_continuous(breaks = seq(0,100,10), expand = c(0,0)) +
  labs(y="Percentage (%)") +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.key.width = unit(0.5,units = "cm"),
        legend.key.height = unit(0.5,units = "cm"),
        legend.text = element_text(size=13),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14, color="black", face = "bold"),
        axis.text = element_text(size=13, color="black")
  )
plot_1_2 <- plot_grid(plot1, plot2, ncol=2, rel_widths = c(1,1.07))
ggsave(plot_1_2, filename = "Figure S3.pdf", width = 15, height = 8)


####Figure S4####
tree1 <- read.tree("best_PGM_IQT_UFBoot_run5_GTR2+FO+R3.nwk")
tree2 <- read.tree("snp_phylogeny.nwk")
association <- cbind(sort(tree1$tip.label), sort(tree2$tip.label))
colnames(association) <- c("tree1", "tree2")
graph <- cophylo(tree1, tree2, assoc=association)

pdf(file="Figure S4.pdf", height=8, width=6)
plot(graph, link.type="straight", link.lwd=2, link.lty="solid", 
     link.col=make.transparent("red", 0.5),
     length.line=0.7, font=3, fsize=0.6, space=20
)
dev.off()


####Figure S5####
cog_data <- read.xlsx("KEGG_COG_histogram.xlsx",
                       colNames = T, sheet = "cog_detail")
cog_data$type <- factor(cog_data$type, levels=c("Core genome","Accessory genome","Unique genes"))
cog_plot <- ggplot(cog_data, aes(x=major,y=value, fill=type)) +
  geom_hline(yintercept=seq(5,15,5),color="black",size=0.5,linetype="dashed") +
  geom_bar(stat = "identity", width = 0.75, position = "dodge") +
  scale_y_continuous(breaks = seq(0,15,5), limits = c(0, 17), expand = c(0,0)) +
  scale_fill_manual(values = c("#FF0000", "#FFD700", "#00EEEE")) +
  labs(x="COG categories", y="Percentage (%)") +
  theme_bw() +
  theme(panel.background=element_rect(fill="transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size=15, color="black"),
        axis.text.x = element_text(angle=30, hjust=1, size=15, color="black"),
        axis.title = element_text(size=16, color="black", face = "bold"),
        legend.position = "left", 
        legend.background = element_blank(), 
        legend.title = element_blank(),
        legend.spacing.y = unit(5, "mm"),
        legend.text = element_text(size=14)
  )
ggsave(cog_plot, filename="Figure S5.pdf", height = 10, width = 15)


####THD calculation####
dna <- read.dna("snp_seq.fa", "fasta")
dna_len <- 9566
# SNP distance matrix
H <- dist.dna(dna, "raw", pairwise.deletion = T, as.matrix = T) * dna_len
# parameter definition
time <- 200
m <- 4000000    # no. of markers
mu <- 1e-7      # substitution rate per marker
# calculate THD index
thd <- thd(H, time, m, mu)*10000
lab <- labels(dna)
dat_thd <- data.frame(lab, thd)
dat_characteristic <- read.xlsx("basic_characteristic.xlsx")
data <- merge(dat_characteristic, dat_thd, by.x="SampleID", by.y="lab")
# the median of THD index
median(data$thd)
# between-group comparison of THD index
kruskal.test(thd ~ Region, data)
kruskal.test(thd ~ Lineage, data)
# data were saved for pan-genome association analysis
write.csv(data, file="data.csv")
