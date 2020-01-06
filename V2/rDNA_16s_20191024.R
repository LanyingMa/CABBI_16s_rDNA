
library(phyloseq)
library(ggplot2)
library("data.table")
library(plyr)
library("vegan")
library("DESeq2")
library(ggpubr)
library("VennDiagram")
library(openxlsx)
library(dada2)
library(biomformat)
library(reshape2)
library(cowplot)
library(phylosmith)
library(RColorBrewer)
library(devtools)
library(pairwiseAdonis)

setwd("/Users/lanyingma/OneDrive/Documents/repos/CABBI_16sRDNA")
mis.noBulk <- readRDS( "mis_noBulk.rds")
sample_data(mis.noBulk)$Year <- as.factor(sample_data(mis.noBulk)$Year)
sample_data(mis.noBulk)$Nitrogen <- as.factor(sample_data(mis.noBulk)$Nitrogen)
sample_data(mis.noBulk)$Plot <- as.factor(sample_data(mis.noBulk)$Plot)
##trying different NMDS or PCoA
mis.noBuk.N0 <- prune_samples(sample_data(mis.noBulk)$Nitrogen=="0", mis.noBulk)
mis.noBuk.N200 <- prune_samples(sample_data(mis.noBulk)$Nitrogen=="200", mis.noBulk)
mis.noBuk.N400 <- prune_samples(sample_data(mis.noBulk)$Nitrogen=="400", mis.noBulk)

##NMDS
NMDS.N0 <- plot_ordination(mis.noBuk.N0, ordinate(mis.noBuk.N0,"NMDS"), color = "Year")+stat_ellipse()
ggsave(NMDS.N0, file="figures/NMDS_mis_noBUlk_N0.pdf", units="in", width=6, height=6)
NMDS.N0.pcoa <- plot_ordination(mis.noBuk.N0, ordinate(mis.noBuk.N0,"PCoA"), color = "Year")+stat_ellipse()
ggsave(NMDS.N0.pcoa, file="figures/PCoA_mis_noBUlk_N0.pdf", units="in", width=6, height=6)
NMDS.N0.rda <- plot_ordination(mis.noBuk.N0, ordinate(mis.noBuk.N0,"RDA"), color = "Year")##nope
NMDS.N0.CCA <- plot_ordination(mis.noBuk.N0, ordinate(mis.noBuk.N0,"CCA"), color = "Year")##nope
NMDS.N0.UniFrak <- plot_ordination(mis.noBuk.N0, ordinate(mis.noBuk.N0,"RDA"), color = "Year")##nope


NMDS.N200 <- plot_ordination(mis.noBuk.N200, ordinate(mis.noBuk.N200,"NMDS"), color = "Year")+stat_ellipse()
ggsave(NMDS.N200, file="figures/NMDS_mis_noBUlk_N200.pdf", units="in", width=6, height=6)
NMDS.N200.pcoa <- plot_ordination(mis.noBuk.N200, ordinate(mis.noBuk.N200,"PCoA"), color = "Year")+stat_ellipse()
ggsave(NMDS.N200.pcoa, file="figures/PCoA_mis_noBUlk_N200.pdf", units="in", width=6, height=6)
NMDS.N400 <- plot_ordination(mis.noBuk.N400, ordinate(mis.noBuk.N400,"NMDS"), color = "Year")+stat_ellipse()
ggsave(NMDS.N400, file="figures/NMDS_mis_noBUlk_N400.pdf", units="in", width=6, height=6)
NMDS.N400.pcoa <- plot_ordination(mis.noBuk.N400, ordinate(mis.noBuk.N400,"PCoA"), color = "Year")+stat_ellipse()
ggsave(NMDS.N400.pcoa, file="figures/PCoA_mis_noBUlk_N400.pdf", units="in", width=6, height=6)

sample_data(mis.noBulk)$Plotsub = factor(get_variable(mis.noBulk, "Plot") %in% c("3","4","5","6","7","9","12","13","14"))
mis.noBulk.subplot <- prune_samples(sample_data(mis.noBulk)$Plotsub=="TRUE", mis.noBulk)
pcoa.mis.noBulk.subplot <- plot_ordination(mis.noBulk.subplot, ordinate(mis.noBulk.subplot,"PCoA"), color = "Year", shape="Nitrogen")
ggsave(pcoa.mis.noBulk.subplot, file="figures/pcoa_mis_noBulk_subplot.pdf", units="in", width=6, height=6)
nmds.mis.noBulk.subplot <- plot_ordination(mis.noBulk.subplot, ordinate(mis.noBulk.subplot,"NMDS"), color = "Year", shape="Nitrogen")
ggsave(nmds.mis.noBulk.subplot, file="figures/NMDS_mis_noBulk_subplot.pdf", units="in", width=6, height=6)

## to combine with venn diagram to show who are in each year or who are shared by three years.
mis.noBulk.2015 <- prune_samples(sample_data(mis.noBulk)$Year =="2015", mis.noBulk)
no.zero.mis.noBulk.2015 <- prune_taxa(taxa_sums(mis.noBulk.2015) > 0, mis.noBulk.2015)
list.2015 <- rownames(otu_table(no.zero.mis.noBulk.2015))
temp <- list.2015[!(list.2015 %in% list.2016)]
unique.2015 <- temp[!(temp %in% list.2017)]

mis.noBulk.2016 <- prune_samples(sample_data(mis.noBulk)$Year =="2016", mis.noBulk)
no.zero.mis.noBulk.2016 <- prune_taxa(taxa_sums(mis.noBulk.2016) > 0, mis.noBulk.2016)
list.2016 <- rownames(otu_table(no.zero.mis.noBulk.2016))
temp.2016 <- list.2016[!(list.2016 %in% list.2015)]
unique.2016 <- temp.2016[!(temp.2016 %in% list.2017)]

mis.noBulk.2017 <- prune_samples(sample_data(mis.noBulk)$Year =="2017", mis.noBulk)
no.zero.mis.noBulk.2017 <- prune_taxa(taxa_sums(mis.noBulk.2017) > 0, mis.noBulk.2017)
list.2017 <- rownames(otu_table(no.zero.mis.noBulk.2017))
temp.2017 <- list.2017[!(list.2017 %in% list.2015)]
unique.2017 <- temp.2017[!(temp.2017 %in% list.2016)]

shared.15.16 <- list.2015[list.2015 %in% list.2016]
shared.151617 <- shared.15.16[shared.15.16 %in% list.2017]
shared.151617.phylo <- subset_taxa(mischanthus.noBulk.relative, taxa_names(mischanthus.noBulk.relative) %in% shared.151617)
shared.151617.phylo.glom <-tax_glom(shared.151617.phylo, taxrank="Phylum" )
shared.151617.phylo.glom.melt <- psmelt(shared.151617.phylo.glom )
shared.151617.phylo.glom.melt.mean <- ddply(shared.151617.phylo.glom.melt, c("OTU","Phylum","Year","Nitrogen","Days"),  summarise, mean=mean(Abundance))
colourCount.151617=length(unique((shared.151617.phylo.glom.melt.mean)$Phylum))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
shared.151617.phylo.bar <-  ggplot(data=shared.151617.phylo.glom.melt.mean,  mapping = aes_string(x="as.factor(Year)", y="mean"))+
  geom_bar(aes( fill=Phylum), stat="identity", position="stack")+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(values = getPalette(colourCount.151617))+
  theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size =16),
        axis.title.y = element_text(angle=90, size =16),
        panel.background=element_rect(fill="white", colour="white"),
        axis.title = element_text(size = 10),
        legend.text = element_text(size=10),
        legend.key.size = unit(4, 'mm'),
        legend.title=element_text(size=12),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        strip.background = element_rect(colour = "black", fill = "white",size=2, linetype="solid"),
        strip.text=element_text(size=14)
  )+
  ylim(0.00, 1.00)+
  labs(x="Planting Years of Miscanthus", y="Relative abundance")
ggsave(shared.151617.phylo.bar, file="figures/barofmicrobialCommunityofMischanthus_shared151617.pdf",units="in", width=9, height=9)

unique2015.phy.relative <- subset_taxa(mischanthus.noBulk.relative, taxa_names(mischanthus.noBulk.relative) %in% unique.2015)
unique2015.phy.relative.glom <- tax_glom(unique2015.phy.relative, taxrank="Phylum")
unique2015.phy.melt <- psmelt(unique2015.phy.relative.glom)
unique2015.phy.melt.mean <- ddply(unique2015.phy.melt, c("OTU","Phylum","Year","Nitrogen","Days"),  summarise, mean=mean(Abundance))


unique2016.phy.relative <- subset_taxa(mischanthus.noBulk.relative, taxa_names(mischanthus.noBulk.relative) %in% unique.2016)
unique2016.phy.relative.glom <- tax_glom(unique2016.phy.relative, taxrank="Phylum")
unique2016.phy.melt <- psmelt(unique2016.phy.relative.glom)
unique2016.phy.melt.mean <- ddply(unique2016.phy.melt, c("OTU","Phylum","Year","Nitrogen","Days"), summarise, mean=mean(Abundance))



unique2017.phy.relative <- subset_taxa(mischanthus.noBulk.relative, taxa_names(mischanthus.noBulk.relative) %in% unique.2017)
unique2017.phy.relative.glom <- tax_glom(unique2017.phy.relative, taxrank="Phylum")
unique2017.phy.melt <- psmelt(unique2017.phy.relative.glom)
unique2017.phy.melt.mean <- ddply(unique2017.phy.melt, c("OTU","Phylum","Year","Nitrogen","Days"), summarise, mean=mean(Abundance))


rbind <- rbind(unique2015.phy.melt.mean,unique2016.phy.melt.mean,unique2017.phy.melt.mean)
colourCount=length(unique((rbind)$Phylum))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
unique.15.16.17.bar <- ggplot(data=rbind,  mapping = aes_string(x="as.factor(Year)", y="mean"))+
  geom_bar(aes( fill=Phylum), stat="identity", position="stack")+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(values = getPalette(colourCount))+
  theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size =16),
        axis.title.y = element_text(angle=90, size =16),
        panel.background=element_rect(fill="white", colour="white"),
        axis.title = element_text(size = 10),
        legend.text = element_text(size=10),
        legend.key.size = unit(4, 'mm'),
        legend.title=element_text(size=12),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        strip.background = element_rect(colour = "black", fill = "white",size=2, linetype="solid"),
        strip.text=element_text(size=14)
  )+
  ylim(0.00, 0.15)+
  labs(x="Planting Years of Miscanthus", y="Relative abundance")

shared.unique <- ggarrange(shared.151617.phylo.bar,unique.15.16.17.bar,
          labels=c("Taxa shared by planting year","Taxa unique to each planting year"),
          label.x=0.05,
          label.y=0.95)
ggsave(shared.unique, file="figures/shared_unique_taxaByeachPlantingYear.pdf", units = "in", width=18, height=9)
