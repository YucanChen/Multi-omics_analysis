library("ChIPseeker")
library("org.Mm.eg.db")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("GenomicFeatures")
library(VennDiagram)
library(ggplot2)
library(forcats)
library(magrittr)
library(stats)
library(reshape2)
library(stringr)
library(dplyr)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
diff_peaks <- readPeakFile("diffbind_allDB_neg.bed")
diff_peaks <- list(diff_peaks)
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
peakAnnoList <- lapply(diff_peaks, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Mm.eg.db")
peak_genes = unlist(lapply(peakAnnoList, function(i) as.data.frame(i)$geneId))
peak_diff <- read.table("peak_genes_mRNAadded.txt",sep = '\t')[,1]
write.table(peak_genes,"peak_genes.txt",sep = '\t',row.names = FALSE,col.names = FALSE,quote = FALSE)

# Overlap of H3K36me3-regulated genes and DE-miRNA-regulated genes
DEmiRNA_targets <- read.table("ranked_DEmiRNA_negpairs_down_target_entrez.txt")[,1]
tmp <- read.table("ranked_DEmiRNA_negpairs_up_target_entrez.txt")[,1]
DEmiRNA_targets <- c(DEmiRNA_targets,tmp)
peakcovered_DEmiRNAtargets <- read.table("ranked_peakcovered_DEmiRNA_targets_down_entrez.txt")[,1]
tmp <- read.table("ranked_peakcovered_DEmiRNA_targets_up_entrez.txt")[,1]
peakcovered_DEmiRNAtargets <- c(peakcovered_DEmiRNAtargets,tmp)
Venn.H3K36me3_mRNAs <- peak_diff
Venn.miRNAnegtive_pairs <- DEmiRNA_targets
Venn.Peakcover_DEmiRNA_mRNAs <- peakcovered_DEmiRNAtargets
list_combined2<-list("H3K36me3-regulated"=Venn.H3K36me3_mRNAs,"miRNA-regulated"=Venn.miRNAnegtive_pairs,"peakcovered DE-miRNA-regulated"=Venn.Peakcover_DEmiRNA_mRNAs)
venn.diagram(list_combined2,height = 3000, width = 3000,cat.dist= 0.25,cat.pos=c(0,0,0),
             resolution = 300, imagetype = "png", alpha=c(0.45,0.45,0.45),cex=3,
             fill=c("blue","#CC0099","#990000"), cat.fontface=4,fontfamily=3,cat.cex=3,sub.cex=3,
             filename = "VennDiagram_H3K36me3andmiRNA_3CLASS.png")

# enrichment plot
gProfiler_H3K36me3_BP <- read.csv("gProfiler_H3K36me3regulated_upanddown.csv",header = TRUE)
colnames(gProfiler_H3K36me3_BP)
gProfiler_H3K36me3_BP <- gProfiler_H3K36me3_BP[order(gProfiler_H3K36me3_BP$adjusted_p_value),]
graph_pathway <- gProfiler_H3K36me3_BP[1:10,]
graph_pathway %>%
  mutate(term_name = fct_reorder(term_name, -log10(adjusted_p_value))) %>%
  ggplot(aes(x=term_name,y= -log10(adjusted_p_value)))+
  geom_bar(stat='identity',position = "dodge",width=0.9,fill="#990066")+
  ylim(c(0,2))+
  coord_flip() +
  scale_x_discrete(labels=function(x) str_wrap(x, width=40)) +
  xlab("") +
  ylab("-log10(qvalue)") +
  theme(
    axis.text.x=element_text(color="black",size=rel(1.1)),
    axis.text.y=element_text(color="black", size=rel(1.4)),
    axis.title.x = element_text(color="black", size=rel(1.1)),
    legend.text=element_text(color="black",size=rel(0.9)),
    legend.title = element_text(color="black",size=rel(1.0))
  ) +
  geom_hline(yintercept = -1*log10(0.1),colour = "#333399",lty = 2)+
  theme(panel.background = element_rect(color = 'grey60', fill = 'transparent'))
