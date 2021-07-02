library(openxlsx)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(org.Mm.eg.db)
DEmiRNA_exp <- read.xlsx("svanoval_DESeq2_rmnoval_morethan10_FDRadj.xlsx",sheet = "DE_miRNAs")
miRNA_exp <- read.xlsx("svanoval_DESeq2_rmnoval_morethan10_FDRadj.xlsx",sheet = "DESeq2_result")
miRNA_exp <- miRNA_exp[,c(1,3,4,8,9,17,19)]
colnames(miRNA_exp) <- c("mature_mirbase_id","miRNA_log2FC","miRNA_FC","miRNA_pvalue","miRNA_pvalue","miRNA_sum_readcount","mirbase_id")
peakcovered_premiRNAs <- read.xlsx("peak_covered_premiRNAs.xlsx",sheet = "peakcovered_allpremiRNAs")
peakcovered_premiRNAs <- peakcovered_premiRNAs[,c(1,2,11,13:19)]
colnames(peakcovered_premiRNAs) <- c("mirbase_id","premiRNA_region","width","Conc","Conc_Setd2_knockout","Conc_Setd2_wild.type","H3K36me3_Fold","H3K36me3_pvalue","H3K36me3_FDR","Distance")
association_table <- merge(miRNA_exp,peakcovered_premiRNAs,by="mirbase_id")
peakcovered_DEmiRNAs <- DEmiRNA_exp[(DEmiRNA_exp$mirbase_id %in% peakcovered_premiRNAs$mirbase_id),]
DEmiRNA_negpairs <- read.csv("DEmiRNA_negpairs.csv",header = TRUE)
peakcovered_DEmiRNA_targets <- DEmiRNA_negpairs[which(DEmiRNA_negpairs$mature_mirna_id %in% peakcovered_DEmiRNAs$Row.names),]
ranked_peakcovered_DEmiRNA_targets_down <- peakcovered_DEmiRNA_targets[which(peakcovered_DEmiRNA_targets$target_log2FC < 0),]
ranked_peakcovered_DEmiRNA_targets_down <- ranked_peakcovered_DEmiRNA_targets_down[order(abs(ranked_peakcovered_DEmiRNA_targets_down$target_log2FC),decreasing = TRUE),]
ranked_peakcovered_DEmiRNA_targets_up <- peakcovered_DEmiRNA_targets[which(peakcovered_DEmiRNA_targets$target_log2FC > 0),]
ranked_peakcovered_DEmiRNA_targets_up <- ranked_peakcovered_DEmiRNA_targets_up[order(abs(ranked_peakcovered_DEmiRNA_targets_up$target_log2FC),decreasing = TRUE),]

write.table(unique(ranked_peakcovered_DEmiRNA_targets_down$target_entrez),"ranked_peakcovered_DEmiRNA_targets_down_entrez.txt",sep = "\t",col.names = FALSE,row.names = FALSE)
write.table(unique(ranked_peakcovered_DEmiRNA_targets_up$target_entrez),"ranked_peakcovered_DEmiRNA_targets_up_entrez.txt",sep = "\t",col.names = FALSE,row.names = FALSE)

# Association test
association_table[which(association_table$miRNA_pvalue < 0.05 & abs(association_table$miRNA_log2FC) >= 0.7 & (association_table$H3K36me3_FDR>=0.05)),'diff'] <- 'Only DE-miRNAs'
association_table[which(association_table$miRNA_pvalue < 0.05 & abs(association_table$miRNA_log2FC) >= 0.7 & association_table$H3K36me3_FDR<0.05),'diff'] <- 'DMR-DE-miRNAs'
association_table[which((association_table$miRNA_pvalue >= 0.05 | abs(association_table$miRNA_log2FC) < 0.7) & association_table$H3K36me3_FDR<0.05),'diff'] <- 'Only DMR-miRNAs'
association_table[!(association_table$diff %in% c('Only DE-miRNAs', 'DMR-DE-miRNAs', 'Only DMR-miRNAs')),'diff'] <- 'No-diff_miRNAs'
association_table$H3K36me3_log2FC <- -1*log2(-1*association_table$H3K36me3_Fold)
p1 <- ggplot(association_table, aes(x = H3K36me3_log2FC, y = miRNA_log2FC)) +
  geom_point(aes(color = diff), size = 1.1) +
  scale_colour_manual(limits = c('Only DE-miRNAs', 'DMR-DE-miRNAs', 'Only DMR-miRNAs','No-diff_miRNAs'), values = c('#9900CC', 'red', '#66CC00', 'gray40'), labels = c('Only DE-miRNAs', 'DMR-DE-miRNAs', 'Only DMR-miRNAs','No-diff_miRNAs')) +
  labs(x = 'H3K36me3_log2FC', y = 'miRNA_log2FC')
p1 <- p1 +
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'))+
  geom_hline(yintercept = c(0,1,-1), color = 'gray', linetype = 2,size=0.5)+
  geom_vline(xintercept = c(-1), color = 'gray', linetype = 2)+theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"))
p1 <- p1 + geom_smooth(method = "lm", se=FALSE, color="#3333CC", formula = y ~ x, size = 1)
p1 <- p1 + geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x, size = 1)
p1
model.lm <- lm(miRNA_log2FC ~ H3K36me3_log2FC, data = association_table )
summary(model.lm)
cor.test(association_table$miRNA_log2FC,association_table$H3K36me3_log2FC,method="spearman",alternative = c("two.sided"),conf.level = 0.95,exact=FALSE)

## de No_Diff miRNAs：
new_association_table <- association_table[which(association_table$diff!="No-diff_miRNAs"),]
p2 <- ggplot(new_association_table, aes(x = H3K36me3_log2FC, y = miRNA_log2FC)) +
  geom_point(aes(color = diff), size = 1.1) +
  scale_colour_manual(limits = c('Only DE-miRNAs', 'DMR-DE-miRNAs', 'Only DMR-miRNAs'), values = c('#9900CC', 'red', '#66CC00'), labels = c('Only DE-miRNAs', 'DMR-DE-miRNAs', 'Only DMR-miRNAs','No-diff_miRNAs')) +
  labs(x = 'H3K36me3_log2FC', y = 'miRNA_log2FC')
p2 <- p2 +
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'))+
  geom_hline(yintercept = c(0,1,-1), color = 'gray', linetype = 2, size = 0.5)+
  geom_vline(xintercept = c(-1), color = 'gray', linetype = 2, size = 0.5)+theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"))
p2 <- p2 + geom_smooth(method = "lm", se=FALSE, color="#3333CC", formula = y ~ x, size = 1)
p2 <- p2 + geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x, size = 1)
p2
model.lm2 <- lm(miRNA_log2FC ~ H3K36me3_log2FC, data = new_association_table )
summary(model.lm2)
cor.test(new_association_table$miRNA_log2FC,new_association_table$H3K36me3_log2FC,method="spearman",alternative = c("two.sided"),conf.level = 0.95,exact=FALSE)
write.table(association_table,"association_analysis_H3K36me3andmiRNAexpression.xlsx",sep = '\t',row.names = FALSE,col.names = TRUE)

integenic_associate_new <- new_association_table[which(new_association_table$premiRNA_region=="intergenic"),]
intronic_associate_new <- new_association_table[which(new_association_table$premiRNA_region=="intronic"),]
exonic_associate_new <- new_association_table[which(new_association_table$premiRNA_region=="exonic"),]
cor.test(integenic_associate_new$miRNA_log2FC,integenic_associate_new$H3K36me3_log2FC,method="spearman",alternative = c("two.sided"),conf.level = 0.95,exact=FALSE)
cor.test(intronic_associate_new$miRNA_log2FC,intronic_associate_new$H3K36me3_log2FC,method="spearman",alternative = c("two.sided"),conf.level = 0.95,exact=FALSE)
cor.test(exonic_associate_new$miRNA_log2FC,exonic_associate_new$H3K36me3_log2FC,method="spearman",alternative = c("two.sided"),conf.level = 0.95,exact=FALSE)

# plot drawing:
## intergenic
p2_integenic <- ggplot(integenic_associate_new, aes(x = H3K36me3_log2FC, y = miRNA_log2FC)) +
  geom_point(aes(color = diff), size = 1.1) +
  scale_colour_manual(limits = c('Only DE-miRNAs', 'DMR-DE-miRNAs', 'Only DMR-miRNAs'), values = c('#9900CC', 'red', '#66CC00'), labels = c('Only DE-miRNAs', 'DMR-DE-miRNAs', 'Only DMR-miRNAs','No-diff_miRNAs')) +
  labs(x = 'H3K36me3_log2FC', y = 'miRNA_log2FC')
p2_integenic <- p2_integenic +
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'))+
  geom_hline(yintercept = c(0,1,-1), color = 'gray', linetype = 2, size = 0.5)+
  geom_vline(xintercept = c(-1), color = 'gray', linetype = 2, size = 0.5)+theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"))
p2_integenic <- p2_integenic + geom_smooth(method = "lm", se=FALSE, color="#3333CC", formula = y ~ x, size = 1)
p2_integenic <- p2_integenic + geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x, size = 1)
p2_integenic

## intronic
p2_intronic <- ggplot(intronic_associate_new, aes(x = H3K36me3_log2FC, y = miRNA_log2FC)) +
  geom_point(aes(color = diff), size = 1.1) +
  scale_colour_manual(limits = c('Only DE-miRNAs', 'DMR-DE-miRNAs', 'Only DMR-miRNAs'), values = c('#9900CC', 'red', '#66CC00'), labels = c('Only DE-miRNAs', 'DMR-DE-miRNAs', 'Only DMR-miRNAs','No-diff_miRNAs')) +
  labs(x = 'H3K36me3_log2FC', y = 'miRNA_log2FC')
p2_intronic <- p2_intronic +
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'))+
  geom_hline(yintercept = c(0,1,-1), color = 'gray', linetype = 2, size = 0.5)+
  geom_vline(xintercept = c(-1), color = 'gray', linetype = 2, size = 0.5)+theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"))
p2_intronic <- p2_intronic + geom_smooth(method = "lm", se=FALSE, color="#3333CC", formula = y ~ x, size = 1)
p2_intronic <- p2_intronic + geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x, size = 1)
p2_intronic

## exonic
p2_exonic <- ggplot(exonic_associate_new, aes(x = H3K36me3_log2FC, y = miRNA_log2FC)) +
  geom_point(aes(color = diff), size = 1.1) +
  scale_colour_manual(limits = c('Only DE-miRNAs', 'DMR-DE-miRNAs', 'Only DMR-miRNAs'), values = c('#9900CC', 'red', '#66CC00'), labels = c('Only DE-miRNAs', 'DMR-DE-miRNAs', 'Only DMR-miRNAs','No-diff_miRNAs')) +
  labs(x = 'H3K36me3_log2FC', y = 'miRNA_log2FC')
p2_exonic <- p2_exonic +
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'))+
  geom_hline(yintercept = c(0,1,-1), color = 'gray', linetype = 2, size = 0.5)+
  geom_vline(xintercept = c(-1), color = 'gray', linetype = 2, size = 0.5)+theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"))
p2_exonic <- p2_exonic + geom_smooth(method = "lm", se=FALSE, color="#3333CC", formula = y ~ x, size = 1)
p2_exonic <- p2_exonic + geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x, size = 1)
p2_exonic
