library(openxlsx)
library(biomaRt)
library(clusterProfiler)
library(org.Mm.eg.db)
DEmRNAs <- read.csv("unmerge_DE_gene_results_FC15_pvalue05.csv")
DEmiRNA_predictedtarget <- read.xlsx("multimir_results_data_predicted.xlsx",sheet = "conserved_db")
DEmRNAs_entrez <- bitr(DEmRNAs$gene_name,fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
write.csv(DEmRNAs_entrez,"DEmRNAs_entrez_15_05.csv",row.names = FALSE)
DEmiRNA_exp <- read.xlsx("svanoval_DESeq2_rmnoval_morethan10_FDRadj.xlsx",sheet = "DE_miRNAs")

# DEmiRNA log2FC(KO vs WT), DEmRNA log2FC(KO vs WT)
DEmiRNA_predictedtarget <- subset(DEmiRNA_predictedtarget,DEmiRNA_predictedtarget$target_entrez!="NA")
DEmiRNA_negpairs <- DEmiRNA_predictedtarget[(DEmiRNA_predictedtarget$target_entrez %in% DEmRNAs$entrez_id),]
DEmiRNA_negpairs$miRNA_log2FC <- unlist(lapply(DEmiRNA_negpairs$mature_mirna_id,function(x){DEmiRNA_exp[which(DEmiRNA_exp$Row.names==x),"log2FoldChange"]}))
DEmiRNA_negpairs$target_log2FC <- unlist(lapply(DEmiRNA_negpairs$target_entrez,function(x){DEmRNAs[which(DEmRNAs$entrez_id==x),"log2FoldChange"]}))
DEmiRNA_negpairs <- DEmiRNA_negpairs[which(DEmiRNA_negpairs$miRNA_log2FC*DEmiRNA_negpairs$target_log2FC < 0),]
write.csv(DEmiRNA_negpairs,"part1_miRNAtargets/DEmiRNA_negpairs.csv",row.names = FALSE)

DEmiRNA_negpairs_down <- DEmiRNA_negpairs[which(DEmiRNA_negpairs$target_log2FC < 0),]
DEmiRNA_negpairs_up <- DEmiRNA_negpairs[which(DEmiRNA_negpairs$target_log2FC > 0),]
ranked_DEmiRNA_negpairs_down <- DEmiRNA_negpairs_down[order(abs(DEmiRNA_negpairs_down$target_log2FC),decreasing = TRUE),]
ranked_DEmiRNA_negpairs_up <- DEmiRNA_negpairs_up[order(abs(DEmiRNA_negpairs_up$target_log2FC),decreasing = TRUE),]

write.table(unique(ranked_DEmiRNA_negpairs_down$target_entrez),"ranked_DEmiRNA_negpairs_down_target_entrez.txt",sep = '\n',quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(unique(ranked_DEmiRNA_negpairs_up$target_entrez),"ranked_DEmiRNA_negpairs_up_target_entrez.txt",sep = '\n',quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(ranked_DEmiRNA_negpairs_down,"ranked_DEmiRNA_negpairs_down.xlsx",sep = '\t',quote = FALSE,row.names = FALSE)
write.table(ranked_DEmiRNA_negpairs_up,"ranked_DEmiRNA_negpairs_up.xlsx",sep = '\t',quote = FALSE,row.names = FALSE)
