library("multiMiR")
library(openxlsx)
library(org.Mm.eg.db)
library(biomaRt)
DEmiRNA_list <- read.xlsx("svanoval_DESeq2_rmnoval_morethan10_FDRadj.xlsx",sheet = "DE_miRNAs")

#Predicted targets:
multimir_results.predicted <- get_multimir(org     = 'mmu',
                                 mirna   = DEmiRNA_list$Row.names,
                                 table   = 'predicted',
                                 predicted.cutoff = 10,
                                 predicted.cutoff.type = "p",
                                 predicted.site = "conserved",
                                 summary = TRUE)
write.csv(as.data.frame(multimir_results.predicted@data),"part1_miRNAtargets/multimir_results_data_predicted_10.csv")
write.csv(as.data.frame(multimir_results.predicted@summary),"part1_miRNAtargets/multimir_results_summary_predicted_10.csv")
