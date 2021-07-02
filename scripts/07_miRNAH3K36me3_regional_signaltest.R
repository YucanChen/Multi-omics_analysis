library(openxlsx)
library(ggplot2)
library(ggpubr)

deeptools_table <- read.table("Combined_matrix_sampleclass",sep = '\t')
co_name <- c("chr","start","end","no1","no2","strand",paste("WT1","_",c(-1000:-1,1:1000),sep = ""),paste("WT2","_",c(-1000:-1,1:1000),sep = ""),paste("KO","_",c(-1000:-1,1:1000),sep = ""))
colnames(deeptools_table) <- co_name
group <-  c(rep("peak-covered_TSS",201),rep("peak-covered_premiRNA",216),rep("no_peak-covered_TSS",440),rep("no_peak-covered_premiRNA",586))
deeptools_table$group <- c(rep("peak-covered_TSS",201),rep("peak-covered_premiRNA",216),rep("no_peak-covered_TSS",440),rep("no_peak-covered_premiRNA",586))

table_WT1_peakcovered <- deeptools_table[which(deeptools_table$group=="peak-covered_premiRNA"),7:2006]
table_WT1_peakcovered_TSS <- deeptools_table[which(deeptools_table$group=="peak-covered_TSS"),7:2006]
table_WT1_nopeakcovered <- deeptools_table[which(deeptools_table$group=="no_peak-covered_premiRNA"),7:2006]
table_WT1_nopeakcovered_TSS <- deeptools_table[which(deeptools_table$group=="no_peak-covered_TSS"),7:2006]

table_WT1_peakcovered <- t(table_WT1_peakcovered)
table_WT1_peakcovered_TSS <- t(table_WT1_peakcovered_TSS)
table_WT1_nopeakcovered <- t(table_WT1_nopeakcovered)
table_WT1_nopeakcovered_TSS <- t(table_WT1_nopeakcovered_TSS)

WT1_peakcovered_average <- apply(table_WT1_peakcovered,1,function(x) mean(x))
WT1_peakcovered_TSS_average <- apply(table_WT1_peakcovered_TSS,1,function(x) mean(x))
WT1_nopeakcovered_average <- apply(table_WT1_nopeakcovered,1,function(x) mean(x))
WT1_nopeakcovered_TSS_average <- apply(table_WT1_nopeakcovered_TSS,1,function(x) mean(x))

table_WT2_peakcovered <- deeptools_table[which(deeptools_table$group=="peak-covered_premiRNA"),2007:4006]
table_WT2_peakcovered_TSS <- deeptools_table[which(deeptools_table$group=="peak-covered_TSS"),2007:4006]
table_WT2_nopeakcovered <- deeptools_table[which(deeptools_table$group=="no_peak-covered_premiRNA"),2007:4006]
table_WT2_nopeakcovered_TSS <- deeptools_table[which(deeptools_table$group=="no_peak-covered_TSS"),2007:4006]

table_WT2_peakcovered <- t(table_WT2_peakcovered)
table_WT2_peakcovered_TSS <- t(table_WT2_peakcovered_TSS)
table_WT2_nopeakcovered <- t(table_WT2_nopeakcovered)
table_WT2_nopeakcovered_TSS <- t(table_WT2_nopeakcovered_TSS)

WT2_peakcovered_average <- apply(table_WT2_peakcovered,1,function(x) mean(x))
WT2_peakcovered_TSS_average <- apply(table_WT2_peakcovered_TSS,1,function(x) mean(x))
WT2_nopeakcovered_average <- apply(table_WT2_nopeakcovered,1,function(x) mean(x))
WT2_nopeakcovered_TSS_average <- apply(table_WT2_nopeakcovered_TSS,1,function(x) mean(x))

table_KO_peakcovered <- deeptools_table[which(deeptools_table$group=="peak-covered_premiRNA"),4007:6006]
table_KO_peakcovered_TSS <- deeptools_table[which(deeptools_table$group=="peak-covered_TSS"),4007:6006]
table_KO_nopeakcovered <- deeptools_table[which(deeptools_table$group=="no_peak-covered_premiRNA"),4007:6006]
table_KO_nopeakcovered_TSS <- deeptools_table[which(deeptools_table$group=="no_peak-covered_TSS"),4007:6006]

table_KO_peakcovered <- t(table_KO_peakcovered)
table_KO_peakcovered_TSS <- t(table_KO_peakcovered_TSS)
table_KO_nopeakcovered <- t(table_KO_nopeakcovered)
table_KO_nopeakcovered_TSS <- t(table_KO_nopeakcovered_TSS)

KO_peakcovered_average <- apply(table_KO_peakcovered,1,function(x) mean(x))
KO_peakcovered_TSS_average <- apply(table_KO_peakcovered_TSS,1,function(x) mean(x))
KO_nopeakcovered_average <- apply(table_KO_nopeakcovered,1,function(x) mean(x))
KO_nopeakcovered_TSS_average <- apply(table_KO_nopeakcovered_TSS,1,function(x) mean(x))

# significance test
wilcox.test(WT1_peakcovered_average, WT1_peakcovered_TSS_average, paired = TRUE, alternative = "greater",conf.level = 0.95)
wilcox.test(WT2_peakcovered_average, WT2_peakcovered_TSS_average, paired = TRUE, alternative = "greater",conf.level = 0.95)
wilcox.test(KO_peakcovered_average, KO_peakcovered_TSS_average, paired = TRUE, alternative = "greater",conf.level = 0.95)
