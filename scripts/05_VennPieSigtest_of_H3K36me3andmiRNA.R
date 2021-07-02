library(openxlsx)
library(ggplot2)
library(ggpubr)
library(VennDiagram)
library(stringr)

# venn plot
## preparing the table
miRNA_exp <- read.xlsx("svanoval_DESeq2_rmnoval_morethan10_FDRadj.xlsx",sheet = "DESeq2_result")
miRNA_exp <- miRNA_exp[,c(1,3,4,8,9,10:16,19)]
colnames(miRNA_exp) <- c("mature_mirbase_id","miRNA_log2FC","miRNA_FC","miRNA_pvalue","miRNA_padj","read_count_Setd2KO1","read_count_Setd2KO2","read_count_Setd2KO3","read_count_Setd2KO4","read_count_Setd2WT1","read_count_Setd2WT2","read_count_Setd2WT3","mirbase_id")
express_table <- read.xlsx("core_table_denovel_mirbase.xlsx")
express_table <- express_table[,c(1:2,14:15,17:21)]
express_table$read_count_Sum <- unlist(lapply(1:length(express_table[,1]),function(x){return(sum(express_table[x,3:9]))}))
express_table <- express_table[which(express_table$read_count_Sum >0),]
peakcovered_premiRNAs <- read.xlsx("peak_covered_premiRNAs.xlsx",sheet = "peakcovered_allpremiRNAs")
peakcovered_premiRNAs <- peakcovered_premiRNAs[,c(1,2,11,13:19)]
colnames(peakcovered_premiRNAs) <- c("mirbase_id","premiRNA_region","width","Conc","Conc_Setd2_knockout","Conc_Setd2_wild.type","H3K36me3_Fold","H3K36me3_pvalue","H3K36me3_FDR","Distance")

#delete no_info chr miRNAs and delete duplicate pre-miRNAs:
all_premiRNA<-read.table("all_miRBase_premiRNA_list.txt",stringsAsFactors = FALSE)
colnames(all_premiRNA)<-c("MI_number","pre-miRNA")
all_premiRNA[412,"pre-miRNA"]<-"mmu-mir-12186_2"
getunique_test <- function(test){
  tmp= by(test[,3:9],
          test$mirbase_id,
          function(x) rownames(x)[which.max(rowMeans(x))])
  probes = as.character(tmp)
  print(dim(test))
  test=test[rownames(test) %in% probes,]
  print(dim(test))
  return(test)
}
### de no-expressed pre-miRNAs
new_express_table<-new_express_table[(new_express_table$mirbase_id %in% all_premiRNA$`pre-miRNA`),]
express_table_rm <- express_table[which(express_table$mirbase_id %in% all_premiRNA$`pre-miRNA`),]
H3K36me3_table <- read.xlsx("AnnoMethylation50000_diffbind_allDB.xlsx")

## draw the graph
Venn.H3K36me3_covered_pre_miRNAs <- H3K36me3_table[,1];Venn.H3K36me3_covered_pre_miRNAs<-as.character(Venn.H3K36me3_covered_pre_miRNAs)
Venn.Nopeak_covered_pre_miRNAs <- all_premiRNA[!(all_premiRNA$`pre-miRNA` %in% H3K36me3_table$Name),2]
Venn.Nopeak_covered_pre_miRNAs <- as.character(Venn.Nopeak_covered_pre_miRNAs)
Venn.top_pre_miRNAs <-new_express_table$mirbase_id
Venn.expressed_pre_miRNAs <-new_express_table$mirbase_id
H3K36me3_coverstate<-unlist(lapply(express_table$mirbase_id, function(x){if(x %in% Venn.H3K36me3_covered_pre_miRNAs){return(TRUE)}else{return(FALSE)}}))
express_table$H3K36me3_covered <- H3K36me3_coverstate
write.table(express_table,"miRNA_expressed_table.xlsx",sep = '\t',row.names = FALSE,col.names = TRUE)

list_combined<-list("H3K36me3-covered"=Venn.H3K36me3_covered_pre_miRNAs,"expressed"=Venn.expressed_pre_miRNAs,"Not covered"=Venn.Nopeak_covered_pre_miRNAs)
venn.diagram(list_combined,height = 3000, width = 3000,cat.dist= 0.25,cat.pos=c(0,0,0),
             resolution = 300, imagetype = "png", alpha=c(0.45,0.45,0.45),cex=2.5,
             fill=c("blue","yellow","purple"), cat.fontface=4,fontfamily=3,cat.cex=3,sub.cex=3,
             #main="Intersection of expressed and H3K36me3 peak-covered pre-miRNAs",
             #main.cex = 2.1, main.fontface = 2, main.fontfamily = 3,
             filename = "VennDiagram_expression_H3K36me3_miRNAoverlap.png")

Venn.DMR_pre_miRNAs <- H3K36me3_table[which(H3K36me3_table$FDR<0.05),][,1];Venn.DMR_pre_miRNAs<-as.character(Venn.DMR_pre_miRNAs)
DEmiRNA_list <- read.xlsx("svanoval_DESeq2_rmnoval_morethan10_FDRadj.xlsx",sheet = "DE_miRNAs")
Venn.DE_pre_miRNAs <- DEmiRNA_list[,"mirbase_id"];Venn.DE_pre_miRNAs<-as.character(Venn.DE_pre_miRNAs)
Venn.DE_pre_miRNAs <- unique(Venn.DE_pre_miRNAs)
Venn.analyzed_pre_miRNAs <- unique(express_table$mirbase_id)
Venn.expressed_pre_miRNAs <- unique(new_express_table$mirbase_id)
list_combined2<-list("H3K36me3-covered"=Venn.H3K36me3_covered_pre_miRNAs,"expressed"=Venn.expressed_pre_miRNAs,"significant_DMR"=Venn.DMR_pre_miRNAs,"significant_DE"=Venn.DE_pre_miRNAs)
venn.diagram(list_combined2,height = 3100, width = 3100,cat.dist= 0.28,cat.pos=c(0,0,0,0),
             resolution = 300, imagetype = "png", alpha=c(0.45,0.45,0.45,0.45),cex=2.6,margin=c(0.01,0.01,0.01,0.01),
             fill=c("red","yellow","blue","green"), cat.fontface=4,fontfamily=3,cat.cex=2.3,sub.cex=2.5,
             #main="Intersection of expressed and H3K36me3 peak-covered pre-miRNAs",
             main.cex = 2.1, main.fontface = 2, main.fontfamily = 3,
             filename = "VennDiagram_significantADD.png")

# significance test
test.set <- c(218,57,573,379)
test.dataframe <- matrix(test.set,nrow = 2,dimnames = list(c("expressed","not expressed"),c("peakcovered","not peakcovered")))
test.fisher_result <- fisher.test(test.dataframe, alternative = "greater")

# pie chart
miRNA_genome <- read.xlsx("pre-miRNA_genome_location.xlsx")
Genome.H3K36me3_covered_pre_miRNAs <- miRNA_genome[miRNA_genome$mirbase_id %in% Venn.H3K36me3_covered_pre_miRNAs,]
Genome.Nopeak_covered_pre_miRNAs <- miRNA_genome[miRNA_genome$mirbase_id %in% Venn.Nopeak_covered_pre_miRNAs,]
Genome.expressed_pre_miRNAs <- miRNA_genome[miRNA_genome$mirbase_id %in% Venn.expressed_pre_miRNAs,]
Genome.Noexpressed_pre_miRNAs <- miRNA_genome[!(miRNA_genome$mirbase_id %in% Venn.expressed_pre_miRNAs),]

## Peak-cover vs. No-cover
### Peak-cover
intergenic=Genome.H3K36me3_covered_pre_miRNAs[which(Genome.H3K36me3_covered_pre_miRNAs$region=="intergenic"),]
intronic=Genome.H3K36me3_covered_pre_miRNAs[which(Genome.H3K36me3_covered_pre_miRNAs$region=="intronic"),]
exonic=Genome.H3K36me3_covered_pre_miRNAs[which(Genome.H3K36me3_covered_pre_miRNAs$region=="exonic"),]
x1 <- c(length(intergenic[,1]), length(intronic[,1]), length(exonic[,1]))
df <- data.frame(
  group = c("intergenic","intronic", "exonic"),
  value = x1,
  piepercent=round(100*x1/sum(x1), 1))
labs <- paste0(df$group, " (", df$piepercent, "%)")
ggpie(df,"piepercent", label = labs,size = 1.2,
      lab.pos = "in",lab.font=c(5,"bold","black"),
      fill = "group",color = "white",
      palette = c("#FF66CC","#3366FF","#9933FF"))

### No-peak cover
intergenic=Genome.Nopeak_covered_pre_miRNAs[which(Genome.Nopeak_covered_pre_miRNAs$region=="intergenic"),]
intronic=Genome.Nopeak_covered_pre_miRNAs[which(Genome.Nopeak_covered_pre_miRNAs$region=="intronic"),]
exonic=Genome.Nopeak_covered_pre_miRNAs[which(Genome.Nopeak_covered_pre_miRNAs$region=="exonic"),]
x2 <- c(length(intergenic[,1]), length(intronic[,1]), length(exonic[,1]))
df <- data.frame(
  group = c("intergenic","intronic", "exonic"),
  value = x2,
  piepercent=round(100*x2/sum(x2), 1))
labs <- paste0(df$group, " (", df$piepercent, "%)")
ggpie(df,"piepercent", label = labs,size = 1.2,
      lab.pos = "in",lab.font=c(5,"bold","black"),
      fill = "group",color = "white",
      palette = c("#FF66CC","#3366FF","#9933FF"))

## Expressed vs. No-exressed
### Expressed
intergenic=Genome.expressed_pre_miRNAs[which(Genome.expressed_pre_miRNAs$region=="intergenic"),]
intronic=Genome.expressed_pre_miRNAs[which(Genome.expressed_pre_miRNAs$region=="intronic"),]
exonic=Genome.expressed_pre_miRNAs[which(Genome.expressed_pre_miRNAs$region=="exonic"),]
x3 <- c(length(intergenic[,1]), length(intronic[,1]), length(exonic[,1]))
df <- data.frame(
  group = c("intergenic","intronic", "exonic"),
  value = x3,
  piepercent=round(100*x3/sum(x3), 1))
labs <- paste0(df$group, " (", df$piepercent, "%)")
# Plot the chart.
#png(file = "paperwriting_plots/part2_PIE_Expressed.png")
ggpie(df,"piepercent", label = labs,size = 1.2,
      lab.pos = "in",lab.font=c(5,"bold","black"),
      fill = "group",color = "white",
      palette = c("#FF66CC","#3366FF","#9933FF"))

### No-expressed
intergenic=Genome.Noexpressed_pre_miRNAs[which(Genome.Noexpressed_pre_miRNAs$region=="intergenic"),]
intronic=Genome.Noexpressed_pre_miRNAs[which(Genome.Noexpressed_pre_miRNAs$region=="intronic"),]
exonic=Genome.Noexpressed_pre_miRNAs[which(Genome.Noexpressed_pre_miRNAs$region=="exonic"),]
x4 <- c(length(intergenic[,1]), length(intronic[,1]), length(exonic[,1]))
df <- data.frame(
  group = c("intergenic","intronic", "exonic"),
  value = x4,
  piepercent=round(100*x4/sum(x4), 1))
labs <- paste0(df$group, " (", df$piepercent, "%)")
# Plot the chart.
#png(file = "paperwriting_plots/part2_PIE_Noexpressed.png")
ggpie(df,"piepercent", label = labs,size = 1.2,
      lab.pos = "in",lab.font=c(5,"bold","black"),
      fill = "group",color = "white",
      palette = c("#FF66CC","#3366FF","#9933FF"))

# group bar chart
colors = c("#FF66CC","#3366FF","#9933FF")
group <- c("Peakcover","No-Peakcover","Expressed","No-expressed")
regions <- c("intergenic","intronic", "exonic")
Values <- matrix(c(x1,x2,x3,x4), nrow = 3, ncol = 4, byrow = FALSE)
df2 <- data.frame(region=rep(regions, each=4),
                  group=rep(group,3),
                  premiRNA_number=c(Values[1,],Values[2,],Values[3,]))
ggbarplot(df2, "group", "premiRNA_number",lab.size = 4.3,size=2,width =0.64,
          fill = "region", color = "region", palette = colors,
          label = TRUE, lab.col = "black", lab.pos = "in")

# volcano plot
## change the FC&log2FC to KOvsWT
p_vocano <- ggplot(miRNA_exp, aes(x = miRNA_log2FC, y = -log10(miRNA_pvalue)))+geom_point(size = 0.8)+labs(x = 'miRNA expression KO/WT_log2FC', y = '-log10(pvalue)')+scale_y_sqrt()
p_vocano# + geom_text(aes(y = -log10(pvalue) + .2, label = Row.names))
par( mfrow = c( 1, 2 ) )
miRNA_exp[which((miRNA_exp$miRNA_pvalue < 0.05) & (miRNA_exp$miRNA_FC > 1.7)),'diff_states'] <- 'up'
miRNA_exp[which((miRNA_exp$miRNA_pvalue < 0.05) & (miRNA_exp$miRNA_FC < -1.7)),'diff_states'] <- 'down'
miRNA_exp[!(miRNA_exp$diff_states %in% c('up', 'down')),'diff_states'] <- 'no'
p_vocano <- ggplot(miRNA_exp, aes(x = miRNA_log2FC, y = -log10(miRNA_pvalue))) +scale_y_sqrt()+
  geom_point(aes(color = diff_states), size = 0.8) +
  scale_colour_manual(limits = c('up', 'down', 'no'), values = c('red', 'blue', 'gray40'), labels = c('Up miRNAs', 'Down miRNAs', 'No diff miRNAs')) +
  labs(x = 'miRNA expression KO/WT_log2FC', y = '-log10(pvalue)')
p_vocano <- p_vocano +
  theme(panel.grid.major = element_line(color = '#CCCCCC', size = 0.1), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  geom_vline(xintercept = c(-log2(1.7),log2(1.7)), color = '#666666', linetype = 2, size = 0.5) + 
  geom_hline(yintercept = -log10(0.05), color = '#666666', linetype = 2, size = 0.5)
p_vocano
