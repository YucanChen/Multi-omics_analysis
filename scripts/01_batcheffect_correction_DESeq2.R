library(openxlsx)
library(DESeq2)
library(sva)
newmatrix_miRNAexp <- read.table("core_table.xls",sep='\t',header = TRUE)
newmatrix_miRNAexp <- newmatrix_miRNAexp[,c(1,14:21)]
rownames(newmatrix_miRNAexp) <- newmatrix_miRNAexp[,1]
newmatrix_miRNAexp <- newmatrix_miRNAexp[,-c(1,4)]
newmatrix_miRNAexp <- newmatrix_miRNAexp[rowSums(newmatrix_miRNAexp)>0,]
#attach(newmatrix_miRNAexp)
#newmatrix_miRNAexp<-newmatrix_miRNAexp[which(read_count_Setd2KO1>10 & read_count_Setd2KO2>10 & read_count_Setd2KO3>10 & read_count_Setd2KO4>10 & read_count_Setd2WT1>10 & read_count_Setd2WT2>10 & read_count_Setd2WT3>10),]
#detach(newmatrix_miRNAexp)
newmatrix_miRNAexp <- as.matrix(newmatrix_miRNAexp)
newmatrix_miRNAexp <- log2(newmatrix_miRNAexp+1)

condition <- factor(c("Setd2KO","Setd2KO","Setd2KO","Setd2KO","Setd2WT","Setd2WT","Setd2WT"))
pheno=data.frame("sample"=colnames(newmatrix_miRNAexp),"sample_id"=factor(c(1,2,3,4,5,6,7)),"condition"=factor(c(rep("Setd2KO",4),rep("Setd2WT",3))),"batch"=factor(c(1,1,2,2,1,1,2)))
batch = pheno$batch
mod = model.matrix(~as.factor(condition)+as.factor(batch), data=pheno)
mod0 = model.matrix(~as.factor(batch),data=pheno)
n.sv = num.sv(newmatrix_miRNAexp,mod,method="leek")
svseq <- svaseq(newmatrix_miRNAexp, mod, mod0, n.sv = 3)
colData <- data.frame(row.names=colnames(newmatrix_miRNAexp), condition)
dds <- DESeqDataSetFromMatrix(newmatrix_miRNAexp, pheno, design= ~ condition+batch)

dds$SV1 <- svseq$sv[,1]
dds$SV2 <- svseq$sv[,2]
dds$SV3 <- svseq$sv[,3]

design(dds) <- as.formula(paste("~ SV1 + SV2 + SV3", "+ condition + batch"))
dds <- DESeq(dds)
res <- results(dds,contrast=c("condition","Setd2KO","Setd2WT"))

mcols(res, use.names = TRUE)
summary(res)
res <- res[order(res$padj),]
resdata <- merge (as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata, file="CYC_Setd2KO_vs_Setd2WT.csv",row.names = FALSE)

svanoval_DESeq2_rmnoval_morethan10 <- read.xlsx("svanovalin_DESeq2_rmnoval_morethan10.xlsx")
svanoval_DESeq2_rmnoval_morethan10$padj <- p.adjust(svanoval_DESeq2_rmnoval_morethan10$pvalue,method = "fdr")
write.csv(svanoval_DESeq2_rmnoval_morethan10,"svanoval_DESeq2_rmnoval_morethan10_FDRadj.csv",row.names = FALSE)

svanoval_DESeq2_rmnoval_morethan10 <- read.csv("svanoval_DESeq2_rmnoval_morethan10_FDRadj.csv",header = TRUE)
svanoval_DESeq2_rmnoval_morethan10$padj <- p.adjust(svanoval_DESeq2_rmnoval_morethan10$pvalue,method = "fdr")
write.csv(svanoval_DESeq2_rmnoval_morethan10,"svanoval_DESeq2_rmnoval_morethan10_FDRadj.csv",row.names = FALSE)
