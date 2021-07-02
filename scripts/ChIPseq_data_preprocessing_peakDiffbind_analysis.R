library("biomaRt")
library("DiffBind")
# Diffbind_ColonREP1-REP2_samplesheet.csv: the file that contained all required information for diffbind analysis
IBDDBA <- dba(sampleSheet="Diffbind_ColonREP1-REP2_samplesheet.csv")
IBDDBA <- dba.count(IBDDBA)
dba.peakset(IBDDBA, bRetrieve=TRUE)
IBDDBA <- dba.contrast(IBDDBA, IBDDBA$masks$knockout, IBDDBA$masks$wild_type,
                          "Setd2_knockout", "Setd2_wild-type")
IBDDBA <- dba.analyze(IBDDBA)
IBDDBA.DB <- dba.report(IBDDBA)
pdisplay <- dba.plotPCA(IBDDBA,contrast=1)

# plot drawing
pdf("DiffBind_Plots.pdf")
savewarn <- options("warn")
options(warn=-1)
## PCA plots
dba.plotPCA(IBDDBA,DBA_TREATMENT,label=DBA_ID)
## A PCA plot using only the differentially bound sites (corresponding to Figure 3), using an FDR threshold of 0.05, can be drawn as follows:> 
dba.plotPCA(IBDDBA, contrast=1,label=DBA_ID)
dba.plotPCA(IBDDBA,b3D=T)
## MA plots
par(cex=0.9)
dba.plotMA(IBDDBA)
## Volcano plots
dba.plotVolcano(IBDDBA)
## Boxplots
sum(IBDDBA.DB$Fold<0)
dba.plotBox(IBDDBA)
pvals <- dba.plotBox(IBDDBA, contrast=1, method=IBDDBA$config$AnalysisMethod, 
            th=IBDDBA$config$th, bUsePval=FALSE, 
            bNormalized=TRUE, attribute=DBA_CONDITION, 
            bAll=FALSE, bAllIncreased=FALSE, bAllDecreased=FALSE, 
            bDB=TRUE, bDBIncreased=TRUE, bDBDecreased=TRUE,
            pvalMethod=wilcox.test)
pvals
## Heatmaps
corvals <- dba.plotHeatmap(IBDDBA)
corvals

# output all result peaks
IBDDBA.allDB <- dba.report(IBDDBA, method=DBA_DESEQ2, th=1)
write.csv(IBDDBA.allDB,"diffbind_allDB.csv")
