#!/bin/bash

peakcover_TSS=peak_covered_withTSS.bed
nopeakcover_TSS=no_peakcover_withTSS.bed
miRNATSS_peakcover=miRNA_TSS_peakcover_forplot.bed
miRNATSS_nopeakcover=miRNA_TSS_nopeakcover_forplot.bed

computeMatrix reference-point -p 15 --referencePoint center -b 1000 -a 1000 --binSize 1 -R $miRNATSS_peakcover $peakcover_TSS $miRNATSS_nopeakcover $nopeakcover_TSS -S WT1_H3K36me3_log2ratio.bw WT2_H3K36me3_log2ratio.bw KO1_H3K36me3_log2ratio.bw --outFileName Deeptools_Results/Combined_matrix_regionclass.gz &
computeMatrix reference-point -p 15 --referencePoint center -b 1000 -a 1000 --binSize 1 -S WT1_H3K36me3_log2ratio.bw WT2_H3K36me3_log2ratio.bw KO1_H3K36me3_log2ratio.bw -R $miRNATSS_peakcover $peakcover_TSS $miRNATSS_nopeakcover $nopeakcover_TSS --outFileName Deeptools_Results/Combined_matrix_sampleclass.gz &
wait
exit 0
