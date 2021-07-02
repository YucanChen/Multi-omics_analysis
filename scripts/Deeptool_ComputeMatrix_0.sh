#!/bin/bash

nohup bamCompare -b1 WT_Rep1_aln.bam -b2 WT_Input_Rep1_aln.bam -o WT1_H3K36me3_log2ratio.bw &
nohup bamCompare -b1 colon_ChIP_Rep2_rmDup.bam -b2 colon_Input_Rep2_rmDup.bam -o WT2_H3K36me3_log2ratio.bw &
nohup bamCompare -b1 Setd2KO_Rep1_aln.bam -b2 Setd2KO_Input_Rep1_aln.bam -o KO1_H3K36me3_log2ratio.bw &
exit 0
