# Code chunk 1 ：Fastqc (Result: raw clean reads)
fastqc -o fastqc_test -t 4 colon_ChIP_S5_R1_001.fastq
fastqc -o fastqc_test -t 4 colon_Input_S10_R1_001.fastq

#####batch task#####
# Code chunk 2: fastq2sam (mem)
#!/bin/bash
cd ./Rep2_WT_ChIP_SE_fastq/
ls *.fastq |while read id; do (bwa mem -t 16 -M ~/reference/index/bwa/mm10.fa $id > ./Rep2_WT_ChIP_SE_mapped_sam/`basename $id .fastq`.sam) 2>&1|tee `basename $id .fastq`.sam.tmp && tail `basename $id .fastq`.sam.tmp > ../Rep2_WT_ChIP_SE_mapped_sam/`basename $id .fastq`.sam.bwa_mem.log && rm `basename $id .fastq`.sam.tmp;done
exit 0

#  Code chunk 3: sam2bam (bwa)
#!/bin/bash
cd ./Rep2_WT_ChIP_SE_mapped_sam/
ls *.sam |while read id;do(samtools view -h -bS -F 4 -q 30 -@ 16 $id |samtools sort -@ 10 -O BAM > ./Rep2_WT_ChIP_SE_marked_sorted_bam/`basename $id .sam`.bam); done
cd ./Rep2_WT_ChIP_SE_marked_sorted_bam/
ls *.bam|while read id;do(java -jar ~/picard-2.22.0-0/picard.jar MarkDuplicates \
I=$id \
O=`basename $id .bam`.marked_duplicates.bam \
M=`basename $id .bam`.marked_dup_metrics.txt && samtools flagstat `basename $id .bam`.marked_duplicates.bam > `basename $id .bam`.sorted_bam_stat); done
ls *.marked_duplicates.bam | while read id; do(java -jar ~/picard-2.22.0-0/picard.jar MarkDuplicates \
I=$id \
O=./Rep2_WT_ChIP_SE_rmDup_sorted_bam/`basename $id .marked_duplicates.bam`.rmDup.bam \
M=./Rep2_WT_ChIP_SE_rmDup_sorted_bam/`basename $id .marked_duplicates.bam`.rmDup_metrics.txt REMOVE_DUPLICATES=true && samtools flagstat ./Rep2_WT_ChIP_SE_rmDup_sorted_bam/`basename $id .marked_duplicates.bam`.rmDup.bam>../Rep2_WT_ChIP_SE_rmDup_sorted_bam/`basename $id .marked_duplicates.bam`.rmDup_sorted_bam_stat); done
exit 0

#  Code chunk 4 ：bam call peaks
## <macs2 Call peaks> WT Rep2 example
#!/bin/bash
cd ./Rep2_WT_ChIP_SE_rmDup_sorted_bam/
macs2 callpeak -t colon_ChIP_S5_R1_001.rmDup.bam -c colon_Input_S10_R1_001.rmDup.bam -f BAM -g mm -n Rep2WT_H3K36me3_S5_R1_broad01 -B -q 0.05 --broad --nomodel --scale-to large --outdir ./Rep2_WT_ChIP_SE_H3K36me3_macs2callpeaks_broad01
exit 0

# rmdup bam to bai index
ls *.rmDup.bam| while read id;do samtools index -b $id;done
