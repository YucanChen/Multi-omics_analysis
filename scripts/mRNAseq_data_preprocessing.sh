# Code chunk 1: quality control
## Fastqc
cd /GEOsample_mRNA/
ls GEOmRNA_rawdata|while read id;do gunzip GEOmRNA_rawdata/$id;done
ls GEOmRNA_rawdata|while read id;do fastqc -o fastqc_test -t 4 GEOmRNA_rawdata/$id;done
# generate reports
multiqc fastqc_test -o fastqc_test -i rawdata

## Cut adaptors: Trim galore
cd /GEOsample_mRNA/GEOmRNA_rawdata
ls /GEOsample_mRNA/GEOmRNA_rawdata/*_R1_001.fastq >1
ls /GEOsample_mRNA/GEOmRNA_rawdata/*_R2_001.fastq >2
paste 1 2  > config

#!/bin/bash
cd /GEOsample_mRNA/GEOmRNA_rawdata
cat config |while read id
do
   arr=(${id})
   fq1=${arr[0]}
   fq2=${arr[1]} 
   nohup trim_galore --quality 25 --phred33 --paired --output_dir /GEOsample_mRNA/clean_data --length 75 --stringency 5 $fq1 $fq2
   done
exit 0

## Re-fastqc
#!/bin/bash
cd /GEOsample_mRNA/clean_data/
ls *.fq|while read id;
do
   fastqc -o /GEOsample_mRNA/clean_fastqc_test $id;
   done
multiqc /GEOsample_mRNA/clean_fastqc_test -o /GEOsample_mRNA/clean_fastqc_test/multiqc_data -i clean_data
exit 0

# Code chunk 2: genome mapping
mkdir ~/reference/mm10/index/hisat2
cd ~/reference/mm10/index/hisat2
wget -O mm10.tar.gz https://cloud.biohpc.swmed.edu/index.php/s/mm10/download
tar -zxvf *.tar.gz

## Build Index: Stringtie
extract_splice_sites.py ~/reference/mm10/gtf/Gencode/gencode.vM25.annotation.gtf > ~/reference/mm10/gtf/Gencode/gencode.vM25.annotation.ss
extract_exons.py ~/reference/mm10/gtf/Gencode/gencode.vM25.annotation.gtf > ~/reference/mm10/gtf/Gencode/gencode.vM25.annotation.exon

#!/bin/bash
cd ~/reference/mm10/gtf/Gencode/
module load hisat2/2.1.0-gcc-8.3.0
hisat2-build --ss gencode.vM25.annotation.ss --exon gencode.vM25.annotation.exon ~/reference/mm10/genome.fa gencode.vM25.annotation_tran
exit 0

cd /GEOsample_mRNA/aligned/
ls /GEOsample_mRNA/clean_data/*R1_001_val_1.fq >1 
ls /GEOsample_mRNA/clean_data/*R2_001_val_2.fq >2
paste 1 2  > config

## HISAT2 mm10 mapping (UCSC)
cd /GEOsample_mRNA/aligned/
ls /GEOsample_mRNA/clean_data/*R1_001_val_1.fq >1 
ls /GEOsample_mRNA/clean_data/*R2_001_val_2.fq >2
paste 1 2  > config

#!/bin/bash
cd /GEOsample_mRNA/aligned/
cat config |while read id
do
   arr=(${id})
   fq1=${arr[0]}
   fq2=${arr[1]}
   hisat2 --phred33 -p 8 --dta -x ~/reference/mm10/gtf/Gencode/gencode.vM25.annotation_tran -1 $fq1 -2 $fq2 -S `basename $fq1 R1_001_val_1.fq`hisat2.sam;
done

## Sam2sorted.Bam
ls *.sam |while read id; do samtools sort -@ 8 -o /GEOsample_mRNA/sorted_bam/`basename $id .sam`.bam $id; done
cd /GEOsample_mRNA/sorted_bam/
ls *.bam |while read id; do stringtie -p 8 -G ~/reference/mm10/gtf/Gencode/gencode.vM25.annotation.gtf -o /GEOsample_mRNA/stringtie/`basename $id .bam`.gtf $id; done
exit 0

#!/bin/bash
cd /GEOsample_mRNA/stringtie/
stringtie --merge -p 8 -G ~/reference/mm10/gtf/Gencode/gencode.vM25.annotation.gtf -o stringtie_merged.gtf mergelist.txt
ls /GEOsample_mRNA/sorted_bam/*.bam |while read id; do stringtie -e -B -p 8 -G stringtie_merged.gtf -o /GEOsample_mRNA/stringtie/ballgown/SETD2`basename $id _hisat2.bam`/`basename $id _hisat2.bam`_ballgown.gtf $id; done
exit 0

#!/bin/bash
cd /GEOsample_mRNA/stringtie_unmerge/
ls /GEOsample_mRNA/sorted_bam/*.bam |while read id; 
do 
   stringtie -p 6 -G ~/reference/mm10/gtf/Gencode/gencode.vM25.annotation.gtf -o /GEOsample_mRNA/stringtie_unmerge/ballgown/SETD2`basename $id _hisat2.bam`/`basename $id _hisat2.bam`_ballgown.gtf -A `basename $id _hisat2.bam`.tab -B -e -l `basename $id _hisat2.bam` $id;
  done
exit 0

# Code chunk 3: quality control of bam files (samtools flagstat)
cd /GEOsample_mRNA/sorted_bam/
ls *hisat2sorted.bam|while read id;
do 
   samtools flagstat $id > `basename $id sorted.bam`sorted.flagstat
   done
