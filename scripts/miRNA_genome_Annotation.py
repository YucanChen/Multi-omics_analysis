import os,sys,xlrd,xlwt,linecache, re, operator
import pandas as pd
import numpy as np

mmu_premiRNA=linecache.getlines("./mmu_premiRNA2.gff3")
gencodevM25annotation=linecache.getlines("./gencode.vM25_proteincoding2.gff3")
new_mmu_premiRNA=[]
for i in range(len(mmu_premiRNA)):
    a=mmu_premiRNA[i].replace("\n", "").split('\t')
    new_mmu_premiRNA.append(a)
new_mmu_premiRNA=pd.DataFrame(new_mmu_premiRNA)
new_gencodevM25annotation=[]
for i in range(len(gencodevM25annotation)):
    a=gencodevM25annotation[i].split('\t')
    new_gencodevM25annotation.append(a)
new_gencodevM25annotation=pd.DataFrame(new_gencodevM25annotation)

new_gencodevM25annotation.columns=["chr","source","type","start","end","no1","strand","no2","details1","gene_name","details2"]
new_mmu_premiRNA.columns=["chr","no1","info","start","end","no2","strand","no3","details","names"]
new_gencodevM25annotation['start'] = new_gencodevM25annotation['start'].astype('int')
new_gencodevM25annotation['end'] = new_gencodevM25annotation['end'].astype('int')
new_mmu_premiRNA['start'] = new_mmu_premiRNA['start'].astype('int')
new_mmu_premiRNA['end'] = new_mmu_premiRNA['end'].astype('int')

new_mmu_premiRNA["miRNA_length"]=new_mmu_premiRNA["end"]-new_mmu_premiRNA["start"]
new_mmu_premiRNA["center"]=(new_mmu_premiRNA["end"]+new_mmu_premiRNA["start"])/2
new_gencodevM25annotation["feature_length"]=new_gencodevM25annotation["end"]-new_gencodevM25annotation["start"]
new_gencodevM25annotation["center"]=(new_gencodevM25annotation["end"]+new_gencodevM25annotation["start"])/2

gencodevM25gene=new_gencodevM25annotation[new_gencodevM25annotation["type"]=="gene"]
gencodevM25exon=new_gencodevM25annotation[new_gencodevM25annotation["type"]=="exon"]

positive_list = []
type_list = []
gene_list = []
for i in range(820,1227):
    test = new_mmu_premiRNA.iloc[i:(i+1),:]
    for row1 in test.itertuples():
        subexon = gencodevM25exon[(gencodevM25exon["chr"] == getattr(row1, 'chr')) & (gencodevM25exon["strand"] == getattr(row1, 'strand'))]
        Distance = abs(subexon["center"] - getattr(row1, 'center'))
        sum_radius = (subexon["feature_length"] + getattr(row1, 'miRNA_length')) / 2
        compare1 = sum_radius - Distance
        compare1=list(compare1)
        positive_list = [n for n in compare1 if n > 0]
        if len(positive_list) > 0:
            type_list.append("exonic")
            gene_list.append(subexon.iloc[compare1.index(max(compare1)),9].replace('\n',''))
        else:
            for row2 in test.itertuples():
                subgene = gencodevM25gene[(gencodevM25gene["chr"] == getattr(row2, 'chr')) & (gencodevM25gene["strand"] == getattr(row2, 'strand'))]
                Distance = abs(subgene["center"] - getattr(row2, 'center'))
                sum_radius = (subgene["feature_length"] + getattr(row2, 'miRNA_length')) / 2
                compare2 = sum_radius - Distance
                compare2=list(compare2)
            positive_list = [n for n in compare2 if n > 0]
            if len(positive_list) > 0:
                type_list.append("intronic")
                gene_list.append(subgene.iloc[compare2.index(max(compare2)),9].replace('\n',''))
            else:
                type_list.append("intergenic")
                gene_list.append("None")

mmu_premiRNA_list3=new_mmu_premiRNA.loc[820:1227,["chr","start","end","details","names","strand"]]

f=open("./pre_miRNA_type_list3_new.txt","w")
for i in range(len(type_list)):
    f.write(str(list(mmu_premiRNA_list3.iloc[i,:])).replace('[','').replace(']','').replace(', ','\t').replace("'",'')+ '\t' +type_list[i]+ '\t' + gene_list[i]+'\n')
f.close()
