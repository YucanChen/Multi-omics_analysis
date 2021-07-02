import os, sys, xlrd, xlwt, linecache, re, operator
import pandas as pd
import numpy as np

gencodevM25annotation=linecache.getlines("./gencode.vM25_proteincoding2.gff3")
new_gencodevM25annotation=[]
for i in range(len(gencodevM25annotation)):
    a=gencodevM25annotation[i].split('\t')
    new_gencodevM25annotation.append(a)
new_gencodevM25annotation=pd.DataFrame(new_gencodevM25annotation)
new_gencodevM25annotation.columns=["chr","source","type","start","end","no1","strand","no2","details1","gene_name","details2"]
new_gencodevM25annotation=new_gencodevM25annotation[new_gencodevM25annotation["type"]=="gene"]

H3K36me3_peak=pd.read_csv("./diffbind_allDB.csv")
H3K36me3_peak=H3K36me3_peak[H3K36me3_peak["Fold"]<=0]

new_gencodevM25annotation['start'] = new_gencodevM25annotation['start'].astype('int')
new_gencodevM25annotation['end'] = new_gencodevM25annotation['end'].astype('int')
H3K36me3_peak['start'] = H3K36me3_peak['start'].astype('int')
H3K36me3_peak['end'] = H3K36me3_peak['end'].astype('int')

H3K36me3_peak["center"]=(H3K36me3_peak["end"]+H3K36me3_peak["start"])/2
new_gencodevM25annotation["feature_length"]=(new_gencodevM25annotation["end"]-new_gencodevM25annotation["start"]+1)
new_gencodevM25annotation["center"]=(new_gencodevM25annotation["end"]+new_gencodevM25annotation["start"])/2

positive_list = []
mapping_list = []
for i in range(len(H3K36me3_peak)):
    test = H3K36me3_peak.iloc[i:(i+1),:]
    for row in test.itertuples():
        subgene = new_gencodevM25annotation[(new_gencodevM25annotation["chr"] == getattr(row, 'seqnames'))]
        Distance = abs(subgene["center"] - getattr(row, 'center'))
        sum_radius = (subgene["feature_length"] + getattr(row, 'width')) / 2
        compare1 = sum_radius - Distance
        compare1=list(compare1)
        positive_list = [n for n in compare1 if n > 0]
        Distance=list(Distance)
        compare2 = [n for n in Distance if abs(n) < 50000]
        if len(positive_list) > 0:
            for i in range(len(compare1)):
                if compare1[i] > 0:
                    mapping_list.append(str(list(row)).replace('[','').replace(']','').replace('\'','')+','+str(list(subgene.iloc[i,:])).replace('[','').replace(']','').replace('\\n','').replace('\'',''))
        if len(compare2) > 0:
            for i in range(len(Distance)):
                if abs(Distance[i]) < 50000:
                    mapping_list.append(str(list(row)).replace('[','').replace(']','').replace('\'','')+','+str(list(subgene.iloc[i,:])).replace('[','').replace(']','').replace('\\n','').replace('\'',''))
        if (len(positive_list) == 0) & (len(compare2) == 0):
            mapping_list.append(str(list(row)).replace('[','').replace(']','').replace('\'','')+','+"None")

mapping_list = list(np.unique(mapping_list))
               
f=open("./Part4_mRNAprotein_H3K36me3peakmapping50kb_list.txt","a") #add within 50kb
f.write(str(0)+','+str(list(H3K36me3_peak.columns)).replace('[','').replace(']','').replace('\'','')+','+str(list(new_gencodevM25annotation.columns)).replace('[','').replace(']','').replace('\'','')+'\n')
for i in range(len(mapping_list)):
    f.write(mapping_list[i] + '\n')
f.close()
