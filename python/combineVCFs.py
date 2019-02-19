import sys
import matplotlib
import pandas as pd
import numpy as np

def addVCInfo(row,df):
    if(df[ ( \
            (df['CHROM']==row['CHROM']) & \
            (df['POS']==row['POS']) & \
            (df['REF']==row['REF']) & \
            (df['ALT']==row['ALT']) )].shape[0] == 1):
        return ";"+str(1)
    else:
        return ";"+str(0)

##################################
### input variables
##################################

SAMPLE=sys.argv[1]
OUT_DIR=sys.argv[2]
ENV=sys.argv[3]

file1 = OUT_DIR+"/varscan/"+SAMPLE+".varscan.filter.vcf.txt"
df_f1 = pd.read_csv(file1,sep='\t')
print(df_f1.shape)
file2 = OUT_DIR+"/mutect/"+SAMPLE+".mutect.filter.vcf.txt"
df_f2 = pd.read_csv(file2,sep='\t')
print(df_f2.shape)
file3 = OUT_DIR+"/strelka/"+SAMPLE+".strelka.filter.vcf.txt"
df_f3 = pd.read_csv(file3,sep='\t')
print(df_f3.shape)

df_all = pd.concat([df_f1,df_f2,df_f3])
print(df_all.shape)
df_all.drop_duplicates(subset=['CHROM','POS','ID','REF','ALT'],keep='first',inplace=True)
print(df_all.shape)
# df_all.to_csv(OUT_DIR+"/temp",sep='\t',index=False)

for indx,row in df_all.iterrows():
    temp = row['INFO']
    temp += addVCInfo(row,df_f1)
    temp += addVCInfo(row,df_f2)
    temp += addVCInfo(row,df_f3)
    df_all.loc[indx,'INFO'] = temp



df_all['POS'] = df_all['POS'].astype(int)
df_all.to_csv(OUT_DIR+"/"+SAMPLE+".vc.combine.before.vep.vcf.txt",sep='\t',index=False,header=False)
