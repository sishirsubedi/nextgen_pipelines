import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib.pyplot import figure


def calcTiTvRatio(df):

    AtoG = df[ (df.REF=='A') & (df.ALT=='G') ].shape[0]
    GtoA = df[ (df.REF=='G') & (df.ALT=='A') ].shape[0]
    CtoT = df[ (df.REF=='C') & (df.ALT=='T') ].shape[0] ## very common
    TtoC = df[ (df.REF=='T') & (df.ALT=='C') ].shape[0]

    transition = AtoG + GtoA + CtoT + TtoC

    AtoC = df[ (df.REF=='A') & (df.ALT=='C') ].shape[0]
    CtoA = df[ (df.REF=='C') & (df.ALT=='A') ].shape[0]

    AtoT = df[ (df.REF=='A') & (df.ALT=='T') ].shape[0]
    TtoA = df[ (df.REF=='T') & (df.ALT=='A') ].shape[0]

    GtoC = df[ (df.REF=='G') & (df.ALT=='C') ].shape[0]
    CtoG = df[ (df.REF=='C') & (df.ALT=='G') ].shape[0]

    GtoT = df[ (df.REF=='G') & (df.ALT=='T') ].shape[0]
    TtoG = df[ (df.REF=='T') & (df.ALT=='G') ].shape[0]

    transversions = AtoC + CtoA + AtoT + TtoA + GtoC + CtoG + GtoT + TtoG

    return float(transition/transversions)


def calcIndelRatio(df):
    insertion = df[(df.REF.str.len()==1)& (df.ALT.str.len()>1)].shape[0]
    deletion = df[(df.REF.str.len()>1)& (df.ALT.str.len()==1)].shape[0]

    return float(insertion/deletion)


SAMPLE=sys.argv[1]
OUT_DIR=sys.argv[2]
ENV=sys.argv[3]

VEP_FILE=OUT_DIR+"/"+SAMPLE+".vc.combine.after.vep.vcf.txt"


header=[]
with open(VEP_FILE) as myfile:
    header = [next(myfile) for x in range(3)]
columns_line=str(header[2].strip().split(':')[1]).split('|')

df_vep = pd.read_csv(VEP_FILE,sep='\t',skiprows=4,header=None)
df_vep.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
df_vep[columns_line] = df_vep['INFO'].str.split('|',expand=True)
df_vep.drop(['INFO'],axis=1,inplace=True)
df_vep= df_vep[df_vep['GIVEN_REF']==df_vep['USED_REF']]

costum_line=['ND','TD','VARSCAN','MUTECT','STRELKA','CSQ']
df_vep[costum_line] = df_vep[' Allele'].str.split(';',expand=True)
df_vep.drop([' Allele'],axis=1,inplace=True)


for col in ['IMPACT','Consequence']:
    figure(figsize=(16,8))
    fig=sns.countplot(y=df_vep[col], data=df_vep)
    plt.tight_layout()
    plt.savefig(col)
    plt.close()

df_vep.to_csv(OUT_DIR+"/"+SAMPLE+".vc.combine.parsed.vep.vcf.txt",sep='\t',index=False)


print(calcTiTvRatio(df_vep))
print(calcIndelRatio(df_vep))
