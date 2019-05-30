import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import warnings
warnings.filterwarnings("ignore")

##################################
### input variables
##################################

SAMPLE=sys.argv[1]
OUT_DIR=sys.argv[2]
ENV=sys.argv[3]
DEPTH=sys.argv[4]
NALF=sys.argv[5]
TALF=sys.argv[6]

##################################
### parameters
##################################




def calcSum(row):
    vals= row.split(',')
    return int(vals[0])+int(vals[1])

def calcADJ_AF(normal_af,tumor_af):
    if tumor_af > normal_af :
        return 1.0-(normal_af/tumor_af)
    else:
        return 0.0

def calcADJ_QSS(tumor_qss,tumor_depth):
    if tumor_depth > 0.0 :
        return tumor_qss/tumor_depth
    else:
        return 0.0




df_main = pd.DataFrame()
cols=[]
for chr in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']:
    df=pd.read_csv(OUT_DIR+"chr"+chr+"_gatk.output.vcf.filter.txt",sep='\t')
    df['NORMAL.DP'] = df['NORMAL.AD'].apply(calcSum)
    df['TUMOR.DP'] = df['TUMOR.AD'].apply(calcSum)
    # df = df[(df['NORMAL.DP']>=20) & (df['TUMOR.DP']>=20 )]
    df_main = df_main.append(df,ignore_index=True)
    if chr =='1':
        cols=df.columns

df_main.columns=cols
df_main.to_csv(OUT_DIR+SAMPLE+".mutect.combine.vcf.txt",sep='\t',index=False)

df_both= pd.read_csv(OUT_DIR+SAMPLE+".mutect.combine.vcf.txt",sep='\t')

print(df_both.shape)

###filtering:
filter_indx = df_both[df_both['CHROM'].str.contains("GL|gl")].index
df_both.drop(filter_indx, inplace=True)

filter_indx = df_both[df_both['CHROM'].str.contains("chrM")].index
df_both.drop(filter_indx, inplace=True)


print('After GL/gl:'+ str(df_both.shape[0]))

## add normal and tumor depth
df_both['NORMAL.DP'] = df_both['NORMAL.AD'].apply(calcSum)
df_both['TUMOR.DP'] = df_both['TUMOR.AD'].apply(calcSum)
df_both['ADJ_AF'] = [calcADJ_AF(x[0],x[1]) for x in zip(df_both['NORMAL.AF'],df_both['TUMOR.AF'])]
df_both['ADJ_QSS'] = [calcADJ_QSS(x[0],x[1]) for x in zip(df_both['TUMOR.QSS'].str.replace(',','' ).astype(float),df_both['TUMOR.DP'])]

ADJ_QSS_FILTER=30
NLOD_FILTER=10
TLOD_FILTER=10
NORMAL_DEPTH=int(DEPTH)
TUMOR_DEPTH=int(DEPTH)
NORMAL_ALT_ALLELE_FREQ=(float(NALF)/100.0)
TUMOR_ALT_ALLELE_FREQ=(float(TALF)/100.0)
NORMAL_STRAND_BIAS_REF_F1R2=0
NORMAL_STRAND_BIAS_REF_F2R1=0
TUMOR_STRAND_BIAS_ALT_F1R2=0
TUMOR_STRAND_BIAS_ALT_F1R2=0


print( "Before filter", df_both.shape)
df_both =df_both[ ( \
                     (df_both['NORMAL.DP'] >= NORMAL_DEPTH ) & \
                     (df_both['TUMOR.DP'] >= TUMOR_DEPTH ) & \
                     (df_both['NORMAL.AF'] <= NORMAL_ALT_ALLELE_FREQ ) & \
                     (df_both['TUMOR.AF'] >= TUMOR_ALT_ALLELE_FREQ ) & \
                     (df_both['NLOD'] >= NLOD_FILTER ) & \
                     (df_both['TLOD'] >= NLOD_FILTER ) & \
                     (df_both['NORMAL.REF_F1R2'] >= NORMAL_STRAND_BIAS_REF_F1R2 ) & \
                     (df_both['NORMAL.REF_F2R1'] >= NORMAL_STRAND_BIAS_REF_F2R1 ) & \
                     (df_both['TUMOR.ALT_F1R2'] >= TUMOR_STRAND_BIAS_ALT_F1R2 ) & \
                     (df_both['TUMOR.ALT_F2R1'] >= TUMOR_STRAND_BIAS_ALT_F1R2 ) & \
                     (df_both['ADJ_QSS'] >= ADJ_QSS_FILTER ) \
                )]

print( "After filter", df_both.shape)

# df_both.to_csv("01_gatk_filter_combine.vcf.txt",sep='\t',index=False)


##################################
### format for VEP
##################################
vcf=[]
for indx,row in df_both.iterrows():
    info=str(row['NORMAL.DP'])+';'+str(row['TUMOR.DP'])
    vcf.append([row['CHROM'],row['POS'],'.',row['REF'],row['ALT'],'.','.',info])

df_vcf=pd.DataFrame(vcf)
df_vcf.columns=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']
df_vcf.to_csv(OUT_DIR+SAMPLE+".mutect."+DEPTH+"_"+NALF+"_"+TALF,sep='\t',index=False)
