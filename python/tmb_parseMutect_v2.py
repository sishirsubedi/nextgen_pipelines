import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.stats as st
import pandas as pd
import seaborn as sns
import numpy as np
import warnings
warnings.filterwarnings("ignore")

##################################
### input variables
##################################

SAMPLE="COLO-829-T_COLO-829BL-N"
OUT_DIR="/home/environments/ngs_test/nextseqAnalysis/tmbAssay/190128_NS500761_0287_AHLW2JBGX9/COLO-829-T/Paired/COLO-829-T_COLO-829BL-N/mutect/"
ENV="test"
DEPTH="10"
NALF="10"
TALF="10"

##################################
### parameters
##################################




def calcDepthSum(row):
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
    df['NORMAL.DP'] = df['NORMAL.AD'].apply(calcDepthSum)
    df['TUMOR.DP'] = df['TUMOR.AD'].apply(calcDepthSum)
    df_main = df_main.append(df,ignore_index=True)
    if chr =='1':
        cols=df.columns

df_main.columns=cols
df_main.to_csv(OUT_DIR+SAMPLE+".mutect.combine.vcf.txt",sep='\t',index=False)

df_both= pd.read_csv(OUT_DIR+SAMPLE+".mutect.combine.vcf.txt",sep='\t')

# print(df_both.shape)

###filtering:
filter_indx = df_both[df_both['CHROM'].str.contains("GL|gl")].index
df_both.drop(filter_indx, inplace=True)

filter_indx = df_both[df_both['CHROM'].str.contains("chrM")].index
df_both.drop(filter_indx, inplace=True)


# print('After GL/gl:'+ str(df_both.shape[0]))

## adjustments for filtering
df_both['ADJ_AF'] = [calcADJ_AF(x[0],x[1]) for x in zip(df_both['NORMAL.AF'],df_both['TUMOR.AF'])]
# df_both['ADJ_QSS'] = [calcADJ_QSS(x[0],x[1]) for x in zip(df_both['TUMOR.QSS'].str.replace(',','' ).astype(float),df_both['TUMOR.DP'])]


def expoCutVal(df_list):
    mn = df_list.mean()
    rate = 1./mn
    var = (1./rate)**2
    return [mn,rate, var]


NLOD_FILTER=2.55
TLOD_FILTER=7.85
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
                     (df_both['TLOD'] >= TLOD_FILTER ) & \
                     (df_both['NORMAL.REF_F1R2'] >= NORMAL_STRAND_BIAS_REF_F1R2 ) & \
                     (df_both['NORMAL.REF_F2R1'] >= NORMAL_STRAND_BIAS_REF_F2R1 ) & \
                     (df_both['TUMOR.ALT_F1R2'] >= TUMOR_STRAND_BIAS_ALT_F1R2 ) & \
                     (df_both['TUMOR.ALT_F2R1'] >= TUMOR_STRAND_BIAS_ALT_F1R2 )
                )]

print( "After filter", df_both.shape)
#
# # df_both.to_csv("01_gatk_filter_combine.vcf.txt",sep='\t',index=False)
#
#
# ##################################
# ### format for VEP
# ##################################
# vcf=[]
# for indx,row in df_both.iterrows():
#     info=str(row['NORMAL.DP'])+';'+str(row['TUMOR.DP'])
#     vcf.append([row['CHROM'],row['POS'],'.',row['REF'],row['ALT'],'.','.',info])
#
# df_vcf=pd.DataFrame(vcf)
# df_vcf.columns=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']
# df_vcf.to_csv(OUT_DIR+SAMPLE+".mutect."+str(DEPTH)+"_"+str(NALF)+"_"+str(TALF),sep='\t',index=False)
#
#
# df_ref = df_vcf
#
# design= "/home/hhadmin/exome_pipeline/agilentCre/cre_design.bed"
# df_design = pd.read_csv(design,sep='\t')
# df_design.columns=['CHROM','START','END']
#
#
# filter=[]
# for indx,row in df_ref.iterrows():
#     if ( df_design[ (df_design['CHROM'] == row['CHROM']) & (df_design['START'] <= row['POS']) & (df_design['END'] >= row['POS']) ].shape[0]>=1):
#         filter.append([row['CHROM'],row['POS'],row['ID'],row['REF'],row['ALT'],row['QUAL'],row['FILTER'],row['INFO']])
#
# df_filter=pd.DataFrame(filter)
# df_filter.columns=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']
# df_filter.to_csv(OUT_DIR+SAMPLE+".mutect."+str(DEPTH)+"_"+str(NALF)+"_"+str(TALF)+".crefilter",sep='\t',index=False)
#
#
#
#
#
# OUT_DIR="/home/environments/ngs_test/nextseqAnalysis/tmbAssay/190128_NS500761_0287_AHLW2JBGX9/COLO-829-T/Paired/"
# file1 = OUT_DIR+SAMPLE+"/varscan/"+SAMPLE+".varscan."+DEPTH+"_"+NALF+"_"+TALF+".crefilter"
# df_f1 = pd.read_csv(file1,sep='\t')
# file2 = OUT_DIR+SAMPLE+"/mutect/"+SAMPLE+".mutect."+DEPTH+"_"+NALF+"_"+TALF+".crefilter"
# df_f2 = pd.read_csv(file2,sep='\t')
# file3 = OUT_DIR+SAMPLE+"/strelka/"+SAMPLE+".strelka."+DEPTH+"_"+NALF+"_"+TALF+".crefilter"
# df_f3 = pd.read_csv(file3,sep='\t')
#
#
#
# df_all = pd.concat([df_f1,df_f2,df_f3])
# df_all.drop_duplicates(subset=['CHROM','POS','ID','REF','ALT'],keep='first',inplace=True)
# df_all.reset_index(inplace=True)
# # print("Unique variants-",df_all.shape)
#
# def addVCInfo(row,df):
#     if(df[ ( \
#             (df['CHROM']==row['CHROM']) & \
#             (df['POS']==row['POS']) & \
#             (df['REF']==row['REF']) & \
#             (df['ALT']==row['ALT']) )].shape[0] == 1):
#         return 1
#     else:
#         return 0
#
#
# for indx,row in df_all.iterrows():
#     temp=[]
#     temp.append(addVCInfo(row,df_f1))
#     temp.append(addVCInfo(row,df_f2))
#     temp.append(addVCInfo(row,df_f3))
#     if sum(temp) >=2:
#         df_all.at[indx,'INFO'] += ';'+str(temp[0]) + ';'+str(temp[1]) + ';'+str(temp[2])
#     else:
#         df_all.drop(indx,inplace=True)
#
# df_all['POS'] = df_all['POS'].astype(int)
#
#
# # df_all=df_all[['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']]
# df_snp_mut_1 = df_all[['CHROM','POS','REF','ALT']]
#
# # print("Unique variants 2/3-",df_all.shape)
#
# # print(df_snp_mut_1.head(2))
# df_snp_2 = pd.read_csv("/home/hhadmin/exome_pipeline/RefStd/reference_standard_creDesign.txt",sep='\t')
# # df_snp_2 = pd.read_csv("/home/hhadmin/exome_pipeline/RefStd/colo_old.txt",sep='\t')
#
# df_snp_mut_2 = df_snp_2[['CHROM','POS','REF','ALT']]
# # print(df_snp_mut_2.head(2))
# # print("refstd total:")
# # print(df_snp_mut_2.shape)
#
# df_join = pd.merge(left=df_snp_mut_2, right=df_snp_mut_1, on=['CHROM','POS','REF','ALT'],how='outer',indicator=True)
#
# # df_join.to_csv(OUTPUT_DIR_SAMPLE+'/'+SAMPLE_1+".2VENN.OUT.txt",sep='\t',index=False)
#
# print(df_join[df_join['_merge']=='left_only'].shape[0],df_join[df_join['_merge']=='both'].shape[0],df_join[df_join['_merge']=='right_only'].shape[0])
