import sys
import pandas as pd

dir=sys.argv[1]
df_cardiac = pd.read_csv(sys.argv[2],header=None)
out_file = sys.argv[3]
df_cardiac.columns = ['chr','start','end','gene','refseq']

df_main = pd.DataFrame()
cols=[]

for indx,row in df_cardiac.iterrows():

    chr=row['chr']
    gene=row['gene']
    df=pd.read_csv(dir+"."+chr+"_"+gene+".vcf",sep='\t')
    df = df[~df.iloc[:,12].isnull()]
    df['AL1.AD'] = [x.split(',')[0] for x in df.iloc[:,12].values]
    df['AL2.AD'] = [x.split(',')[1] for x in df.iloc[:,12].values]


    df_main = df_main.append(df,ignore_index=True)

    if chr =="chr1" and gene=="ACTN2":
        cols=df.columns

df_main.columns=cols
df_main['POS'] = df_main['POS'].astype(int)
df_main.to_csv(out_file,sep='\t',index=False,header=False)
