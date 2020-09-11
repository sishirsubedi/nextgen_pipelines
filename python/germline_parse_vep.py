import sys
import pandas as pd

VEP_FILE=sys.argv[1]
OUT_FILE=sys.argv[2]

header=[]
with open(VEP_FILE) as myfile:
    header = [next(myfile) for x in range(3)]
columns_line=str(header[2].strip().split(':')[1]).split('|')

df_vep = pd.read_csv(VEP_FILE,sep='\t',skiprows=4,header=None)
df_vep.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
df_vep["HGMD_INFO"]= [x.split(";")[0] for x in df_vep.INFO]
df_vep["VEP_INFO"]= [x.split(";")[1] for x in df_vep.INFO]

df_vep = df_vep[['CHROM', 'POS', 'REF', 'ALT','VEP_INFO']]

rowwise=[]
for indx,row in df_vep.iterrows():
    transcripts = row.VEP_INFO.split(',')
    for transcript in transcripts:
            transcript_info = []
            transcript_info.append(row.CHROM)
            transcript_info.append(row.POS)
            transcript_info.append(row.REF)
            transcript_info.append(row.ALT)
            for t_item in transcript.split('|'):
                transcript_info.append(t_item)
            rowwise.append(transcript_info)
new_columns = ['CHROM', 'POS', 'REF', 'ALT']
for x in columns_line: new_columns.append(x.replace(" ","").replace('">',''))

df_final = pd.DataFrame(rowwise,columns=new_columns)

df_final = df_final[['CHROM', 'POS', 'REF', 'ALT',
       'Feature',
       'CANONICAL',
       'SIFT',
       'PolyPhen',
       'HGVSc',
       'HGVSp']]
df_final.to_csv(OUT_FILE,sep="\t",index=False)
