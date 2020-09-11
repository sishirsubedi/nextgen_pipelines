import sys
import pandas as pd

def createInfo (df):
    info_column = []
    for indx,row in df.iterrows():
        info={}
        info["DP"]=row[8]
        info["GT"]=row[10]
        info["DP.REF"]=row[15]
        info["DP.ALT"]=row[16]
        info["ALT_TRANSCRIPT"]= str(row[18]) +":"+str(row[19])
        info["GENE"]= row[20]
        info["HGMD_TRANSCRIPT"]= row[21]
        info["STRAND"]= row[22]
        info["TOTAL_EXONS"]= row[23]
        info["GENOMIC_START_END"]= str(row[24]) +":"+str(row[25])
        info["CDS_START_END"]= str(row[26]) +":"+str(row[27])
        info["ALL_TRANSCRIPT_START_END"]= str(row[28]) +":"+str(row[29])
        # transcript_starts = [int(x) for x in row[28].split(",") if x != '']
        # info["ALT_EXON"]= transcript_starts.index(row[18])+1
        info_column.append(info)
    return info_column

infile=sys.argv[1]
outfile = sys.argv[2]
df = pd.read_csv(infile,sep='\t',header=None)
# df = pd.read_csv("NA12878-N.mutect.combine.v3.vcf.txt",sep='\t',header=None)

info_column = createInfo(df)
df = df.iloc[:,[0,1,3,4,5,7,6]]
df["INFO"] = info_column
df.to_csv(outfile,sep='\t',index=False,header=False)
