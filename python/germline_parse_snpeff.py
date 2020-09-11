import sys
import pandas as pd
import ast

def updateRefSeq(df):
    vep_selected = []
    for indx,row in df.iterrows():
        hgmd = ast.literal_eval(row.HGMD_INFO)['HGMD_TRANSCRIPT']
        new_vep = []
        for vep_entry in row.VEP_INFO.split(','):
            if hgmd == vep_entry.split('|')[6]:
                new_vep.append(vep_entry)
                break
        if len(new_vep)==0:
            for vep_entry in row.VEP_INFO.split(','):
                if hgmd.split('.')[0] == vep_entry.split('|')[6].split('.')[0]:
                    new_vep.append(vep_entry)
                    break
            if len(new_vep)==0:
                vep_selected.append("SNPEFF:Not Found")
            else:
                vep_selected.append(new_vep)
        else:
            vep_selected.append(new_vep)
    return [x[0] for x in vep_selected]

def updateAltTranscriptPosition(df):
    ## check if variant lies in +/- 3 bp from start or end of transcript
    variant_position_group = []
    for indx,row in df.iterrows():
        start=int(row["ALT_TRANSCRIPT_START"])
        end=int(row["ALT_TRANSCRIPT_END"])
        pos=int(row["pos"])

        variant_position_start = abs(pos-start)
        variant_position_end = abs(end-start)

        if variant_position_start <= 3:
            if row["STRAND"] == "FORWARD":
                variant_position_group.append("START")
            else:
                variant_position_group.append("END")


        elif variant_position_end <= 3:
            if row["STRAND"] == "FORWARD":
                variant_position_group.append("END")
            else:
                variant_position_group.append("START")

        else:
            variant_position_group.append("MIDDLE")

    return variant_position_group


def updateNegativeStrandStartEndPosition(df):
    ## REVERSE start or end of transcript FOR NEG STRAND
    for indx,row in df.iterrows():
        if row["STRAND"] == "REVERSE":
            start=int(row["ALT_TRANSCRIPT_START"])
            end=int(row["ALT_TRANSCRIPT_END"])
            df.ix[indx,"ALT_TRANSCRIPT_START"] = end
            df.ix[indx,"ALT_TRANSCRIPT_END"] = start
    return df



VEP_FILE=sys.argv[1]
OUT_FILE=sys.argv[2]
OUT_FILE_ALL=sys.argv[3]

header=[]
with open(VEP_FILE) as myfile:
    header = [next(myfile) for x in range(3)]
columns_line=str(header[2].strip().split(':')[1]).split('|')

df_snpeff = pd.read_csv(VEP_FILE,sep='\t',skiprows=5,header=None)
df_snpeff.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
df_snpeff["HGMD_INFO"]= [x.split(";")[0] for x in df_snpeff.INFO]
df_snpeff["VEP_INFO"]= [x.split(";")[1] for x in df_snpeff.INFO]

#### select HGMD transcript and expand SNPEFF annotation
df_snpeff["VEP_SELECTED"]= updateRefSeq(df_snpeff)
df_snpeff[columns_line] = df_snpeff['VEP_SELECTED'].str.split('|',expand=True)

#### select expand HGMD/DEPTH annotation
df_snpeff["HGMD_INFO2"] = df_snpeff["HGMD_INFO"].map(lambda d :ast.literal_eval(d))
df_snpeff = df_snpeff.join(pd.DataFrame( df_snpeff["HGMD_INFO2"].to_dict()).T)
df_snpeff.drop([ 'INFO','HGMD_INFO', 'VEP_INFO', 'VEP_SELECTED','HGMD_INFO2'],axis=1,inplace=True)

### adjust columns to match as VEP output
### add alt frequecy and type(snp/indel)
df_snpeff["altFreq"] = [round((float(x)/float(y))*100,2) for x,y in zip(df_snpeff["DP.ALT"],df_snpeff.DP)]
df_snpeff["type"] = [ "snp" if len(x)==1 and len(y)==1 else "indel" for x,y in zip(df_snpeff.REF,df_snpeff.ALT)]
df_snpeff.columns = ['chr', 'pos', 'ID', 'ref', 'alt', 'quality', 'FILTER', 'Allele',
       'consequence', 'impact', 'Gene_Name', 'gene',
       'Feature_Type', 'Feature_ID', 'Transcript_BioType', 'exon',
       'HGVSc', 'HGVSp', 'cDNA.pos/cDNA.length',
       'CDS.pos/CDS.length', 'AA.pos/AA.length', 'Distance',
       'ERRORS/WARNINGS/INFO', 'ALL_TRANSCRIPT_START_END',
       'ALT_TRANSCRIPT_START_END', 'CDS_START_END', 'readDepth', 'altReadDepth', 'DP.REF', 'GENE',
       'GENOMIC_START_END', 'GT', 'HGMD_TRANSCRIPT', 'STRAND', 'TOTAL_EXONS','altFreq','type']


df_snpeff = df_snpeff[['gene','exon','chr', 'pos',  'ref', 'alt', 'impact', 'type','quality', 'altFreq', 'readDepth', 'altReadDepth', 'consequence','HGVSc', 'HGVSp',
'FILTER', 'Allele','Gene_Name', 'Feature_Type', 'Feature_ID', 'Transcript_BioType',
'cDNA.pos/cDNA.length','CDS.pos/CDS.length', 'AA.pos/AA.length', 'Distance', 'ERRORS/WARNINGS/INFO', 'ALL_TRANSCRIPT_START_END',
'ALT_TRANSCRIPT_START_END', 'CDS_START_END', 'DP.REF', 'GENE','GENOMIC_START_END', 'GT', 'HGMD_TRANSCRIPT', 'STRAND', 'TOTAL_EXONS','ID']]

df_snpeff["quality"] = [round(float(x),2)for x in df_snpeff.quality]
df_snpeff["HGVSc"] = [str(x)+':'+str(y) for x,y in zip(df_snpeff.Feature_ID,df_snpeff.HGVSc)]
df_snpeff["HGVSp"] = [str(x)+':'+str(y) for x,y in zip(df_snpeff.Feature_ID,df_snpeff.HGVSp)]

##### change ucsc zero based coordinate to one based coordinate by +1
df_snpeff["STRAND"] = ["REVERSE" if x=="-" else "FORWARD" for x in df_snpeff.STRAND]
df_snpeff["ALT_TRANSCRIPT_START"] = [int(x.split(':')[0])+1 for x in df_snpeff.ALT_TRANSCRIPT_START_END]
df_snpeff["ALT_TRANSCRIPT_END"] = [int(x.split(':')[1]) for x in df_snpeff.ALT_TRANSCRIPT_START_END]
df_snpeff["ALT_VARIANT_POSITION"] = updateAltTranscriptPosition(df_snpeff)
df_snpeff = updateNegativeStrandStartEndPosition(df_snpeff)


df_snpeff[['gene','exon','chr', 'pos',  'ref', 'alt', 'impact', 'type','quality','altFreq', 'readDepth', 'altReadDepth', 'consequence',
'HGVSc', 'HGVSp','STRAND','ALT_TRANSCRIPT_START','ALT_TRANSCRIPT_END','ALT_VARIANT_POSITION']].to_csv(OUT_FILE,index=False,sep='\t',header=None)

###save a file with all info for later use
df_snpeff.drop([ 'ERRORS/WARNINGS/INFO','ALL_TRANSCRIPT_START_END','Distance','Gene_Name'],axis=1,inplace=True)
df_snpeff.Allele = [x.replace("ANN=","") for x in df_snpeff.Allele]
df_snpeff.to_csv(OUT_FILE_ALL,index=False)
