
### get cardiac genes list from medical director
 1_cardiac_genes.csv
## get unique genes
sort < 1_cardiac_genes.csv | uniq > 1_cardiac_genes_unq.csv

#### get ucsc genes with exons for these unique cardiac genes
#############################################################
import pandas as pd

def get_coding_coordinates(cds,cde,es,ee):
    exon_start = es.split(',')
    exon_end = ee.split(',')
    exon_pairs = [x for x in zip(exon_start,exon_end)]
    prime5 = []
    prime3 = []
    for indx in range(len(exon_pairs)-1):
        c_pair = exon_pairs[indx]
        if cds>=int(c_pair[0]) and cds<int(c_pair[1]):
            prime5.append(indx)
            prime5.append((cds,c_pair[1]))
            for indx2 in range(indx,len(exon_pairs)-1):
                c_pair = exon_pairs[indx2]
                if cde>int(c_pair[0]) and cde<=int(c_pair[1]):
                    prime3.append(indx2)
                    prime3.append((c_pair[0],cde))
                    break
            break

    coding_pairs = []

    if prime5[0] == prime3[0]:
        coding_pairs.append((cds,cde))
    else:
        coding_pairs.append(prime5[1])
        running_index = prime5[0] + 1
        while prime3[0] != running_index :
            coding_pairs.append(exon_pairs[running_index])
            running_index += 1
        coding_pairs.append(prime3[1])
    return coding_pairs


#### use HGMD crawler program to get  refseq from 1_cardiac_genes_unq.csv
### separate genes that needs manaual curation:
## -1_cardiac_genes_unq_withHMGD_RefSeq.csv
## - 1_cardiac_genes_unq_withHMGD_RefSeq_manual.csv

cardiac = pd.read_csv("1_cardiac_genes_unq_withHMGD_RefSeq.csv")
HMGD_RefSeq = cardiac['HMGD_RefSeq'].unique()


df = pd.read_csv("db_ucsc_genes.vcf",sep='\t')
df['ucsc_refseq'] = [str(x).split('.')[0] for x in df['name']]
df = df[df['ucsc_refseq'].isin(HMGD_RefSeq)]


df['coding_start'] = ''
df['coding_end'] = ''
for indx,row in df.iterrows():
    if int(row['exonCount']) == 1:
        df.ix[indx,'coding_start'] = row['cdsStart']
        df.ix[indx,'coding_end'] = row['cdsEnd']

    elif int(row['cdsStart']) < int(row['cdsEnd']):
        try :
            coding_pairs = get_coding_coordinates(row['cdsStart'],row['cdsEnd'], row['exonStarts'], row['exonEnds'])
            df.ix[indx,'coding_start'] = ';'.join([str(elem[0]) for elem in coding_pairs])
            df.ix[indx,'coding_end'] = ';'.join([str(elem[1]) for elem in coding_pairs])
        except :
            df.ix[indx,'coding_start'] = "ERROR"
            df.ix[indx,'coding_end'] = "ERROR"


df.columns = ['#bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart',
       'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'Gene',
       'cdsStartStat', 'cdsEndStat', 'exonFrames', 'ucsc_refseq',
       'coding_start', 'coding_end']

join = pd.merge(cardiac,df,on=['Gene'],how='left',indicator=True)

join.to_csv("3_cardiac_ucsc_join.csv",index=False)


df = pd.read_csv("3_cardiac_ucsc_join.csv")
df.head(2)


temp = []
for indx,row in df.iterrows():
    cdstart = row['coding_start'].split(';')
    cdend = row['coding_end'].split(';')
    for indx2 in range(len(cdstart)):
        temp.append([row['Gene'],row['HMGD_RefSeq'],row['chrom'],row['strand'],row['exonCount'],row['txStart'],row['txEnd'],row['cdsStart'],row['cdsEnd'],int(cdstart[indx2]),int(cdend[indx2])])
df_filt = pd.DataFrame(temp)
df_filt.columns = ['Gene','HMGD_RefSeq','chr','strand','exonCount','txStart','txEnd','cdsStart','cdsEnd','coding_start','coding_end']
df_filt = df_filt[['chr','coding_start','coding_end','Gene','HMGD_RefSeq','strand','exonCount','txStart','txEnd','cdsStart','cdsEnd']]
df_filt.sort_values(by=['chr','coding_start'],inplace=True)
df_filt.to_csv("4_ucsc_hgmd_cardiac_exons_rowwise.txt",sep='\t',index=False,header=None) #header none to bedtools intersect

### add MT genes, CAVIN4, PPA2, and TCRL gene info from curated HGMD 1_cardiac_genes_unq_withHMGD_RefSeq_manual_input.txt
#### curated output file is 4_ucsc_hgmd_cardiac_exons_rowwise.txt
####### intersect with alignment file to get coverage
use ${SAMPLE}.depth.bed
awk '{ if ( ( $1 !~ /\_/  )  ) {print $1 "\t" $2 "\t" $2+1 "\t" $3} }' NA12878-N.depth.bed > NA12878-N.depth.filter2.bed
/opt/bedtools2/bin/bedtools intersect -a  NA12878-N/NA12878-N.depth.filter2.bed -b 4_ucsc_hgmd_cardiac_exons_rowwise.txt -wa -wb > NA12878-N/NA12878-N.cardiac_intersect.bed


#### use above intersect file to obtain exonwise coverage result 

file_in = 'NA12878-N/NA12878-N.cardiac_intersect.bed'
file_out = 'NA12878-N/NA12878-N.cardiac_result.txt'


file_in = 'Exome46/Exome46-N.cardiac_intersect.bed'
file_out = 'Exome46/Exome46-N.cardiac_result.txt'

import pandas as pd
df = pd.read_csv(file_in,sep='\t',header=None)
df.sort_values(by=[0, 1],inplace=True)

def searchGap(df):
    table = []
    current = 0
    exon = 1
    tag = "temp"
    for indx,row in df.iterrows():

        if df.iloc[indx,3] >= 10 :
            tag = "high"
        else:
            tag = "low"

        start = df.iloc[indx,1]
        end = df.iloc[indx,2]

        if indx == 0:
            table.append([df.iloc[indx,0],start,end,tag,df.iloc[indx,3],exon,df.iloc[indx,5],df.iloc[indx,6],df.iloc[indx,7],
            df.iloc[indx,8],df.iloc[indx,9],df.iloc[indx,10],df.iloc[indx,11],df.iloc[indx,12],df.iloc[indx,13],df.iloc[indx,14]])
        else:
            if indx+1 < df.shape[0]:

                new_tag = "new_temp"
                if df.iloc[indx+1,3] >= 10 :
                    new_tag = "high"
                else:
                    new_tag = "low"

                if end == df.iloc[indx+1,1] and tag == new_tag: # next nt is adjacent and tag is also same then just change coordinate
                    table[current][2]=df.iloc[indx+1,2]
                elif end == df.iloc[indx+1,1] and tag != new_tag: # next nt is adjacent but tag is NOT same then change coordinate and tag BUT EXON is same
                    table.append([df.iloc[indx+1,0],df.iloc[indx+1,1],df.iloc[indx+1,2],
                    new_tag,df.iloc[indx+1,3],exon,df.iloc[indx+1,5],df.iloc[indx+1,6],df.iloc[indx+1,7],
                    df.iloc[indx+1,8],df.iloc[indx+1,9],df.iloc[indx+1,10],df.iloc[indx+1,11],df.iloc[indx+1,12],
                    df.iloc[indx+1,13],df.iloc[indx+1,14]])
                    current += 1
                elif end != df.iloc[indx+1,1]: # next nt is NOT adjacent then EXON is NOT same, tag is new tag
                    exon += 1
                    table.append([df.iloc[indx+1,0],df.iloc[indx+1,1],df.iloc[indx+1,2],
                    new_tag,df.iloc[indx+1,3],exon,df.iloc[indx+1,5],df.iloc[indx+1,6],df.iloc[indx+1,7],
                    df.iloc[indx+1,8],df.iloc[indx+1,9],df.iloc[indx+1,10],df.iloc[indx+1,11],df.iloc[indx+1,12],
                    df.iloc[indx+1,13],df.iloc[indx+1,14]])
                    current += 1

    return table



result = []
for gene in df[7].unique():
    df_genewise = df[df[7]==gene]
    df_genewise.reset_index(drop=True,inplace=True)
    res = searchGap(df_genewise)
    for row in res:
        result.append(row)

df_result = pd.DataFrame(result)
df_result.columns = ['chr','coding_start','coding_end','cov_status','coverage','exon_number','exon_start','exon_end','Gene','HMGD_RefSeq','strand','exonCount','txStart','txEnd','cdsStart','cdsEnd']
df_result.to_csv(file_out,sep='\t',index=False)
