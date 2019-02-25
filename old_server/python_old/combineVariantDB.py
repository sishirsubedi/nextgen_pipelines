import pandas as pd
import cyvcf2
import numpy as np
import optparse
from pymongo import MongoClient

client = MongoClient()
db = client.ngs_test
cosmic= db.cosmic
clinvar=db.clinvar
g1000=db.genome1000
oncokb=db.oncokb

infile='/home/environments/ngs_test/protonAnalysis/Auto_sn247770016_sn247770016-68-092717KPNeuro_145_203/IonXpress_018/variantCaller_out.586/TSVC_variants.split.vep.parse.newVarView.filter.txt'

def combineDB(infile, outfile):
    sample = pd.read_csv(infile,sep='\t',header=None)

    result=[]
    for indx,row in sample.iterrows():
        # print(row)
        chr=row[2]
        pos=row[3]
        ref=row[4]
        alt=row[5]

        cosmic_cursor = list(cosmic.find({'chr':int(chr[3:]),'pos':pos, 'ref':ref,'alt':alt},projection={'_id': 0,'cosmic-id':1}))
        cosmic_ids=[]
        for v in  cosmic_cursor:cosmic_ids.append(v['cosmic-id'])

        clinvar_cursor = list(clinvar.find({'chr':int(chr[3:]),'pos':pos, 'ref':ref,'alt':alt},projection={'_id': 0,'clinvar-id':1}))
        clinvar_ids=[]
        for v in  clinvar_cursor:clinvar_ids.append(v['clinvar-id'])

        g1000_cursor = list(g1000.find({'chr':chr,'pos':pos, 'ref':ref,'alt':alt},projection={'_id': 0,'g1000-id':1}))
        g1000_ids=[]
        for v in  g1000_cursor:g1000_ids.append(v['g1000-id'])


        result.append([chr,pos,ref,alt,cosmic_ids,clinvar_ids,g1000_ids])
    df_result=pd.DataFrame(result)
    df_result.columns=['chr','pos','ref','alt','cosmic_id','clinvar_id','g1000_id']


try:
    parser = optparse.OptionParser()
    parser.add_option('-i', '--infile', help = 'provide input file')
    parser.add_option('-o', '--outfile', help = 'provide output file')
    options,args = parser.parse_args()
    infile = options.infile
    outfile = options.outfile
    combineDB(infile, outfile)

except TypeError:
	print ("python combineVariantDB.py -help for help")
