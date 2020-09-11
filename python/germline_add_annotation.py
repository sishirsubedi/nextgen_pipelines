import sys
import ast
import pandas as pd
import mysql.connector
from mysql.connector import Error
from mysql.connector import errorcode
import subprocess

def bashCommunicator(command):
    process = subprocess.Popen([command],shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        print("Process failed %s %d %s %s" % (command,process.returncode, stdout, stderr))
    else:
        return [x for x in stdout.split("\n")]

def getConnection(db_host,db,user,passwd):
        try:
            return mysql.connector.connect(host=db_host,database=db, user=user, password=passwd)
        except mysql.connector.Error as error :
            print("Failed to connect database:{}".format(error))
            connection.rollback()

def getGnomadAF(connection,chr,pos,ref,alt):
    try:
        cursor = connection.cursor(prepared=True)
        query = "select AF from db_gnomad_r211_lf where chr = %s and pos = %s and ref = %s and alt = %s "
        cursor.execute(query,(chr,pos,ref,alt))
        af = []
        for row in cursor:
            af.append([el.decode('utf-8') if type(el) is bytearray else el for el in row])
        return af

    except mysql.connector.Error as error :
        print("Failed to get information from database:{}".format(error))
        connection.rollback()

def add_protein_prot2hg(connection, gene, feature,position):
    try:
        cursor = connection.cursor(prepared=True)
        query = "select protein_id,type,feature,note,prot_start,prot_end from db_prot2hg where gene = %s and gene_id = %s and  prot_start <= %s and prot_end >= %s limit 1 "
        cursor.execute(query,(gene,feature,position,position))
        domain_info = []
        for row in cursor:
            domain_info.append([el.decode('utf-8') if type(el) is bytearray else el for el in row])
        return domain_info

    except mysql.connector.Error as error :
        print("Failed to get information from database:{}".format(error))
        connection.rollback()

def add_protein_nextprot(connection, pos, feature):
    try:
        cursor = connection.cursor(prepared=True)
        query = "select domain from db_nextprot where hg19_start <= %s and hg19_end >= %s and transcript_id=%s "
        cursor.execute(query,(pos,pos,feature))
        domain = []
        for row in cursor:
            domain.append([el.decode('utf-8') if type(el) is bytearray else el for el in row])
        return domain

    except mysql.connector.Error as error :
        print("Failed to get information from database:{}".format(error))
        connection.rollback()

def add_protein_pfam(connection, pos,transcript):
    try:
        cursor = connection.cursor(prepared=True)
        query = "select hg19_knownGene_proteinID,hg19_pfamDesc_description, hg19_scopDesc_description from db_ucsc_pfam where hg19_knownGene_cdsStart <= %s and hg19_knownGene_cdsEnd >= %s and hg19_kgXref_refseq = %s limit 1"
        cursor.execute(query,(pos,pos,transcript))
        domain_info = []
        for row in cursor:
            domain_info.append([el.decode('utf-8') if type(el) is bytearray else el for el in row])
        return domain_info

    except mysql.connector.Error as error :
        print("Failed to get information from database:{}".format(error))
        connection.rollback()

def add_protein_uniprot(row, aa_change, gene_name):
    try:
        cursor = connection.cursor(prepared=True)
        query = "select AC,Source_DB_ID  from db_uniprot where Variant_AA_Change = %s and Gene_Name = %s limit 1"
        cursor.execute(query,(aa_change, gene_name))
        domain_info = []
        for row in cursor:
            domain_info.append([el.decode('utf-8') if type(el) is bytearray else el for el in row])
        return domain_info

    except mysql.connector.Error as error :
        print("Failed to get information from database:{}".format(error))
        connection.rollback()

def add_prediction_revel(connection, pos,ref,alt):
    try:
        cursor = connection.cursor(prepared=True)
        query = "select REVEL from db_revel where hg19_pos = %s and ref = %s and alt = %s limit 1"
        cursor.execute(query,(pos,ref,alt))
        score = []
        for row in cursor:
            score.append([el.decode('utf-8') if type(el) is bytearray else el for el in row])
        return score

    except mysql.connector.Error as error :
        print("Failed to get information from database:{}".format(error))
        connection.rollback()

def group_revel(row):
    #####grouping  revel scores based on published paper
    ## low - 0.0-0.25
    ## mid - 0.25-0.5
    ## high - 0.5-1.0
    ## https://sites.google.com/site/revelgenomics/
    ###########################################################
    if row.revel !="":
        score = float(row.revel)
        if score < 0.25:
            return "low ("+str(score)+")"
        elif score > 0.25 and score <0.5:
            return "mid ("+str(score)+")"
        elif score > 0.5:
            return "high ("+str(score)+")"
    else:
        return row.revel

def group_cadd(row):
    #####grouping  revel scores based on published paper
    ## low - 0-10
    ## mid - 10-20
    ## high - >20
    ## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6323892/
    ###########################################################
    if row.cadd_phred !="":
        score = float(row.cadd_phred)
        if score < 10:
            return "low ("+str(score)+")"
        elif score > 10 and score <20:
            return "mid ("+str(score)+")"
        elif score > 20:
            return "high ("+str(score)+")"
    else:
        return row.cadd_phred

def fixHGVSprotein_id(df):
    nm_ids = [x.split(":")[0] for x in df["HGVSp_x"]]
    nm_variant = [x.split(":")[1] for x in df["HGVSp_x"]]
    np_ids = []
    for np_id in df["HGVSp_y"]:
        try:
            np_ids.append(np_id.split(":")[0])
        except:
            np_ids.append("")

    adjusted_hgvsp = []
    for indx,row in enumerate(np_ids):
        if row != "": ## if np_id from vep is not empty them use id from vep and variant from snpeff
            adjusted_hgvsp.append(row+":"+str(nm_variant[indx]))
        else:
            if str(nm_variant[indx]) != "":
                adjusted_hgvsp.append(str(nm_ids[indx])+":"+str(nm_variant[indx]))
            else: ### if variant from snpeff is empty then just use vep output
                adjusted_hgvsp.append(df.ix[indx,"HGVSp_y"])

    return adjusted_hgvsp


################### INPUT
SAMPLE=sys.argv[1]
OUT_DIR=sys.argv[2]
environment = sys.argv[3]
db_host = sys.argv[4]
db = sys.argv[5]
user = sys.argv[6]
passwd = sys.argv[7]
################################################################################
SNPEFF_FILE = OUT_DIR+"/"+SAMPLE+".mutect.snpeff.parse.vcf"
VEP_FILE = OUT_DIR+"/"+SAMPLE+".mutect.vep.parse.vcf"

CADD_SNP = "/home/doc/ref/db_cadd/gnomad.genomes.r2.1.1.snv.tsv.gz"
CADD_INDEL = "/home/doc/ref/db_cadd/gnomad.genomes.r2.1.1.indel.tsv.gz"

############################################################################


#### read snpeff output
df_snpeff = pd.read_csv(SNPEFF_FILE,sep='\t',header=None)
df_snpeff.columns = ['gene','exon','CHROM', 'POS',  'REF', 'ALT', 'impact', 'type','quality',
'altFreq', 'readDepth', 'altReadDepth', 'consequence',
'HGVSc', 'HGVSp','STRAND','ALT_TRANSCRIPT_START','ALT_TRANSCRIPT_END','ALT_VARIANT_POSITION']
df_snpeff["Feature"] = [x.split(':')[0] for x in df_snpeff["HGVSp"]]


####

connection = getConnection(db_host,db,user,passwd)

global_af = []
for indx,row in df_snpeff.iterrows():
     global_af.append(getGnomadAF(connection,row.CHROM,row.POS,row.REF,row.ALT))

df_snpeff["Gnomad_GAF"] = [x[0][0] if len(x)>=1 else 0.0 for x in global_af]
# df_snpeff["Gnomad_GAF"] = df_snpeff["Gnomad_GAF"].astype(float)
# df_snpeff = df_snpeff[df_snpeff["Gnomad_GAF"]<1.0]


prot2hg = []
for indx,row in df_snpeff.iterrows():
    try:
        protein_position = row.HGVSp.split(":")[1].split(".")[1]
        protein_p =""
        for v in protein_position:
            if v.isdigit():
                protein_p +=v
        prot2hg.append(add_protein_prot2hg(connection,row.gene,row.Feature.split('.')[0],int(protein_p)))
    except:
        prot2hg.append([[""]])
entries = [x[0] if len(x)>=1 else "" for x in prot2hg]
df_snpeff["protein_id"] = [x[0] if len(x)>1 else "" for x in entries]
df_snpeff["protein_type"] = [x[1] if len(x)>1  else ""  for x in entries]
df_snpeff["protein_domain"] = [x[2] if len(x)>1  else "" for x in entries]
df_snpeff["protein_info"] = [x[3] if len(x)>1  else "" for x in entries]
df_snpeff["protein_start"] = [x[4] if len(x)>1  else "" for x in entries]
df_snpeff["protein_end"] = [x[5] if len(x)>1 else "" for x in entries]



nextprot = []
for indx,row in df_snpeff.iterrows():
     nextprot.append(add_protein_nextprot(connection,row.POS,row.Feature.replace(".","_")))
df_snpeff["NextProt_domain"] = [x[0][0] if len(x)>=1 else "" for x in nextprot]

pfam = []
for indx,row in df_snpeff.iterrows():
     pfam.append(add_protein_pfam(connection,row.POS,row.Feature.split('.')[0]))
entries = [x[0] for x in pfam]
df_snpeff["UniProt_ID"] = [x[0] for x in entries]
df_snpeff["Pfam"] = [x[1].split(",")[0].replace('"','').replace("n/a","").replace("\r","") for x in entries]
df_snpeff["Scoop"] = [x[2].split(",")[0].replace('"','').replace("n/a","").replace("\r","") for x in entries]

uniprot = []
for indx,row in df_snpeff.iterrows():
     uniprot.append(add_protein_uniprot(connection,row.HGVSp.split(':')[1],row.gene))
df_snpeff["UniProt_Variant"] = [x[0][0] if len(x)>=1 else "" for x in uniprot]
df_snpeff["ExPasy"] = [x[0][1] if len(x)>=1 else "" for x in uniprot]


revel = []
for indx,row in df_snpeff.iterrows():
     revel.append(add_prediction_revel(connection,row.POS,row.REF,row.ALT))
df_snpeff["revel"] = [x[0][0] if len(x)>=1 else "" for x in revel]
df_snpeff["revel"] = df_snpeff.apply(group_revel,axis=1)



cadd = []
for indx,row in df_snpeff.iterrows():
    cmd = ""
    if len(row.REF)==1 and len(row.REF)==1:
        cmd = "/opt/tabix-0.2.6/tabix %s %s:%s-%s" %(CADD_SNP,row.CHROM.replace("chr",""), row.POS, row.POS)
    else:
        cmd = "/opt/tabix-0.2.6/tabix %s %s:%s-%s" %(CADD_INDEL,row.CHROM.replace("chr",""), row.POS, row.POS)

    cadd_output = bashCommunicator(cmd)
    if cadd_output != None:
        if len(cadd_output)>1:
            cadd_output = cadd_output[0].split('\t')
            cadd.append([cadd_output[4],cadd_output[5]])
        else:
            cadd.append([""])
    else:
        cadd.append([""])

df_snpeff["cadd_phred"] = [x[1] if len(x)==2 else "" for x in cadd]
df_snpeff["cadd_phred"] = df_snpeff.apply(group_cadd,axis=1)


#### read vep output for sift and polyphen
df_vep = pd.read_csv(VEP_FILE,sep='\t')
df_snpeff_vep = pd.merge(df_snpeff,df_vep,on=["CHROM","POS","REF","ALT","Feature"])

### get HGSVp from vep
df_snpeff_vep["HGVSp_x"] = fixHGVSprotein_id(df_snpeff_vep)
df_snpeff_vep.drop(["HGVSc_y","HGVSp_y"],axis=1,inplace=True)
df_snpeff_vep.rename(columns={'HGVSc_x': 'HGVSc', 'HGVSp_x': 'HGVSp'},inplace=True)


###write output
df_snpeff_vep.to_csv(OUT_DIR+"/"+SAMPLE+".mutect.annotated.vcf",sep="\t",index=False)
