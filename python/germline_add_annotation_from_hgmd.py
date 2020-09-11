import sys
import ast
import pandas as pd
import mysql.connector
from mysql.connector import Error
from mysql.connector import errorcode
import json 
import datetime

def getConnection(db_host,db,user,passwd):
        try:
            return mysql.connector.connect(host=db_host,database=db, user=user, password=passwd)
        except mysql.connector.Error as error :
            print("Failed to connect database:{}".format(error))
            connection.rollback()

def get_hgmd_id(connection,chr,pos,ref,alt):
    try:
        cursor = connection.cursor(prepared=True)
        query = "select id from hgmd_hg19_vcf where chrom = %s and pos = %s and ref = %s and alt = %s limit 1"
        cursor.execute(query,(chr,pos,ref,alt))
        id = []
        for row in cursor:
            id.append([el.decode('utf-8') if type(el) is bytearray else el for el in row])
        return id

    except mysql.connector.Error as error :
        print("Failed to get information from database:{}".format(error))
        connection.rollback()

def get_hgmd_info(connection,id,table):
    try:
        cursor = connection.cursor(prepared=True)
        
        query=""
            
        if table == "mutation":
            query = "select \
            mutation.disease,mutation.gene,mutation.base,mutation.amino,mutation.codon,mutation.tag,\
            mutation.pmid,mutation.comments,mutation.new_date\
            from mutation where mutation.acc_num = %s "
        
        elif table == "splice":
            query = "select \
            splice.disease,splice.gene,splice.ivs,splice.type,splice.base,splice.location,splice.tag,\
            splice.pmid,splice.comments, splice.new_date\
            from splice where splice.acc_num = %s "
        
        elif table == "prom":
            query = "select \
            prom.disease,prom.gene,prom.base,prom.location,prom.locref,prom.tag,\
            prom.pmid,prom.comments, prom.new_date\
            from prom where prom.acc_num = %s "

        elif table == "deletion":
            query = "select \
            deletion.disease, deletion.gene,deletion.deletion,deletion.codon,deletion.tag,\
            deletion.pmid,deletion.comments, deletion.new_date\
            from deletion where deletion.acc_num = %s "

        elif table == "insertion":
            query = "select \
            insertion.disease,insertion.gene,insertion.insertion,insertion.codon,insertion.nucleotide,insertion.tag,\
            insertion.pmid,insertion.comments, insertion.new_date\
            from insertion where insertion.acc_num = %s "

        elif table == "indel":
            query = "select \
            indel.disease,indel.gene,indel.wildtype,indel.insertion,indel.codon,indel.tag,\
            indel.pmid,indel.comments, indel.new_date\
            from indel where indel.acc_num = %s "

        elif table == "grosdel":
            query = "select \
            grosdel.disease,grosdel.gene,grosdel.descr,grosdel.cdna,grosdel.tag,\
            grosdel.pmid,grosdel.comments, grosdel.new_date\
            from grosdel where grosdel.acc_num = %s "


        elif table == "grosins":
            query = "select \
            grosins.disease,grosins.gene,grosins.type,grosins.descr,grosins.cdna,grosins.tag,\
            grosins.pmid,grosins.comments, grosins.new_date\
            from grosins where grosins.acc_num = %s "


        elif table == "complex":
            query = "select \
            complex.disease,complex.gene,complex.descr,complex.tag,\
            complex.pmid,complex.comments, complex.new_date\
            from complex where complex.acc_num = %s "

        elif table == "amplet":
            query = "select \
            amplet.disease,amplet.gene,amplet.chrom,amplet.amplet,amplet.norcopy,amplet.patcopy,amplet.location,amplet.tag,\
            amplet.pmid,amplet.comments, amplet.new_date\
            from amplet where amplet.acc_num = %s "

        elif table == "func_anotat":
            query = "select \
            func_anotat.site_id,func_anotat.direction,func_anotat.wildtype,func_anotat.mutant,\
            func_anotat.pmid,func_anotat.comments\
            from func_anotat where func_anotat.acc_num = %s "

        elif table == "extrarefs":
            query = "select \
            extrarefs.pmid,extrarefs.comments,extrarefs.reftag \
            from extrarefs where extrarefs.acc_num = %s "


        cursor.execute(query,(id,))
        mutation = []
        for row in cursor:
            mutation.append([el.decode('utf-8') if type(el) is bytearray else el for el in row])
        return mutation

    except mysql.connector.Error as error :
        print("Failed to get information from database:{}".format(error))
        connection.rollback()

def get_hgmd_columns():
    return {
    "mutation":"mutation_disease,mutation_gene,mutation_base,mutation_amino,mutation_codon,mutation_tag,\
            mutation_pmid,mutation_comments,mutation_new_date",
    "splice":"splice_disease,splice_gene,splice_ivs,splice_type,splice_base,splice_location,splice_tag,\
            splice_pmid,splice_comments, splice_new_date",        
    "prom":"prom_disease,prom_gene,prom_base,prom_location,prom_locref,prom_tag,\
            prom_pmid,prom_comments, prom_new_date",
    "deletion":"deletion_disease, deletion_gene,deletion_deletion,deletion_codon,deletion_tag,\
            deletion_pmid,deletion_comments, deletion_new_date",
    "insertion":"insertion_disease,insertion_gene,insertion_insertion,insertion_codon,insertion_nucleotide,insertion_tag,\
            insertion_pmid,insertion_comments, insertion_new_date",
    "indel":"indel_disease,indel_gene,indel_wildtype,indel_insertion,indel_codon,indel_tag,\
            indel_pmid,indel_comments, indel_new_date",
    "grosdel":"grosdel_disease,grosdel_gene,grosdel_descr,grosdel_cdna,grosdel_tag,\
            grosdel_pmid,grosdel_comments, grosdel_new_date",
    "grosins":"grosins_disease,grosins_gene,grosins_type,grosins_descr,grosins_cdna,grosins_tag,\
            grosins_pmid,grosins_comments, grosins_new_date",
    "complex":"complex_disease,complex_gene,complex_descr,complex_tag,\
            complex_pmid,complex_comments, complex_new_date",
    "amplet":"amplet_disease,amplet_gene,amplet_chrom,amplet_amplet,amplet_norcopy,amplet_patcopy,amplet_location,amplet_tag,\
            amplet_pmid,amplet_comments, amplet_new_date",
    "func_anotat":"func_anotat_site_id,func_anotat_direction,func_anotat_wildtype,func_anotat_mutant,\
            func_anotat_pmid,func_anotat_comments",
    "extrarefs":"extrarefs_pmid,extrarefs_comments,extrarefs_reftag"
    }

def update_hgmd_info(hgmd_dict):
    hgmd_columns = get_hgmd_columns()
    for table in hgmd_columns.keys(): 
        hgmd_table_row = get_hgmd_info(hgmd_connection,hgmd_id,table)
        if hgmd_table_row !=[]:
            for indx,col in enumerate(hgmd_columns[table].replace(" ","").split(",")):
                hgmd_dict[col] = str(hgmd_table_row[0][indx]).replace(",","-").replace(":","-")

################### INPUT
SAMPLE=sys.argv[1]
OUT_DIR=sys.argv[2]
environment = sys.argv[3]
db_host = sys.argv[4]
db = sys.argv[5]
user = sys.argv[6]
passwd = sys.argv[7]
################################################################################
VCF_FILE = OUT_DIR+"/"+SAMPLE+".mutect.annotated.vcf"
############################################################################
df_vcf = pd.read_csv("VCF_FILE",sep='\t')


hgmd_connection = getConnection(db_host,"hgmd_pro",user,passwd)
hgmd_column = []
for indx,row in df_vcf.iterrows():
    hgmd_id = get_hgmd_id(hgmd_connection, row.CHROM.replace("chr",""),row.POS,row.REF,row.ALT)
    hgmd_dict = {}
    hgmd_dict["variant_number"]=str(indx+1)
    hgmd_dict["variant_id"]=row.CHROM+"_"+str(row.POS)+"_"+row.REF+"_"+row.ALT
    if hgmd_id != []:
        hgmd_id = hgmd_id[0][0]
        hgmd_dict["hgmd_information"] = "Present"
        hgmd_dict["hgmd_id"] = hgmd_id
        update_hgmd_info(hgmd_dict)
        hgmd_column.append(hgmd_dict)

f = open("test.txt", "w")
for variant in hgmd_column:
    f.write(str(variant)+'\n')
f.close()
