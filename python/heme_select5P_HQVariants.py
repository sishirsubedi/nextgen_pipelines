import sys
import pandas as pd
import optparse
import mysql.connector
from mysql.connector import Error
from mysql.connector import errorcode

def filterVariants(path,sampleName,mode):

    vc1_file = path+sampleName+".combine_VC2of3.vcf"
    df_vc1 = pd.read_csv(vc1_file,sep='\t',skiprows=1)
    main_columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTE', 'INFO', 'FORMAT','sample']
    df_vc1.columns=main_columns


    df_vc1['vf'] = [ x[1] for x in df_vc1['sample'].str.split(':')]
    df_vc1['vf'] = df_vc1['vf'].astype(float)


    if mode == "Low":
        df_vc1 = df_vc1[( df_vc1['vf']>=0.05) & (df_vc1['vf']<0.1) ]
    elif mode == "High":
        df_vc1 = df_vc1[( df_vc1['vf']>=0.1)]

    df_vc1.reset_index(drop=True,inplace=True)

    return df_vc1

def getConnection(db_host,db,user,passwd):
    try:
        return mysql.connector.connect(host=db_host,database=db, user=user, password=passwd)
    except mysql.connector.Error as error :
        print("Failed to connect database:{}".format(error))
        connection.rollback()

def checkDB(connection,tab,chr,pos,ref,alt):
    try:
        found = 0
        cursor = connection.cursor(buffered=True)
        query =  "SELECT EXISTS ( select * from "+ tab + " where chr = %s and pos = %s and ref = %s and alt = %s limit 1 )"
        cursor.execute(query, (chr,pos,ref,alt,))

        for r in cursor : found = r[0]

        return found

    except mysql.connector.Error as error :
        print("Failed :{}".format(error))
        connection.rollback()



def dbFilter5Pvariants(connection,df_all):
    db_cosmic = "db_cosmic_grch37v86"
    db_clinvar = "db_clinvar_42019"
    db_index = []
    for indx,row in df_all.iterrows():
        if checkDB(connection,db_cosmic,row['CHROM'], row['POS'], row['REF'], row['ALT']):
            db_index.append(indx)
        elif checkDB(connection,db_clinvar,row['CHROM'], row['POS'], row['REF'], row['ALT']):
            db_index.append(indx)
    return df_all.iloc[db_index,:]

def select5Pvariants(connection,path,sampleName,outFile):

    df_10Above = filterVariants(path,sampleName,"High")

    df_5to10 = filterVariants(path,sampleName,"Low")

    df_5to10_filter = dbFilter5Pvariants(connection,df_5to10)

    ###merge both
    df_filter = pd.concat([df_10Above, df_5to10_filter], ignore_index=True)

    df_filter.to_csv(path+outFile,sep='\t',header=False,index=False)

try:
    parser = optparse.OptionParser()
    parser.add_option('-u', '--user', help = ' input user to connect to database')
    parser.add_option('-p', '--password', help = ' input user password to connect to database')
    parser.add_option('-d', '--database', help = ' input database')
    parser.add_option('-b', '--dbhost', help = ' input database host')
    parser.add_option('-f', '--filepath', help = ' input filepath')
    parser.add_option('-s', '--sampleName', help = ' input sampleName')
    parser.add_option('-o', '--outfile', help = ' input outfile')


    options,args = parser.parse_args()
    user = options.user
    passwd = options.password
    db = options.database
    db_host = options.dbhost

    path = options.filepath
    sampleName = options.sampleName
    outFile = options.outfile

    connection = getConnection(db_host,db,user,passwd)
    select5Pvariants(connection,path,sampleName,outFile)

    if(connection.is_connected()):
        connection.close()

except TypeError:
	print ("python heme_select5P_HQVariants.py -help for help")
