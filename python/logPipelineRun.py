import sys
import os
import glob
import optparse
import mysql.connector #/opt/python3/bin/pip3 install mysql-connector
from mysql.connector import Error
from mysql.connector import errorcode

class Pipeline:
    def __init__(self, elements):
        self.queueID = elements[0]
        self.runID = elements[1]
        self.sampleName = elements[2]
        self.assay = elements[3]
        self.instrument = elements[4]
        self.timeSubmitted = elements[5]


def logIntoDB(queueID,user,passwd,db):

    pipeline=getPipelineInfo(queueID,user, passwd,db)

    # for attr, value in pipeline.__dict__.items():
    #     print(attr, value)

    #### calculate the file size
    vcf_file=''
    if pipeline.instrument == 'miseq':
        vcf_file = glob.glob(os.path.join('/home', pipeline.instrument, '*_'+pipeline.runID+'_*', 'Data','Intensities', 'BaseCalls', 'Alignment' ,pipeline.sampleName+'_*.vcf' ))
    elif pipeline.instrument == 'proton':
        vcf_file = glob.glob(os.path.join('/home', pipeline.instrument, '*'+pipeline.runID, 'plugin_out','variantCaller_out*',pipeline.sampleName,'TSVC_variants.vcf'))
    elif pipeline.instrument == 'nextseq': ##check fastq here
        vcf_file = glob.glob(os.path.join('/home', pipeline.instrument, '*_'+pipeline.runID+'_*', 'out1',pipeline.sampleName+'*_R1_001.fastq.gz'))
    filesize = os.stat(vcf_file[0]).st_size / 1000 # store in KB

    ## calculate the total run timeUpdated
    runtime = getPipelineRunTime(queueID,user, passwd,db)

    updatePipelineLog(user,passwd,db, pipeline.queueID, filesize, runtime)


def getPipelineInfo(queueID,user,passwd,db):

    try:
        connection = mysql.connector.connect(host='localhost',database=db, user=user, password=passwd)
        cursor = connection.cursor(prepared=True)
        query = "select pipelineQueue.queueID, samples.runID, samples.sampleName, assays.assayName, instruments.instrumentName, pipelineQueue.timeSubmitted from pipelineQueue join samples on samples.sampleID=pipelineQueue.sampleID join assays on assays.assayID = samples.assayID join instruments on instruments.instrumentID = samples.instrumentID where pipelineQueue.queueID= %s"
        cursor.execute(query, (queueID,))
        pipeline = None;
        for row in cursor:
            pipeline = Pipeline([el.decode('utf-8') if type(el) is bytearray else el for el in row])
        return pipeline

    except mysql.connector.Error as error :
        print("Failed to update record to database:{}".format(error))
        connection.rollback()

    finally:
        if(connection.is_connected()):
            connection.close()

def getPipelineRunTime(queueID,user,passwd,db):

    try:
        connection = mysql.connector.connect(host='localhost',database=db, user=user, password=passwd)
        cursor = connection.cursor(prepared=True)
        query = "SELECT TIME_TO_SEC(TIMEDIFF(t2.timeUpdated,t1.timeUpdated))  as difference from pipelineStatus as t1 inner join pipelineStatus as t2 on  t1.queueID = t2.queueID where t1.queueID=%s and t1.plStatus='started' and t2.plStatus='pipelineCompleted'"
        cursor.execute(query, (queueID,))
        timediff=0
        for row in cursor:
            timediff=row[0]
        return timediff

    except mysql.connector.Error as error :
        print("Failed to update record to database:{}".format(error))
        connection.rollback()

    finally:
        if(connection.is_connected()):
            connection.close()

def updatePipelineLog(user,passwd,db, queueID, filesize, runtime):

    try:
        connection = mysql.connector.connect(host='localhost',database=db, user=user, password=passwd)
        cursor = connection.cursor(prepared=True)
        query = "INSERT INTO pipelineLogs (queueID,fileSizeKB,runTimeSecs) VALUES(%s,%s,%s) "
        cursor.execute(query, (queueID, filesize, runtime))
        connection.commit()

    except mysql.connector.Error as error :
        print("Failed to update record to database:{}".format(error))
        connection.rollback()

    finally:
        if(connection.is_connected()):
            connection.close()


try:
    parser = optparse.OptionParser()
    parser.add_option('-q', '--queueID', help = ' input queueID to log pipeline run')
    parser.add_option('-u', '--user', help = ' input user to connect to database')
    parser.add_option('-p', '--password', help = ' input user password to connect to database')
    parser.add_option('-d', '--database', help = ' input database')
    options,args = parser.parse_args()
    queueID = options.queueID
    user = options.user
    passwd = options.password
    db = options.database
    logIntoDB(queueID,user,passwd,db)

except TypeError:
	print ("python logPipelineRun.py -help for help")
