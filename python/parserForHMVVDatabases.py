import pandas as pd
import cyvcf2
import numpy as np
import optparse

aaDic= {
"A" :	"Ala",
"C"	:	"Cys",
"D"	:	"Asp",
"E"	:	"Glu",
"F"	:	"Phe",
"G"	:	"Gly",
"H"	:	"His",
"I"	:	"Ile",
"K"	:	"Lys",
"L"	:	"Leu",
"M"	:	"Met",
"N"	:	"Asn",
"O"	:	"Pyl",
"P"	:	"Pro",
"Q"	:	"Gln",
"R"	:	"Arg",
"S"	:	"Ser",
"T"	:	"Thr",
"U"	:	"Sec",
"V"	:	"Val",
"W"	:	"Trp",
"Y"	:	"Tyr",
"*"	:	"Ter",
"X"	:	"Unknown"
}

def countLetter(s):
  d=l=0
  for c in s:
      if c.isdigit():
          d=d+1
      elif c.isalpha():
          l=l+1
      else:
          pass
  return l

def fileParser(infile, outfile, source):

    print ( "Processing file - " , infile )
    print ( "source is  - " , source )

    if source == 'cosmic':
        df = pd.read_csv(infile)
        filter=[]
        for indx,row in df.iterrows():
            info_line=row['INFO'].split(';')
            gene = info_line[0].split('=')[1]
            strand = info_line[1].split('=')[1]
            if info_line[2] == 'SNP':
                cds = info_line[3].split('=')[1]
                aa = info_line[4].split('=')[1]
                count = info_line[5].split('=')[1]
            else:
                cds = info_line[2].split('=')[1]
                aa = info_line[3].split('=')[1]
                count = info_line[4].split('=')[1]
            filter.append([row['ID'],row['CHROM'],row['POS'],row['REF'], row['ALT'],gene,strand,cds,aa,count])
        df_filter=pd.DataFrame(filter)
        df_filter.columns=['cosmic-id','chr','pos','ref','alt','gene','strand','cds','aa','count']
        df_filter.to_csv(outfile,index=False)

    elif source == 'clinvar':
        origin_code = {"0":"unknown","1":"germline", "2":"somatic",
        "4":"inherited", "8":"paternal", "16":"maternal","32":"de-novo",
        "64" :"biparental", "128" :"uniparental", "256":"not-tested",
         "512" : "tested-inconclusive", "1073741824" :"other"}

        df = pd.read_csv(infile,sep='\t',dtype={"CHROM": str, "POS": int, "ID":str, "REF":str,"ALT":str, "INFO":str})
        filter=[]
        for indx,row in df.iterrows():
            info_line=row['INFO'].split(';')
            CLNDN=''
            CLNSIG=''
            MC=''
            ORIGIN=''
            for pairs in info_line:
                code=pairs.split('=')[0]
                if code == 'CLNDN':
                    CLNDN=pairs.split('=')[1].split('|')[0]
                    if len(CLNDN.split(','))>1:
                        CLNDN= CLNDN.split(',')[0]
                elif code == 'CLNSIG':
                    CLNSIG=pairs.split('=')[1].split(',')[0]
                elif code == 'MC':
                    MC=pairs.split('=')[1].split('|')[1].split(',')[0]
                elif code == 'ORIGIN':
                    try:
                        ORIGIN=origin_code[pairs.split('=')[1]]
                    except:
                        ORIGIN='other'
            filter.append([row['ID'],'chr'+str(row['CHROM']),row['POS'],row['REF'], row['ALT'],CLNDN,CLNSIG,MC,ORIGIN])

        df_filter=pd.DataFrame(filter)
        df_filter.columns=['clinvar-id','chr','pos','ref','alt','CLNDN','CLNSIG','MC','ORIGIN']
        df_filter.to_csv(outfile,sep='\t',index=False)

    elif source == 'genome1000':
        df = pd.read_csv(infile,sep='\t')
        df.columns=["chr","pos","g1000-id","ref","alt","altCount","totalCount","altGlobalFreq","americanFreq","asianFreq","afrFreq", "eurFreq"]
        df = df[["g1000-id","chr","pos","ref","alt","altCount","totalCount","altGlobalFreq","americanFreq","asianFreq","afrFreq", "eurFreq"]]
        df.to_csv(outfile,index=False)

    elif source == 'oncokb':

        df = pd.read_csv(infile,sep='\t',encoding='latin-1')
        df = df.iloc[:,[0,1,3,5,6,7]]
        lf_pc=[]
        for indx,row in df.iterrows():
          pc=row['Protein Change']
          if countLetter(pc) ==2 and len(pc)>2:
              lf_pc.append(aaDic[pc[0]]+pc[1:len(pc)-1]+aaDic[pc[len(pc)-1]])
          else:
              lf_pc.append(pc)
        df['Protein Change LF'] = lf_pc
        df.to_csv(outfile,index=False)

    elif source == 'civic':
        df = pd.read_csv(infile,sep='\t')
        df = df[['gene', 'variant', 'variant_origin' ,'variant_civic_url']]

        lf_pc=[]
        for indx,row in df.iterrows():
          pc=row['variant']
          if countLetter(pc) ==2 and len(pc)>2 and '/' not in pc and pc[1].isdigit():
              lf_pc.append(aaDic[pc[0]]+pc[1:len(pc)-1]+aaDic[pc[len(pc)-1]])
              # print(pc, aaDic[pc[0]]+pc[1:len(pc)-1]+aaDic[pc[len(pc)-1]])
          else:
              lf_pc.append(pc)
        df['variant_LF'] = lf_pc
        df.to_csv(outfile,index=False)

    elif source == 'gnomad_lf':

        codes=['AF','AF_afr', 'AF_amr', 'AF_asj', 'AF_eas','AF_fin','AF_nfe','AF_sas','AF_oth','AF_male','AF_female']

        vcf = cyvcf2.VCF(infile)

        with open(outfile,"w") as w:
            header = ['gnomad-id','chr','pos','ref','alt','AF','AF_afr', 'AF_amr', 'AF_asj', 'AF_eas','AF_fin','AF_nfe','AF_sas','AF_oth','AF_male','AF_female']
            bunch=100000
            count=1
            filter=[]
            w.writelines('\t'.join(str(j) for j in header) + '\n')
            for row in vcf:
                current =[]
                current.append(row.ID)
                current.append('chr'+str(row.CHROM))
                current.append(row.POS)
                current.append(row.REF)
                current.append(str(row.ALT).replace(']','').replace('[','').replace("'",''))

                for code in codes:
                    if str(row.INFO.get(code)) == 'None':
                        current.append(0.0)
                    else:
                        current.append(round(row.INFO.get(code)*100,5))

                filter.append(current)

                if len(filter) == bunch:
                    print("processed " , count , " hundred thousand rows.")
                    w.writelines('\t'.join(str(j) for j in i) + '\n' for i in filter)
                    filter=[]
                    count += 1
            w.writelines('\t'.join(str(j) for j in i) + '\n' for i in filter)
            w.close()


    elif source == 'gnomad':

        vcf = cyvcf2.VCF(infile)

        with open(outfile,"w") as w:
            header = ['gnomad-id','chr','pos','ref','alt','AF']
            bunch=100000
            count=1
            filter=[]
            w.writelines('\t'.join(str(j) for j in header) + '\n')
            for row in vcf:
                current =[]
                current.append(row.ID)
                current.append('chr'+str(row.CHROM))
                current.append(row.POS)
                current.append(row.REF)
                current.append(str(row.ALT).replace(']','').replace('[','').replace("'",''))
                if str(row.INFO.get('AF')) == 'None':
                    current.append(0.0)
                else:
                    current.append(round(row.INFO.get('AF')*100,5))
                filter.append(current)

                if len(filter) == bunch:
                    print("processed " , count , " hundred thousand rows.")
                    w.writelines('\t'.join(str(j) for j in i) + '\n' for i in filter)
                    filter=[]
                    count += 1
            w.writelines('\t'.join(str(j) for j in i) + '\n' for i in filter)
            w.close()


    elif source =='pmkb':
        df=pd.read_csv(infile,dtype={"Gene": str, "Tumor Type(s)": str, "Tissue Type(s)":str, "Variant(s)":str})
        df.columns=['gene','tumor_type','tissue_type','variant']
        filter=[]
        for indx,row in df.iterrows():

            if str(row['variant']) == 'nan':
                continue

            variants = row['variant'].split(',')

            for v in variants:

                variant_sep= v.replace(row['gene'],'').replace(' ','')

                variant_lf=""
                if countLetter(variant_sep) ==2 and len(variant_sep)>2 and '/' not in variant_sep and variant_sep[1].isdigit():
                    variant_lf=aaDic[variant_sep[0]]+variant_sep[1:len(variant_sep)-1]+aaDic[variant_sep[len(variant_sep)-1]]
                else:
                    variant_lf=variant_sep


                tumors=''
                if str(row['tumor_type']) != 'nan':
                    tumor_type= row['tumor_type'].split(',')
                    if len(tumor_type)>2:
                        tumors = ','.join(tumor_type[0:2])
                    else:
                        tumors=tumor_type[0]

                tissues=''
                if str(row['tissue_type']) != 'nan':
                    tissue_type=row['tissue_type'].split(',')
                    if len(tissue_type)>1:
                        tissues = ','.join(tissue_type[0:2])
                    else:
                        tissues = tissue_type[0]

                filter.append([row['gene'],tumors,tissues,variant_sep,variant_lf])

                df_filter=pd.DataFrame(filter)
                df_filter.columns=['gene','tumor_type','tissue_type','variant','variant_lf']
                df_filter.to_csv(outfile,sep='\t',index=False)

    elif source == 'cardiac':

        import requests
        from bs4 import BeautifulSoup


        cardiac_genes =["ABCC9","ACTC1","ACTN2","ANK2","ANKRD1","CALR3","CASQ2","CAV3","CRYAB","CSRP3","CTF1","DES","DSC2","DSG2","DSP","DTNA", \
        "EMD","FHL1","FHL2","FXN","GLA","JPH2","JUP","KCNE1","KCNE2","KCNH2","KCNJ2","KCNQ1","KLF10","LAMA4","LAMP2","LDB3","LMNA","MYBPC3","MYH6",\
        "MYH7","MYL2","MYL3","MYLK2","MYOM1","MYOZ2","MYPN","NEXN","OBSCN","PDLIM3",\
        "PKP2","PLN","PRKAG2","RBM20","RYR2","SCN5A","SGCD","SRI","TAZ","TCAP","TGFB3","TMEM43","TNNC1","TNNI3","TNNT2","TPM1","TRIM54","TRIM55","TRIM63","TTN","TTR","VCL"]


        for gene in cardiac_genes:
            print("processing...."+gene)
            url = 'http://www.cardiodb.org/acgv/acgv_exac_variants.php?gene='+gene+'&vartype=all'
            page = requests.get(url)
            soup = BeautifulSoup(page.content, 'html.parser')
            tbl = soup.find("table",{"id":"gene_variants"})
            df = pd.read_html(str(tbl),flavor='html5lib')[0]
            df = df.iloc[1:,[1,2,3,4]]
            df.to_csv(gene+".csv",index=False)

        df_combine = pd.DataFrame()
        for gene in cardiac_genes:
            df = pd.read_csv(gene+".csv")
            df.columns=["genomic_coordinate","cds_variant","protein_variant","variant_type"]
            df["gene"]=gene
            df = df[["gene","genomic_coordinate","cds_variant","protein_variant","variant_type"]]

            if gene=="ABCC9":
                df_combine = df
            else:
                df_combine = df_combine.append(df)

        df_combine.to_csv(outfile,sep='\t',index=False)

try:
    parser = optparse.OptionParser()
    parser.add_option('-i', '--infile', help = 'provide input file')
    parser.add_option('-o', '--outfile', help = 'provide output file')
    parser.add_option('-s', '--source', help = 'provide file source')
    options,args = parser.parse_args()
    infile = options.infile
    outfile = options.outfile
    source = options.source
    fileParser(infile, outfile, source)

except TypeError:
	print ("python parserForHMVVDatabases.py -help for help")
