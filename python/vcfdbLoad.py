import sys
import optparse
import re
from operator import itemgetter

cln_dic = {"0":"Uncertain_significance", "1":"not_provided" , "2":"benign",
"3":"LIkely_benign", "4":"Likely_Pathogenic", "5":"Pathogenic", "6":"drug_response",
"7":"histocompatibility", "255":"other"}

def brca_her_load(infile,outfile,inref,sampleID):
	ref_dic = vcf_to_dic(inref)
	infile = open(infile, 'r')
	outfile = open(outfile, 'w')
	outfile.write('#' + ' '.join(sys.argv) + '\n')
	outfile.write('@chr\tpos\tsnpID\tref\talt\tfrequency\ttype\tgenotype\tmulti\tfixed\tbase_corrected\tsample\n')
	sampleID = sampleID.strip()
	for line in infile:
		if line.startswith('#'):
			continue
		else:
			line = line.split()
			genotype = line[9].split(':')[0]
			if genotype == "./." or genotype == ".|." or genotype == "0/0" or genotype == "0|0":
				continue
			vcf_db = vcf_to_db(line, ref_dic)
			for record in vcf_db:
				#output record
				outfile.write('\t'.join(record) + '\t' + sampleID + '\n')
	infile.close()
	outfile.close()

def correct_base(ref, alt):
	ref_right = ref[1:]
	alt_right = alt[1:]
	base_corrected = "N"
	if len(ref_right) < len(alt_right):
		#insertion
		regex = ref_right + "$"
		m = re.search(regex, alt_right)
		if m:
			# in case of a simple insertion event, get the inserted bases
			insertion = alt_right[0:m.start()]
			ref = ref[0]
			alt = alt[0] + insertion
			base_corrected = "Y"
	elif len(ref_right) > len(alt_right):
		#deletion
		regex = alt_right + "$"
		m = re.search(regex, ref_right)
		if m:
			# in case of a simple deletion event, get the deleted bases
			deletion = ref_right[0:m.start()]
			ref = ref[0] + deletion
			alt = alt[0]
			base_corrected = "Y"
	else:
		#check if it is simple SNP case
		if ref_right == alt_right:
			ref = ref[0]
			alt = alt[0]
			base_corrected = "Y"
	return(ref, alt, base_corrected)
	
	
def vcf_to_db(oneVcf, ref_dic):
	multi = "N"  #record if it is a multiple dbSNP on one coordinate
 	vcf_db=[]  # output matrix, one uniq dbSNP per row
	if not oneVcf[0].startswith('chr'):
		oneVcf[0] = 'chr' + oneVcf[0]
	#grab info column and parse into dictionary
	info_dic = {}
	info = oneVcf[7].split(";")
	for field in info:
		field = field.split("=")
		try:
			info_dic[field[0]] = field[1]
		except IndexError:
			continue
	dbSNP = oneVcf[2]
	alt = oneVcf[4]
	ref = oneVcf[3]
	chr = oneVcf[0]
	pos = oneVcf[1]
	genotype = oneVcf[9].split(':')[0]
	if len(dbSNP.split(";")) > 1 or len(alt.split(",")) > 1:
		#multiple dbSNP or alt alleles per line
		multi = "Y"
		snpIDs = dbSNP.split(";")
		altAlleles = alt.split(",")
		altCounts = info_dic["FAO"].split(",")
		n = 0
		snpN = len(snpIDs) - 1
		altN = len(altAlleles) - 1
		altCountN = len(altCounts) - 1
		maxN = max(snpN, altN, altCountN)
		while n <= maxN:
			#split multiple dbSNP or alt Alleles
			refAllele = ref
			base_corrected = "N"
			fixed = "N" #record if the dbSNP id is reassigned from the clinvar file
			snpID = snpIDs[min(n,snpN )]
			altAllele = altAlleles[min(n,altN )]
			altCount = altCounts[min(n,altCountN )]
			frequency = float(altCount)/float(info_dic["FDP"])  * 100
			frequency = str(round(frequency, 2))
			#change ref and alt column for ins/del to match clinvar
			if (len(refAllele) > 1 and len(altAllele) > 1):
				(refAllele, altAllele, base_corrected) = correct_base(refAllele, altAllele)
			refDicKey = chr + "|" + pos + "|" + refAllele + "|" + altAllele
			#print refDicKey
			#check if the dbSNP ID assignment is the same as in the reference file
			if refDicKey in ref_dic:
				#print "coordinates found in ref"
				if not snpID in ref_dic[refDicKey]:
					#print "dbSNP assignment not correct"
					snpID = ref_dic[refDicKey][0]
					fixed = "Y"
			if len(altAllele) == 1 and len(refAllele) == 1:
				type = "SNP"
			elif len(altAllele) < len(refAllele):
				type = "del"
			elif len(altAllele) > len(refAllele):
				type = "ins"
			else:
				type = "other"
			vcf_db.append([chr,pos,snpID,refAllele,altAllele,frequency,type,genotype,multi,fixed,base_corrected])
			n += 1
	else:
		#single dbSNP and alt allele per line
		fixed = "N" #record if the dbSNP id is reassigned from the clinvar file
		base_corrected = "N"
		#change ref and alt column for ins/del to match clinvar
		if (len(ref) > 1 and len(alt) > 1):
			base_corrected = "Y"
			(ref, alt, base_corrected) = correct_base(ref, alt)
		refDicKey = chr + "|" + pos + "|" + ref + "|" + alt
		if refDicKey in ref_dic:
			if not dbSNP in ref_dic[refDicKey]:
				dbSNP = ref_dic[refDicKey][0]
				fixed = "Y"
		frequency = float(info_dic["FAO"])/float(info_dic["FDP"])  * 100
		frequency = str(round(frequency, 2))
		if len(alt) == 1 and len(ref) == 1:
			type = "SNP"
		elif len(alt) < len(ref):
			type = "del"
		elif len(alt) > len(ref):
			type = "ins"
		else:
			type = "other"
		vcf_db.append([chr,pos,dbSNP,ref,alt,frequency,type,genotype,multi,fixed,base_corrected])
	return vcf_db
	
def vcf_to_dic(vcf):
	infile = open(vcf, 'r')
	outdic = {}
	for line in infile:
		if line.startswith('#'):
			continue
		else:
			line = line.split()
			if not line[0].startswith('chr'):
				line[0] = 'chr' + line[0]
			alts = line[4].split(',')
			for alt in alts:
				key = line[0] + '|' + line[1] + '|' + line[3] + '|' + alt
				if not key in outdic:
					outdic[key] = [line[2]]
				else:
					outdic[key].append(line[2])
	infile.close()
	return outdic

def clinvar_load(infile, outfile):
	infile = open(infile, 'r')
	outfile = open(outfile, 'w')
	outfile.write('#' + ' '.join(sys.argv) + '\n')
	outfile.write('@chr\tpos\tdbSNP\tref\talt\tgene\tclnsig\tclnsig_db\tclnsig_disease\tclnallele\tcaf\tcommon\n')
	for line in infile:
		if line.startswith('#'):
			continue
		else:
			line = line.split()
			if not line[0].startswith('#'):
				line[0] = 'chr' + line[0]
			#parse info column into dic
			info_dic = {}
			info = line[7].split(";")
			for field in info:
				field = field.split("=")
				try:
					info_dic[field[0]] = field[1]
				except IndexError:
					continue
			#assign values to relevant parameters from info dic
			info_list = ["GENEINFO", "CLNSIG", "CLNDSDB", "CLNDBN", "CLNALLE","CAF","COMMON"]
			info_value = []
			for info in info_list:
				try:
					info_value.append(info_dic[info])
				except KeyError:
					info_value.append("NA")
			
			alts = line[4].split(',')
			n = 0
			caf = info_value[5]
			clnSig = info_value[1]
			clnAllele = info_value[4]
			clnDb = info_value[2]
			clnDisease = info_value[3]
			clnSigmaxN = len(clnSig.split(',')) - 1
			cafmaxN = len(caf.split(',')) - 1
			clnDbMaxN = len(clnDb.split(',')) - 1
			clnDiseaseMaxN = len(clnDisease.split(',')) - 1
			clnAlleleMaxN = len(clnAllele.split(',')) - 1
			while n < len(alts):
				#loop through the multiple alt alleles
				if clnSig != "NA":
					info_value[1] = clnSig.split(',')[min(n,clnSigmaxN)]
					#replace sig tag with values
					sig_keys = re.findall('[0-9]+', info_value[1])
					sig=''
					for sig_key in sig_keys:
						sig = sig + '|' + cln_dic[sig_key]
					info_value[1] = sig
				if caf != "NA":
					info_value[5] = caf.split(',')[0] + ',' + caf.split(',')[n+1]
				if  clnAllele != "NA":
					info_value[4] = clnAllele.split(',')[min(n,clnAlleleMaxN)]
				if clnDb != "NA":
					info_value[2] = clnDb.split(',')[min(n, clnDbMaxN)]
				if clnDisease != "NA":
					info_value[3] = clnDisease.split(',')[min(n,clnDiseaseMaxN)]
				outfile.write('\t'.join(line[0:4]) + '\t' + alts[n] + '\t')
				outfile.write('\t'.join(info_value) + '\n')
				n += 1
	infile.close()
	outfile.close()


def eff_to_db(infile, outfile):
	sig_dic = {"HIGH":1, "MODERATE":2, "LOW":3, "MODIFIER":4}
	infile = open(infile, 'r')
	outfile = open(outfile, 'w')
	outfile.write('#' + ' '.join(sys.argv) + '\n')
	outfile.write('@chr\tpos\tsnpID\tref\talt\tfrequency\ttype\tgenotype\tmulti\tfixed\tbase_corrected\tsample\tgeneID\tchangeType\tcDNAchange\tproteinChange\n')
	for line in infile:
		if line.startswith('#') or line.startswith('@'):
			continue
		else:
			line = line.split()
			genotype=line[7].split(';')[0]
			effect=line[7].split(';')[1]
			effect=effect.split('=')[1]
			effect_list = effect.split(',')
			effect_dic = {}
			for effect in effect_list:
				#put all effects into a dictionary
				#with geneID as dictionary keys
				#each effect is a list in a matrix
				effect = effect.split('|')
				geneID=effect[4]
				if geneID in effect_dic:
					effect_dic[geneID].append([sig_dic[effect[2]], geneID, effect[1], effect[9], effect[10]])
				else:
					effect_dic[geneID] = []
					effect_dic[geneID].append([sig_dic[effect[2]], geneID, effect[1], effect[9], effect[10]])
			output=[]
			for geneID in effect_dic:
				#select the most severe effect from each geneID
				effect_dic[geneID].sort(key=itemgetter(0))
				output.append(effect_dic[geneID][0])
			#select most severe effect among genes
			#when multiple genes are equally severe, output them all
			output.sort(key=itemgetter(0))
			effect_out = map(str,output[0])
			outfile.write('\t'.join(line[0:7]) + '\t' + genotype + '\t' + '\t'.join(line[8:]) + '\t' + '\t'.join(effect_out[1:]) + '\n')
			n = 1
			while n < len(output):
				if output[n][0] == output[0][0]:
					effect_out = map(str,output[n])
					outfile.write('\t'.join(line[0:7]) + '\t' + genotype + '\t' + '\t'.join(line[8:]) + '\t' + '\t'.join(effect_out[1:]) + '\n')
				else:
					break
				n += 1
	infile.close()
	outfile.close()
	
			
	
try:
	if sys.argv[1] == 'eff_to_db':
		parser = optparse.OptionParser()
		parser.add_option('-I', '--infile',
		help = 'input file')
		parser.add_option('-o', '--outfile',
		help = 'output file')
		options,args = parser.parse_args()
		infile = options.infile
		outfile = options.outfile
		eff_to_db(infile, outfile)
	if sys.argv[1] == 'brca_her_load':
		parser = optparse.OptionParser()
		parser.add_option('-I', '--infile',
		help = 'input file')
		parser.add_option('-o', '--outfile',
		help = 'output file')
		parser.add_option('-r', '--inref',
		help = 'input reference vcf')
		parser.add_option('-s', '--sampleID',
		help = 'sample ID in the output')
		options,args = parser.parse_args()
		infile = options.infile
		outfile = options.outfile
		inref = options.inref
		sampleID = options.sampleID
		brca_her_load(infile, outfile, inref, sampleID)
	if sys.argv[1] == 'clinvar_load':
		parser = optparse.OptionParser()
		parser.add_option('-I', '--infile',
		help = 'input file')
		parser.add_option('-o', '--outfile',
		help = 'output file')
		options,args = parser.parse_args()
		infile = options.infile
		outfile = options.outfile
		clinvar_load(infile, outfile)
except IndexError:
		print 'Functions:'
		print ' brca_her_load(): take an input vcf file and output a file for database import'
		print ' clinvar_load(): convert clinvar vcf file into database load format'
		print ' eff_to_db(): take the snpEff output and produce a database loading file'