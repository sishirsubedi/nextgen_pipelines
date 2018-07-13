import sys
import optparse
import re
from operator import itemgetter
from vcfdbLoad import correct_base
from subprocess import check_output
from subprocess import call

#define class
class vcf(object):
	def __init__(self, record):
		self.record = record.split()
		self.info = self.record[7].split(';')
		self.info_dic = {}
		for record in self.info:
			#put info column in dictionary
			record = record.split('=')
			if not record[0] in self.info_dic and len(record) > 1:
				self.info_dic[record[0]] = record[1]
			else:
				continue
	def chr(self):
		if self.record[0].startswith('chr'):
			return self.record[0]
		else:
			return 'chr'+self.record[0]
	def pos(self):
		return self.record[1]
	def id(self):
		#list of ids
		return self.record[2].split(';')
	def ref(self):
		return self.record[3]
	def alt(self):
		#list of alt alleles
		return self.record[4].split(',')
	def fao(self):
		#list of alt allele counts
		return self.info_dic["FAO"].split(',')
	def fdp(self):
		#read depth
		return self.info_dic["FDP"]
	def genotype(self):
		#genotype
		return self.record[9].split(':')[0]

	
class cosmic(object):
	def __init__(self, record):
		self.record = record.split()
		self.info = self.record[7].split(';')
		self.info_dic = {}
		for record in self.info:
			#put info column in dictionary
			record = record.split('=')
			if not record[0] in self.info_dic:
				self.info_dic[record[0]] = record[1]
			else:
				continue
	def chr(self):
		if self.record[0].startswith('chr'):
			return self.record[0]
		else:
			return 'chr'+self.record[0]
	def pos(self):
		return self.record[1]
	def ref(self):
		return self.record[3]
	def alt(self):
		return self.record[4]
	def gene(self):
		return self.info_dic["GENE"]
	def cds(self):
		return self.info_dic["CDS"]
	def aa(self):
		return self.info_dic["AA"]


		
def vcf_to_db(oneVcf):
	#oneVcf is vcf class object
	multi = "N"  #record if it is a multiple cosmic on one coordinate
 	vcf_db=[]  # output matrix, one uniq cosmic per row
	n = 0
	genotype = oneVcf.genotype()
	ref = oneVcf.ref()
	altAlleles = oneVcf.alt()  #list of alt alleles
	altAlleleCounts = oneVcf.fao()  #list of alt allele counts
	ids = oneVcf.id()     #list of ids
	if len(altAlleles) > 1:
		multi = "Y"
	for alt in altAlleles:
		try:
			#get some info for the current alt allele
			#if number of alt alleles exceeds number of ids or number of alt allele counts
			#stop and return the current records
			id = ids[n]
			altAlleleCount = altAlleleCounts[n]
		except IndexError:
			break
		if altAlleleCount == "0":
			#no counts for this alt allele
			n += 1
			continue
		#current alt allele has reads and id
		#change ref and alt column for ins/del to match cosmic
		base_corrected = "N"
		if (len(ref) > 1 and len(alt) > 1):
			(ref, alt, base_corrected) = correct_base(ref, alt)
		frequency=float(altAlleleCount)/float(oneVcf.fdp())
		vcf_db.append([oneVcf.chr(),oneVcf.pos(),id,ref,alt,genotype, altAlleleCount, oneVcf.fdp(), frequency,multi,base_corrected])
		n += 1
	return vcf_db
	
def gene50vcfLoad(vcfInfile, dbOutfile):
	vcfInfile = open(vcfInfile, 'r')
	dbOutfile = open(dbOutfile, 'w')
	for line in vcfInfile:
		if line.startswith('#'):
			continue
		else:
			
			oneVcf = vcf(line)
			genotype = oneVcf.genotype()
			if re.match('0.0', genotype) or re.match('\..\.', genotype):
				#if genotype is wildtype or not called, continue
				continue
			else:
				#convert vcf to dbload format as a matrix
				vcf_db = vcf_to_db(oneVcf)
				for record in vcf_db:
					record = map(str, record)
					dbOutfile.write('\t'.join(record) + '\n')
	vcfInfile.close()
	dbOutfile.close()
	
def cosmic_to_dic(cosmic):
	infile = open(cosmic, 'r')
	cosmic_dic = {}
	print "reading cosmic file"
	for line in infile:
		if not line.startswith('#'):
			line = line.split('\t')
			chr = line[0]
			if not chr.startswith('chr'):
				chr = "chr" + chr
			pos = line[1]
			id = line[2]
			ref = line[3]
			alt = line[4]
			info = line[7].split(';')
			info_dic = {}
			for record in info:
				record = record.split('=')
				try:
					info_dic[record[0]] = record[1]
				except IndexError:
					#no "=" in this information
					continue
			try:
				gene = info_dic["GENE"]
			except KeyError:
				gene = "."
			try:
				cDNA = info_dic["CDS"]
			except KeyError:
				cDNA = "."
			try:
				protein = info_dic["AA"]
			except KeyError:
				protein = "."
			cosmic_dic[id] = [chr, pos, ref, alt, cDNA, protein, gene]
	infile.close()
	print 'finished reading cosmic'
	return cosmic_dic
	
			
def effect(no_cosmic, outfile, report_no, barcode):
	print 'getting effect info for novel mutations'
	temp = open('/home/niyunyun/temp/temp.txt', 'w')
	record_dic = {}  #mutation record
	for record in no_cosmic:
		temp.write('\t'.join([record[0], record[1], record[3], record[4]]) + '\n')
		record_dic[(record[0], record[1], record[3], record[4])] = record  #mutation record
	temp.close()
	print 'running snpeff'
	tempEff = open('/home/niyunyun/temp/temp.eff.txt', 'w')
	command = ['java', '-Xmx4g', '-jar', '/opt/software/snpEff-4.1b/snpEff.jar', 'GRCh37.75', '/home/niyunyun/temp/temp.txt']
	call(command, stdout=tempEff)
	tempEff.close()
	print 'finished running snpeff'
	effect_file = open('/home/niyunyun/temp/temp.eff.txt', 'r')
	for line in effect_file:
		if line.startswith('#'):
			continue
		else:
			line = line.split('\t')
			eff = line[7]
			effects = eff.split('=')[1].split(',')  #list of effects
			sig_dic = {"HIGH":1, "MODERATE":2, "LOW":3, "MODIFIER":4}
			eff_matrix = []  #each effect is a list in this matrix
			for effect in effects:
				effect = effect.split('|')
				effect[2] = sig_dic[effect[2]]  # convert effect significance to a number to be sorted later
				eff_matrix.append(effect)
			eff_matrix.sort(key=itemgetter(2))
			sig_effect = eff_matrix[0]  #most severe effect
			if sig_effect[9] == "":
				sig_effect[9] = "."
			if sig_effect[10] == "":
				sig_effect[10] = "."
			if sig_effect[3] == "":
				sig_effect[3] = "."
			mutation = record_dic[(line[0], line[1], line[2], line[3])]  #get full mutation record for this mutation
			output =  mutation + [sig_effect[3], sig_effect[9], sig_effect[10]] + [report_no, barcode] 
			outfile.write('\t'.join(output) + '\n')
	return(outfile)
	
	
	


'''def gene50excelLoad(infile, cosmic, outfile, report_no, barcode):
	infile = open(infile, 'r')
	cosmic_dic = cosmic_to_dic(cosmic)
	outfile = open(outfile, 'w')
	no_cosmic = []
	for line in infile:
		if line.startswith('Chrom'):
			#header
			continue
		else:
			line = line.split('\t')
			if line[4] == "Heterozygous" or line[4] == "Homozygous":
				#called
				if line[11] != "---":
					#has a cosmic ID, use the cosmic coordinates
					chr = cosmic_dic[line[11]][0]
					pos = cosmic_dic[line[11]][1]
					ref = cosmic_dic[line[11]][2]
					alt = cosmic_dic[line[11]][3]
					cDNA = cosmic_dic[line[11]][4]   # cDNA change
					protein = cosmic_dic[line[11]][5]   #codon change
					gene = cosmic_dic[line[11]][6]     #gene
					depth = line[18]      #read depth
					altCount = line[24]     #alt allele count
					frequency = line[6]    #alt allele frequency
					genotype = line[4]
					cosmic = line[11]
					outfile.write('\t'.join([chr, pos, cosmic, ref, alt, genotype, depth, altCount, frequency, gene, cDNA, protein,report_no, barcode]) + '\n')
				else:
					#no cosmic ID, use the vcf coordinates
					chr = line[0]
					pos = line[14]
					ref = line[15]
					alt = line[16]
					depth = line[18]
					altCount = line[24]
					frequency = line[6]
					genotype = line[4]
					cosmic = line[11]
					#correct base if necessary
					(ref,alt, base_corrected) = correct_base(ref,alt)
					no_cosmic.append([chr, pos, cosmic, ref, alt, genotype, depth, altCount, frequency])
	#get effect of all mutations without cosmic ID
	outfile = effect(no_cosmic, outfile, report_no, barcode)
	infile.close()
	outfile.close()'''
	
	
try:
	if sys.argv[1] == "gene50vcfLoad":
		parser = optparse.OptionParser()
		parser.add_option('-I', '--vcfInfile',
		help = 'input vcf file')
		parser.add_option('-o', '--dbOutfile',
		help = 'output file')
		options,args = parser.parse_args()
		vcfInfile = options.vcfInfile
		dbOutfile = options.dbOutfile
		gene50Load(vcfInfile, dbOutfile)
		
	if sys.argv[1] == 'gene50excelLoad':
		parser = optparse.OptionParser()
		parser.add_option('-I', '--infile',
		help = 'input excel file')
		parser.add_option('-o', '--outfile',
		help = 'output file')
		parser.add_option('-c', '--cosmic',
		help = 'input cosmic file')
		parser.add_option('-r', '--report_no',
		help = 'ion torrent report number')
		parser.add_option('-b', '--barcode',
		help = 'sequencing barcode (three digits)')
		options,args = parser.parse_args()
		infile = options.infile
		outfile = options.outfile
		cosmic = options.cosmic
		report_no = options.report_no
		barcode = options.barcode
		gene50excelLoad(infile, cosmic, outfile, report_no, barcode)

except ValueError:
	print 'Available functions:'
	print ' gene50vcfLoad(): take vcf file from 50 gene results and output a file for dbloading'
	print ' gene50excelLoad(): take excel output of ion torrent and output a file for dbloading'
	
