import sys
import optparse
import re

def cosmic2db(cosmic, outfile):
	cosmic = open(cosmic, 'r')
	outfile = open(outfile, 'w')
	for line in cosmic:
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
			synonymousP = re.compile(r"^p\.([A-Z])[0-9]+\1$")
			nonsynonymousP = re.compile(r"^p\.([A-Z])[0-9]+[A-Z]$")
			if synonymousP.search(protein):
				coding = "synonymous"
			elif nonsynonymousP.search(protein):
				coding = "non-synonymous"
			elif protein == "p.?":
				coding = "non-coding"
			else:
				coding = "other"
			outfile.write('\t'.join([chr, pos, id, ref, alt, cDNA, protein, coding, gene]) + '\n')
	cosmic.close()
	outfile.close()
	
def dbSNP2db(dbsnp, outfile):
	dbsnp = open(dbsnp, 'r')
	outfile = open(outfile, 'w')
	for line in dbsnp:
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
				property = info_dic["SAO"]
				if property == "0":
					property = "unspecified"
				if property == "1":
					property = "germline"
				if property == "2":
					property = "somatic"
				if property == "3":
					property = "both"
			except KeyError:
				property = "."
			outfile.write('\t'.join([chr, pos, id, ref, alt, property]) + '\n')
	dbsnp.close()
	outfile.close()
	
	
try:
	if sys.argv[1] == 'cosmic2db':
		parser = optparse.OptionParser()
		parser.add_option('-I', '--cosmic',
		help = 'input cosmic vcf file')
		parser.add_option('-o', '--outfile',
		help = 'output file')
		options,args = parser.parse_args()
		cosmic = options.cosmic
		outfile = options.outfile
		cosmic2db(cosmic, outfile)
		
	if sys.argv[1] == 'dbSNP2db':
		parser = optparse.OptionParser()
		parser.add_option('-I', '--dbsnp',
		help = 'input cosmic vcf file')
		parser.add_option('-o', '--outfile',
		help = 'output file')
		options,args = parser.parse_args()
		dbsnp = options.dbsnp
		outfile = options.outfile
		dbSNP2db(dbsnp, outfile)

except IndexError:
	print 'Available functions:'
	print ' cosmic2db(): convert cosmic vcf to db loading file'
	print ' dbSNP2db(): convert dbSNP vcf to db loading file'

