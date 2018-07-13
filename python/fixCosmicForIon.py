import sys
import optparse
import re
import commonTool



class fixCosmic:
	def __init__(self,cosmic,outfile):
		self.infile = cosmic
		self.outfile = outfile
	def fix(self):
		infile = self.infile
		outfile = self.outfile
		infile = open(infile, 'r')
		outfile = open(outfile, 'w')
		outfile.write("##fileformat=VCFv4.1\n")
		outfile.write("##allowBlockSubstitutions=true\n")
		outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
		currentPosition = ('','')
		records={}
		for line in infile:
			if line.startswith('#'):
				continue
			else:
				record = line.split('\t')
				if (record[0], record[1]) != currentPosition:
					#new position
					self.writeOut(records, outfile)
					currentPosition = (record[0], record[1])
					records = {}
					records = self.append(records, record)
				else:
					records = self.append(records, record)
		infile.close()
		outfile.close()
	
	
	def writeOut(self, records, outfile):
		for key in records:
			gene = key[4]
			#if transcript variant, check if the same mutation is represented in full length transcript
			if re.search('_', gene):
				geneName = gene.split('_')[0]
				if (key[0], key[1], key[2], key[3], geneName) in records:
					continue
			for record in records[key]:
				out = '\t'.join(record)
				outfile.write(out)
			
			
	
	def append(self, records, record):
		gene = commonTool.vcf(record).parseInfo("GENE")
		key = (record[0], record[1], record[3], record[4], gene)
		if key in records:
			records[key].append(record)
		else:
			records[key] = [record]
		return records
		











try:
	parser = optparse.OptionParser()
	parser.add_option('-I', '--cosmic',
	help = 'input cosmic vcf file')
	parser.add_option('-o', '--outfile',
	help = 'output file')
	options,args = parser.parse_args()
	cosmic = options.cosmic
	outfile = options.outfile
	result = fixCosmic(cosmic, outfile)
	result.fix()
		
	

except TypeError:
	print 'This script removes variant records annotated to be within a transcript variant' 
	print 'when the same variant is already represented in a canonical transcript'
	print 'type -h for running help'