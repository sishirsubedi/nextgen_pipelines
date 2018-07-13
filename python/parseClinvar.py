import commonTool
import sys
import optparse

sao={"0":"Unspecified", "1":"Germline", "2":"Somatic", "3":"Both"}
sig = {"0":"Uncertain", "1":"Not-Provided", "2":"Benign", "3":"Likely-Benign", "4":"Likely-pathogenic", "5":"Pathogenic", "6":"Drug-response", "7":"Histocompatibility", "255":"Other"}

def parseClinvar(infile, outfile):
	infile = open(infile, 'r')
	outfile = open(outfile, 'w')
	for line in infile:
		if line.startswith("#"):
			outfile.write(line)
		else:
			record = line.split()
			if not ("," in record[4] or "," in record[6]):
				##single alt allele, and single clnical sig allele
				if record[6] == "0":
					record[6] = "ref"
				elif record[6] == "1":
					record[6] = "alt"
				record = parseTag(record)
				outfile.write("\t".join(record))
				outfile.write("\n")
			else:
				altAlleles = record[4].split(",")
				clnAlleles = record[6].split(",")
				clnSig = record[7].split(",")
				clnAcc = record[8].split(",")
				n = 0
				for clnAllele in clnAlleles:
					if clnAllele == "0":
						index = 0
						allele = "ref"
					elif clnAllele == "-1":
						n += 1
						continue
					else:
						index = int(clnAllele) - 1
						allele = "alt"
					try:
						output= [record[0], record[1], record[2], record[3],altAlleles[index], record[5], allele, clnSig[n], clnAcc[n]]
					except IndexError:
						print line
						sys.exit()
					output = parseTag(output)
					outfile.write("\t".join(output))
					outfile.write("\n")
					n += 1
	infile.close()
	outfile.close()
	
def parseTag(record):
	record[5] = sao[record[5]]
	clinsig = record[7].split("|")
	wordSig = []
	for numberSig in clinsig:
		try:
			wordSig.append(sig[numberSig])
		except KeyError:
			print record
			sys.exit()
	record[7] = "|".join(wordSig)
	return record
		

if len(sys.argv) <= 1:
	print "For usage: python parseClinvar.py -h "
else:
	parser = optparse.OptionParser()
	parser.add_option('-I', '--infile',
	help = 'input vcf file')
	parser.add_option('-o', '--outfile',
	help = 'output file')
	options,args = parser.parse_args()
	infile = options.infile
	outfile = options.outfile
	parseClinvar(infile, outfile)