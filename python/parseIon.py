import commonTool
import sys
import optparse

def parseIon(infile, outfile):
	infile = open(infile, 'r')
	outfile = open(outfile, 'w')
	outfile.write("#chr\tpos\tid\tref\talt\tAF\tFAO\tFDP\n")
	for line in infile:
		if not line.startswith("#"):
			record = line.split("\t")
			if not "," in record[4]:
				vcf = commonTool.vcf(record)
				AF = vcf.parseInfo("AF")
				FAO = vcf.parseInfo("FAO")
				FDP = vcf.parseInfo("FDP")
				output = [record[0], record[1], record[2], record[3], record[4], AF, FAO, FDP]
				outfile.write("\t".join(output))
				outfile.write("\n")
			else:
				vcf = commonTool.vcf(record)
				altAlleles = record[4].split(",")
				AF = vcf.parseInfo("AF").split(",")
				FAO = vcf.parseInfo("FAO").split(",")
				FDP = vcf.parseInfo("FDP")
				n = 0
				while n < len(altAlleles):
					try:
						newRecord = [record[0], record[1], record[2], record[3], altAlleles[n], AF[n], FAO[n], FDP]
					except IndexError:
						print line
						sys.exit()
					newVcf = commonTool.vcf(newRecord)
					output = newVcf.leftAlign()
					outfile.write("\t".join(output))
					outfile.write("\n")
					n += 1
	infile.close()
	outfile.close()
					


if len(sys.argv) <= 1:
	print "For usage: python parseIon.py -h "
else:
	parser = optparse.OptionParser()
	parser.add_option('-I', '--infile',
	help = 'input vcf file')
	parser.add_option('-o', '--outfile',
	help = 'output file')
	options,args = parser.parse_args()
	infile = options.infile
	outfile = options.outfile
	parseIon(infile, outfile)