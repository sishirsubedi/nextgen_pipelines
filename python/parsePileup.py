import sys
import optparse
import re
import scipy.stats as stats
import commonTool

			
def parsePileup(infile, outfile):
	infile = open(infile, 'r')
	outfile = open(outfile, 'w')
	outfile.write("#chr\tpos\tdepth\tref+\tref-\talt+\talt-\tfisherBiasP\n")
	for line in infile:
		if not line.startswith("#"):
			record = line.split()
			altAlleles = record[4].split(",")
			oneVcf = commonTool.vcf(record)
			forwardCounts = oneVcf.parseInfo("ADF").split(",")
			refForward = forwardCounts[0]
			reverseCounts = oneVcf.parseInfo("ADR").split(",")
			refReverse = reverseCounts[0]
			depth = oneVcf.parseInfo("DP")
			for alt in altAlleles:
				if not alt == "<*>" or len(altAlleles) == 1:
					index = altAlleles.index(alt) + 1
					forward = forwardCounts[index]
					reverse = reverseCounts[index]
					biasRatio, biasP = stats.fisher_exact([[int(refForward), int(forward)],[int(refReverse), int(reverse)]])
					coordinate = [record[0], record[1], ".", record[3], alt]
					coordinateLeft = commonTool.vcf(coordinate).leftAlign()
					output = [coordinateLeft[0], coordinateLeft[1], coordinateLeft[3], coordinateLeft[4], depth, refForward, refReverse, forward, reverse, str(biasP)]
					outfile.write("\t".join(output))
					outfile.write("\n")
	infile.close()
	outfile.close()
	
if len(sys.argv) > 1:
	parser = optparse.OptionParser()
	parser.add_option('-I', '--vcfInfile',
	help = 'input vcf file')
	parser.add_option('-o', '--outfile',
	help = 'output file')
	options,args = parser.parse_args()
	infile = options.vcfInfile
	outfile = options.outfile
	parsePileup(infile, outfile)
	
else:
	print "for usage: python parsePipeup.py -h"
					