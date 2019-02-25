import commonTool
import sys
import optparse


if len(sys.argv) < 2:
	print("Usage:")
	print("python parseVcfInfo.py -I input.vcf -o output.vcf -f 'comma seperated names of info to parse'")
	print("For more info: python parseVcfInfo.py -h")
	sys.exit()
	
def parseInfo(infile, outfile, fields):
	infile = open(infile, 'r')
	outfile = open(outfile, 'w')
	outfile.write("#chr\tpos\tid\tref\talt\t")
	for field in fields:
		outfile.write("%s\t" %field)
	outfile.write("\n")
	for line in infile:
		if not line.startswith('#'):
			oneVcf = line.split('\t')
			if not oneVcf[0].startswith("chr"):
				oneVcf[0] = "chr" + oneVcf[0]
			altAlleles = oneVcf[4].split(",")
			out = oneVcf[0:4]
			n = 0
			while n < len(altAlleles):
				out.append(altAlleles[n])
				for field in fields:
					info = commonTool.vcf(oneVcf).parseInfo(field)
					if "," in info:
						out.append(info.split(",")[n])
					else:
						out.append(info)
				n += 1
			outfile.write('\t'.join(out))
			outfile.write('\n')
	outfile.close()

if __name__ == "__main__":
	parser = optparse.OptionParser()
	parser.add_option('-I', '--infile',
	help = 'input vcf file')
	parser.add_option('-o', '--outfile',
	help = 'output file')
	parser.add_option('-f', '--fields',
	help = 'comma seperated names of info')
	options,args = parser.parse_args()
	infile = options.infile
	outfile = options.outfile
	fields = options.fields
	fields = fields.split(',')
	parseInfo(infile, outfile, fields)
	

	

		
				


	