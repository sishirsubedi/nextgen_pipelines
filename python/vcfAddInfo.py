import commonTool
import sys
import optparse

def vcfAddInfo(infile, vcfData, outfile,info):
	vcfData = open(vcfData, 'r')
	infile = open(infile, 'r')
	outfile = open(outfile, 'w')
	vcfDic = {}
	for line in vcfData:
		if not line.startswith('#'):
			oneVcf = line.split('\t')
			if not oneVcf[0].startswith('chr'):
				oneVcf[0] = 'chr' + oneVcf[0]
			coordinate = "%s|%s|%s|%s" %(oneVcf[0], oneVcf[1], oneVcf[3], oneVcf[4])
			infoValue = commonTool.vcf(oneVcf).parseInfo(info)
			if coordinate in vcfDic:
				vcfDic[coordinate] += ",%s" %infoValue
			else:
				vcfDic[coordinate] = infoValue
	for line in infile:
		if line.startswith('#'):
			outfile.write(line)
		else:
			line = line.strip()
			oneVcf = line.split('\t')
			if not oneVcf[0].startswith('chr'):
				oneVcf[0] = 'chr' + oneVcf[0]
			coordinate = "%s|%s|%s|%s" %(oneVcf[0], oneVcf[1], oneVcf[3], oneVcf[4])
			if coordinate in vcfDic:
				oneVcf[7] += ";%s=%s" %(info, vcfDic[coordinate])
			else:
				oneVcf[7] += ";%s=None" %info
			outfile.write('\t'.join(oneVcf))
			outfile.write('\n')
			
	vcfData.close()
	infile.close()
	outfile.close()




if __name__ == "__main__":
	parser = optparse.OptionParser()
	parser.add_option('-I', '--infile',
	help = 'input vcf file to be annotated')
	parser.add_option('--vcfData',
	help = 'vcf database where the info is taken from')
	parser.add_option('--info',
	help = 'info field name in vcfData to add to infile')
	parser.add_option('-o', '--outfile',
	help = 'output vcf')
	options,args = parser.parse_args()
	infile = options.infile
	vcfData = options.vcfData
	outfile = options.outfile
	info = options.info
	vcfAddInfo(infile, vcfData, outfile, info)
	