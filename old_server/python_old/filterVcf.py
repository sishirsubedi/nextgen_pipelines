import commonTool
import sys
import optparse
import re
import operator


def filterDetail(vcfInfile, vcfOutfile, sampleCol, field, operation, value):
	infile = open(vcfInfile, 'r')
	outfile = open(vcfOutfile, 'w')
	
	for line in infile:
		if line.startswith('#'):
			outfile.write(line)
		else:
			oneVcf = commonTool.vcf(line.split('\t'))
			fieldValue = oneVcf.parseDetail(field, sampleCol)
			if fieldValue is None:
				infile.close()
				outfile.close()
				return
			if comparison(fieldValue, operation, value) is None:
				print "Error: Aborted"
				return
			elif comparison(fieldValue, operation, value):
				outfile.write(line)
	infile.close()
	outfile.close()
	
def filterInfo(vcfInfile, vcfOutfile, field, operation, value):
	infile = open(vcfInfile, 'r')
	outfile = open(vcfOutfile, 'w')
	
	for line in infile:
		if line.startswith('#'):
			outfile.write(line)
		else:
			oneVcf = commonTool.vcf(line.split('\t'))
			fieldValue = oneVcf.parseInfo(field)
			if fieldValue is None:
				infile.close()
				outfile.close()
				return
			if comparison(fieldValue, operation, value) is None:
				print line
				print "Error: Aborted"
				return
			elif comparison(fieldValue, operation, value):
				outfile.write(line)
	infile.close()
	outfile.close()
			
			

def comparison(fieldValue, operation, value):
	#perform comparison
	#return boolean
	#print operation
	ops = {
			"==": operator.eq,
			"!=": operator.ne,
			"<": operator.lt,
			"<=": operator.le,
			">": operator.gt,
			">=": operator.ge
			}
	if operation not in ops:
		print "specified operation %s not supported" %operation
		return None
	if operation in [">=", ">", "<=", "<"]:
		#print "number comparison"
		try:
			fieldValue = float(fieldValue)
			value = float(value)
		except ValueError:
			print "specified operation %s not valid for %s or %s" %(operation, str(fieldValue), str(value))
			return None
	if ops[operation](fieldValue, value):
		return True
	else:
		return False


if len(sys.argv) > 1:
	parser = optparse.OptionParser()
	parser.add_option('-I', '--vcfInfile',
	help = 'input vcf file')
	parser.add_option('-o', '--vcfOutfile',
	help = 'output file')
	parser.add_option('-n', '--sampleCol',
	help = 'sample column number')
	parser.add_option('-i', '--info', action="store_true", 
	help = 'filter on the info column, not compatible with -n')
	parser.add_option('-f', '--field',
	help = 'field string')
	parser.add_option('-p', '--operation',
	help = 'operations, choose between: >, >=, <, <=, ==, !=; please use double quote to avoid bash confusion')
	parser.add_option('-v', '--value',
	help = 'comparison value')
	options,args = parser.parse_args()
	vcfInfile = options.vcfInfile
	vcfOutfile = options.vcfOutfile
	info = options.info
	sampleCol = options.sampleCol
	field = options.field
	operation = options.operation
	value = options.value
	if sampleCol is None:
		filterInfo(vcfInfile, vcfOutfile, field, operation, value)
	elif info is None:
		filterDetail(vcfInfile, vcfOutfile, sampleCol, field, operation, value)
	else:
		print "Error: -i and -n are mutually exclusive"
	
else:
	print "filter vcf file according to the specified individual detail column"
	print "-h for more help" 