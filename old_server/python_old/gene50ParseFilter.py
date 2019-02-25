import sys
import optparse
import MySQLdb

def parseMutation(input, output, reportNo, barcode,coverage):
	input = open(input, 'r')
	output = open(output, 'w')
	coverage = open(coverage, 'r')
	db = MySQLdb.connect("localhost", "niyunyun", "molSeq3127")
	c = db.cursor()
	command="select * from molseq.gene50_sample_view where report_no = '%s' and barcode = '%s'" %(reportNo, barcode)
	c.execute(command)
	result = c.fetchone()
	if result:
		lastName=result[2]
		firstName=result[3]
		orderNumber = result[4]
		pathNumber = result[5]
		command="select * from molseq.gene50_sample_view where lastName = '%s' and firstName = '%s' and report_no != '%s' and barcode != '%s'" %(lastName, firstName, reportNo, barcode)
		c.execute(command)
		backTracking = c.fetchall()
	else:
		lastName = ''
		firstName = ''
		orderNumber = ''
		pathNumber = ''
		backTracking = None
	output.write("Name: %s,%s\n" %(lastName, firstName))
	output.write("Order Number: %s\n" %orderNumber)
	output.write("Pathology Number: %s\n" %pathNumber)
	if backTracking:
		output.write("This name was tested in the following previous runs\n")
		output.write("reportNumber\tbarcode\torderNumber\tpathologyNumber\n")
		for record in backTracking:
			record = map(str, list(record))
			output.write("%s\t%s\t%s\t%s\n" %(record[0], record[1], record[4], record[5]))
	mutation = 1
	for line in input:
		output.write("#####Mutation %i #####\n"%mutation)
		mutation += 1
		line = line.strip()
		line = line.split('\t')
		annotation = line[20] + ' ' + line[21]
		output.write("Mutation Info: %s:%s,%s\n" %(line[17], line[14], line[15]))
		output.write("Cosmic ID: %s\n" %line[3])
		output.write("dbSNP ID: %s\n" %line[18])
		output.write("Coordinate: %s:%s\n" %(line[1], line[2]))
		output.write("Genotype: %s; Ref: %s; Alt: %s\n" %(line[6], line[4], line[5]))
		output.write("Frequency: %s(%s/%s)\n" %(line[9], line[8], line[7]))
		output.write("Occurances: %s\n" %line[22])
		output.write("Previous annotations:\n")
		annotation = annotation.split()
		characters = 0
		for n in range(0, len(annotation)):
			characters += len(annotation[n])
			characters += 1
			if characters <= 60:
				output.write(annotation[n] + ' ')
			else:
				output.write("\n")
				output.write(annotation[n] + ' ')
				characters = len(annotation[n])
		output.write("\n\n\n")
	output.write('Insufficiently Covered Amplicons\n')
	for line in coverage:
		output.write(line)
	coverage.close()
	input.close()
	output.close()
	
	
if len(sys.argv) > 1:
	parser = optparse.OptionParser()
	parser.add_option('-I', '--input',
	help = 'input filtered mutation file')
	parser.add_option('-o', '--output',
	help = 'output file')
	parser.add_option('-r', '--reportNo',
	help = 'report no')
	parser.add_option('-b', '--barcode',
	help = 'barcode')
	parser.add_option('-c', '--coverage',
	help = 'file with insufficiently covered amplicons')
	options,args = parser.parse_args()
	input = options.input
	output = options.output
	reportNo = options.reportNo
	barcode = options.barcode
	coverage= options.coverage
	parseMutation(input, output, reportNo, barcode, coverage)
else:
	print "python gene50ParseFilter.py -h for help"