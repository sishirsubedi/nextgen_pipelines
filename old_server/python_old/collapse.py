import sys
import optparse

def collapse(infile, outfile, col, separator):
	infile = open(infile, 'r')
	outfile = open(outfile, 'w')
	col = col.split(',')
	col = map(int, col)
	col[:] = [x - 1 for x in col]
	key = ["NA"]*len(col)
	output = []
	for line in infile:
		if line.startswith("#"):
			outfile.write(line)
		else:
			record = line.strip().split("\t")
			newKey = []
			for n in col:
				newKey.append(record[n])
			if newKey != key and output != []:
				key = newKey
				outfile.write("\t".join(output))
				outfile.write("\n")
				output = record
			elif newKey != key and output == []:
				key = newKey
				output = record
			else:
				n = 0
				while n < len(record):
					if not n in col:
						try:
							output[n] = output[n] + separator + record[n]
						except IndexError:
							print output
							print record
							sys.exit()
					n += 1
	infile.close()
	outfile.close()





if len(sys.argv) > 1:
	parser = optparse.OptionParser()
	parser.add_option('-I', '--infile',
	help = 'input file')
	parser.add_option('-o', '--outfile',
	help = 'output file')
	parser.add_option('-c', '--col',
	help = 'comma seperated column numbers to collapse')
	parser.add_option('-s', '--separator',
	help = 'separator')
	options,args = parser.parse_args()
	infile = options.infile
	outfile = options.outfile
	col = options.col
	separator = options.separator
	collapse(infile, outfile, col, separator)
else:
	print "Combine rows of records based on selected columns"
	print "File must be sorted on the selected columns"
	print "For Usage: python collapse.py -h"