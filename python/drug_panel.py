import sys
import optparse



def cutOfftolist(cutOff):
	cutOff = open(cutOff, 'r')
	cutOff_list = []
	for line in cutOff:
		if line.startswith('#'):
			continue
		else:
			line = line.split('\t')
			if len(line) > 1:
				cutOff_list.append((line[0],float(line[1])))
	return cutOff_list
	
def select_positive(data1,data2, cutOff, outfile):
	cutOff_list = cutOfftolist(cutOff)
	data1 = open(data1, 'r')
	data2 = open(data2, 'r')
	n = 1
	outfile = open(outfile, 'w')
	outfile.write('@sampleName\tsampleID\tsampleType\tassay\tValue\tCutoff\n')
	positive = 0
	for line1 in data1:
		line2 = data2.readline()
		if positive == 1:
			#previous sample has positive records
			#output a empty line
			outfile.write('\n')
			positive = 0
		if n == 1:
			n += 1
			continue
		else:
			n +=1
			line1 = line1.split(',')
			line2 = line2.split(',')
			sample = line1[0]
			id = line1[1]
			type = line1[2]
			records = line1[3:] + line2[3:]
			i = 0
			for record in records:
				try:
					record = float(record)
				except ValueError:
					record = 0
				#print 'value is ' + str(record)
				#print 'cutoff is ' + str(cutOff_list[
				if record >= cutOff_list[i][1]:
					#meet cut off 
					positive = 1
					outfile.write('\t'.join([sample, id, type, cutOff_list[i][0], str(record), str(cutOff_list[i][1])]) + '\n')
				i += 1
	data1.close()
	data2.close()
	outfile.close()
	
try:
	if sys.argv[1] == 'select_positive':
		parser = optparse.OptionParser()
		parser.add_option('--data1',
		help = 'input positive data file')
		parser.add_option('--data2',
		help = 'input negative data file')
		parser.add_option('-c', '--cutOff',
		help = 'intput cutoff file')
		parser.add_option('-o', '--outfile',
		help = 'output file')
		options,args = parser.parse_args()
		data1 = options.data1
		data2 = options.data2
		cutOff = options.cutOff
		outfile = options.outfile
		select_positive(data1,data2, cutOff, outfile)
except IndexError:
		print 'Functions:'
		print ' select_positive(): Take input data file and cutoff file, output positive results '