import sys
import datetime
import re

time = datetime.datetime.now().strftime("%m-%d-%Y-%H%M")
infile = sys.argv[1]
outfile = sys.argv[2]

infile = open(infile, 'r')
outfile = open(outfile, 'w')
outfile.write("#time\twell\tsample\tdetector\ttask\tct\tresult\n")

n = 1

for line in infile:
	if n > 28:
		#start analysis
		line = line.strip()
		line = line.split(',')[0:5]
		try:
			ct = float(line[4])
		except ValueError:
			ct = 100.00
		if line[2] == "ABL_JOE":
			if ct <= 34.0:
				line.append('positive')
			else:
				line.append('negative')
		if line[2] == "ALK_Qiagen":
			if ct <= 32.0:
				line.append('positive')
			else:
				line.append('negative')
		outfile.write(time + '\t' + '\t'.join(line) + '\n')
	n += 1
	
infile.close()
outfile.close()
