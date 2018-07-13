import pysam
import re
import sys
import optparse
import os
from operator import itemgetter
import random

dedupedMatrix = []


def checkInput(infile, outfile, vcfFile, variantNumber, minFreq, maxFreq, freqFile):
	#sortedVcfFile = re.sub(r'\.vcf', r'.sort.vcf', vcfFile)
	#sortCommand = "grep -v '#' %s |sort -t$'\t' -k1,1 -k2,2n > %s" %(vcfFile, sortedVcfFile)
	#os.system(sortCommand)
	#vcf = open(sortedVcfFile, 'r')
	global dedupedMatrix
	if vcfFile != None:
		print "input vcf file"
		dedupedMatrix = dedupVcf(vcfFile)
	else:
		print "input frequency file"
		dedupedMatrix = dedupFreq(freqFile)
	print dedupedMatrix
	if len(dedupedMatrix) < int(variantNumber):
		print "Error: After removing overlapping variants, the specified variant number to simulate exceeds supplied variant number in file"
		return False
	if float(minFreq) < 0:
		print "Error: Minimum frequency is 0"
		return False
	if float(maxFreq) > 100:
		print "Error: Maximum frequency is 100"
		return False
	if freqFile != None:
		freqFileHandler = open(freqFile, 'r')
		for line in freqFileHandler:
			if not line.startswith('#'):
				record = line.split('\t')
				if float(record[4]) < 0 or float(record[4]) > 100:
					print "Error: allele frequency should be between 0 and 100 in the frequency file"
					return False
	
	return True
	'''if freqFile != None:
		freqFileHandler = open(freqFile, 'r')
		freqDic = {}
		for line in freqFileHandler:
			record = line.strip().split('\t')
			frequency = float(record[4])
			if not (record[0], record[1], record[2], record[3]) in freqDic:
				freqDic[(record[0], record[1], record[2], record[3])] = frequency
		freqFileHandler.close()
		sortedVcfFileHandler = open(sortedVcfFile, 'r')
		for line in sortedVcfFileHandler:
			record = line.strip().split('\t')
			if not (record[0], record[1], record[3], record[4]) in freqDic:
				print "Error: %s, %s, %s, %s not found in frequency file" %(record[0], record[1], record[3], record[4])
				return False'''
	
		

def dedupVcf(infile):
	infile = open(infile, 'r')
	matrix = []
	dup = False
	for line in infile:
		if not line.startswith('#'):
			#print line
			line = line.split('\t')
			if not line[0].startswith("chr"):
				line[0] = "chr" + line[0]
			length = max(len(line[3]), len(line[4])) - 1
			matrix.append([line[0], int(line[1]), int(line[1]) + length, line[3], line[4]])
	sortedMatrix = sorted(matrix, key=itemgetter(0,1,2))
	print sortedMatrix
	outMatrix = []
	previous = None
	running = []
	previousMax = None
	for variant in sortedMatrix:
		#print variant
		if previous == None:
			#first record
			#print "first"
			previous = variant
			running.append(variant)
		else:
			if variant[0] != previous[0]:
				#Not the same chromosome
				#print "Not the same chromosome"
				outMatrix.append(random.choice(running))
				running = []
				running.append(variant)
				previous = variant
				previousMax = variant[2]
			elif variant[1] > previousMax and len(running) == 1:
				#Not overlapping with previous, and previous not overlapping with anything else
				#print "solo"
				outMatrix.append(running[0])
				running = []
				running.append(variant)
				previous = variant
				previousMax = variant[2]
			elif variant[1] > previousMax and len(running) > 1:
				#Not overlapping with previous, but previous overlapping with others
				#print "gather previous"
				dup = True
				outMatrix.append(random.choice(running))
				running = []
				running.append(variant)
				previous = variant
				previousMax = variant[2]
			else:
				#Over lapping with previous
				#print "append"
				running.append(variant)
				previous = variant
				previousMax = max(variant[2], previousMax)
	#last record
	outMatrix.append(random.choice(running))
	if dup == True:
		print "Warning: Duplicated records exist in the supplied vcf file, randomly seleted one out of the duplicated records to simulate" 
	return outMatrix

def dedupFreq(infile):
	infile = open(infile, 'r')
	matrix = []
	dup = False
	for line in infile:
		if not line.startswith('#'):
			line = line.split('\t')
			if not line[0].startswith("chr"):
				line[0] = "chr" + line[0]
			length = max(len(line[3]), len(line[2])) - 1
			matrix.append([line[0], int(line[1]), int(line[1]) + length, line[2], line[3], line[4]])
	sortedMatrix = sorted(matrix, key=itemgetter(0,1,2))
	outMatrix = []
	previous = None
	running = []
	previousMax = None
	for variant in sortedMatrix:
		#print variant
		if previous == None:
			#first record
			#print "first"
			previous = variant
			running.append(variant)
		else:
			if variant[0] != previous[0]:
				#Not the same chromosome
				#print "Not the same chromosome"
				outMatrix.append(random.choice(running))
				running = []
				running.append(variant)
				previous = variant
				previousMax = variant[2]
			elif variant[1] > previousMax and len(running) == 1:
				#Not overlapping with previous, and previous not overlapping with anything else
				#print "solo"
				outMatrix.append(running[0])
				running = []
				running.append(variant)
				previous = variant
				previousMax = variant[2]
			elif variant[1] > previousMax and len(running) > 1:
				#Not overlapping with previous, but previous overlapping with others
				#print "gather previous"
				dup = True
				outMatrix.append(random.choice(running))
				running = []
				running.append(variant)
				previous = variant
				previousMax = variant[2]
			else:
				#Over lapping with previous
				#print "append"
				running.append(variant)
				previous = variant
				previousMax = max(variant[2], previousMax)
	if dup == True:
		print "Warning: Duplicated records exist in the supplied vcf file, randomly seleted one out of the duplicated records to simulate" 
	return outMatrix	

def editBam(bamFile, variantNumber, outfile):
	global dedupedMatrix
	#print dedupedMatrix
	selectedMatrixDic,selectedMatrixLimit  = prepVariants(dedupedMatrix, variantNumber)
	#print selectedMatrixDic
	#print selectedMatrixLimit
	variantDic = VariantDic(selectedMatrixDic,selectedMatrixLimit)
	bamFile = pysam.AlignmentFile(bamFile, "rb")
	global outfileHandler = pysam.AlignmentFile(outfile, "wb", template = bamFile)
	n = 0
	for read in bamFile.fetch():
		n += 1
		chr = read.reference_name
		if not chr.startswith("chr"):
			chr = "chr" + chr
		start = int(read.reference_start) + 1
		stop = int(read.reference_end)
		variants = variantDic.getVariants(chr, start, stop)
		if variants == None or variants == []:
			outfileHandler.write(read)
		else:
			editRead(read, variants)
	bamFile.close()
	outfile.close()
	
def editRead(read, variants):
	global outfileHandler
			
	
	

class VariantDic:
	def __init__(self, variantDic, variantLimit):
		self.variantDic = variantDic
		self.variantLimit = variantLimit
	def getVariants(self, chr, start, stop):
		variantDic = self.variantDic
		variantLimit = self.variantLimit
		if chr not in variantDic:
			return None
		lower = variantLimit[chr][0]
		higher = variantLimit[chr][1]
		if stop < lower or start > higher:
			return None
		variantInChr = variantDic[chr]
		start = int(start)
		stop = int(stop)
		index = 0
		overlapped = []
		#print "chr: %s, start: %d, stop: %d" %(chr, start, stop)
		readRegion = range(start, stop+1)
		while index < len(variantInChr):
			variantRegion = range(variantInChr[index][1], variantInChr[index][2] + 1)
			if bool(set(readRegion) & set(variantRegion)):
				#overlapping
				print "overlapping variant: %s, %s, %s" %(variantInChr[index][0], variantInChr[index][1], variantInChr[index][2])
				overlapped.append(variantInChr[index])
				#print "overlapping variant: %s, %s, %s" %(variantInChr[index][0], variantInChr[index][1], variantInChr[index][2])
			index += 1
		return overlapped
	'''def getOverlap(self, variants, start, stop, index):
		output = []
		output.append(variants[index])
		#scan lower
		
		front = index - 1
		while front >= 0:
			if variants[front][1] >= start and variants[front][2] <= stop:
				#overlap
				output.append(variants[front])
				front -= 1
		#scan higher
		back = index + 1
		while back <= len(variants) - 1:
			if variants[front][1] >= start and variants[front][2] <= stop:
				#overlap
				output.append(variants[front])
				front -= 1'''
		
	
	
	
def prepVariants(dedupedMatrix, variantNumber):
	global vcfFile
	if vcfFile != None:
		#assign random frequency
		for element in dedupedMatrix:
			freq = random.uniform(float(minFreq), float(maxFreq))
			element.append(freq)
	if variantNumber != "0":
		variantNumber = int(variantNumber)
		random.shuffle(dedupedMatrix)
		selectedMatrix = dedupedMatrix[0:variantNumber]
	else:
		selectedMatrix = dedupedMatrix
	selectedMatrixDic = {}
	for element in selectedMatrix:
		if element[0] not in selectedMatrixDic:
			selectedMatrixDic[element[0]] = []
			selectedMatrixDic[element[0]].append(element)
		else:
			selectedMatrixDic[element[0]].append(element)
	selectedMatrixLimit = {}
	for key in selectedMatrixDic:
		selectedMatrixDic[key] = sorted(selectedMatrixDic[key], key = itemgetter(1,2))
		lower = min(zip(*selectedMatrixDic[key])[1])
		higher = max(zip(*selectedMatrixDic[key])[2])
		selectedMatrixLimit[key] = (lower, higher)
	return selectedMatrixDic, selectedMatrixLimit


if __name__ == "__main__":
	if len(sys.argv) == 1:
		print "For usage: python editBam.py -h"
	else:
		parser = optparse.OptionParser()
		parser.add_option('-b', '--bamFile',
		help = 'input bamfile', metavar = 'FILE', default = None)
		parser.add_option('-o', '--outfile',
		help = 'output bam file', metavar = 'FILE', default = None)
		parser.add_option('-v', '--vcfFile',
		help = 'input vcf file', metavar = 'FILE', default = None)
		parser.add_option('-n', '--variantNumber',
		help = 'number of variants to insert, default: insert all', default = "0", metavar = 'NUMBER')
		parser.add_option('-l', '--minFreq',
		help = 'minimum allele frequency, default: 0.5(0.5%)', default = "0.5", metavar = 'NUMBER')
		parser.add_option('-L', '--maxFreq',
		help = 'maxinum allele frequency, default: 100(100%)', default = "100", metavar = 'NUMBER')
		parser.add_option('-f', '--freqFile',
		help = 'frequency file to use, one line per variant', default = None, metavar = 'FILE')
		options,args = parser.parse_args()
		bamFile = options.bamFile
		outfile = options.outfile
		vcfFile = options.vcfFile
		variantNumber = options.variantNumber
		minFreq = options.minFreq
		maxFreq = options.maxFreq
		freqFile = options.freqFile
		if(freqFile != None and (minFreq != "0.5" or maxFreq != "100")):
			print "Error: -f is not compatible with -l and -L"
			exit()
		if(bamFile == None or outfile == None):
			print "Error: -b, -o are required arguments"
			exit()
		if(freqFile != None and vcfFile != None):
			print "Error: -v and -f are mutually exclusive"
			exit()
		if(freqFile == None and vcfFile == None):
			print "Error: One of -v and -f needs to be specified"
			exit()
		if checkInput(bamFile, outfile, vcfFile, variantNumber, minFreq, maxFreq, freqFile):
			print "check passed"
			editBam(bamFile, variantNumber, outfile)
			
			
'''if __name__ == "__main__":
	infile = sys.argv[1]
	outfile = sys.argv[2]
	outfile = open(outfile, 'w')
	matrix = dedupVcf(infile)
	for variant in matrix:
		variantStr = map(str, variant)
		outfile.write('\t'.join(variantStr))
		outfile.write('\n')
	outfile.close()'''
		