import sys
import optparse
import re
from operator import itemgetter
from subprocess import check_output
from subprocess import call



def correct_base(ref, alt):
	ref_right = ref[1:]
	alt_right = alt[1:]
	base_corrected = "N"
	if len(ref_right) < len(alt_right):
		#insertion
		regex = ref_right + "$"
		m = re.search(regex, alt_right)
		if m:
			# in case of a simple insertion event, get the inserted bases
			insertion = alt_right[0:m.start()]
			ref = ref[0]
			alt = alt[0] + insertion
			base_corrected = "Y"
	elif len(ref_right) > len(alt_right):
		#deletion
		regex = alt_right + "$"
		m = re.search(regex, ref_right)
		if m:
			# in case of a simple deletion event, get the deleted bases
			deletion = ref_right[0:m.start()]
			ref = ref[0] + deletion
			alt = alt[0]
			base_corrected = "Y"
	else:
		#check if it is simple SNP case
		if ref_right == alt_right:
			ref = ref[0]
			alt = alt[0]
			base_corrected = "Y"
	return(ref, alt, base_corrected)
		

	
def cosmic_to_dic(cosmic):
	infile = open(cosmic, 'r')
	cosmic_dic = {}
	print "reading cosmic file"
	for line in infile:
		if not line.startswith('#'):
			line = line.split('\t')
			chr = line[0]
			if not chr.startswith('chr'):
				chr = "chr" + chr
			pos = line[1]
			id = line[2]
			ref = line[3]
			alt = line[4]
			info = line[7].split(';')
			info_dic = {}
			for record in info:
				record = record.split('=')
				try:
					info_dic[record[0]] = record[1]
				except IndexError:
					#no "=" in this information
					continue
			try:
				gene = info_dic["GENE"]
			except KeyError:
				gene = "."
			try:
				cDNA = info_dic["CDS"]
			except KeyError:
				cDNA = "."
			try:
				protein = info_dic["AA"]
			except KeyError:
				protein = "."
			cosmic_dic[id] = [chr, pos, ref, alt, cDNA, protein, gene]
	infile.close()
	print 'finished reading cosmic'
	return cosmic_dic
	
			
def effect(no_cosmic, outfile, report_no, barcode, outfile_dir):
	print 'getting effect info for novel mutations'
	tempFileName = outfile_dir + '/temp.vcf'
	temp = open(tempFileName, 'w')
	tempEffName = outfile_dir + '/temp.eff.vcf'
	record_dic = {}  #mutation record
	for record in no_cosmic:
		temp.write('\t'.join([record[0], record[1], ".", record[3], record[4], ".", "."]) + '\n')
		record_dic[(record[0], record[1], record[3], record[4])] = record  #mutation record
	temp.close()
	print 'running snpeff'
	tempEff = open(tempEffName, 'w')
	command = ['java', '-Xmx4g', '-jar', '/opt/software/snpEff-4.1b/snpEff.jar', 'GRCh37.75', '-nostats', tempFileName]
	call(command, stdout=tempEff)
	tempEff.close()
	print 'finished running snpeff'
	effect_file = open(tempEffName, 'r')
	for line in effect_file:
		if line.startswith('#'):
			continue
		else:
			line = line.split('\t')
			eff = line[7]
			effects = eff.split('=')[1].split(',')  #list of effects
			sig_dic = {"HIGH":1, "MODERATE":2, "LOW":3, "MODIFIER":4}
			eff_matrix = []  #each effect is a list in this matrix
			for effect in effects:
				effect = effect.split('|')
				effect[2] = sig_dic[effect[2]]  # convert effect significance to a number to be sorted later
				eff_matrix.append(effect)
			eff_matrix.sort(key=itemgetter(2))
			sig_effect = eff_matrix[0]  #most severe effect
			if sig_effect[9] == "":
				sig_effect[9] = "."
			if sig_effect[10] == "":
				sig_effect[10] = "."
			if sig_effect[3] == "":
				sig_effect[3] = "."
			mutation = record_dic[(line[0], line[1], line[3], line[4])]  #get full mutation record for this mutation
			output =  mutation + [sig_effect[3], sig_effect[9], sig_effect[10]] + [report_no, barcode] 
			outfile.write('\t'.join(output) + '\n')
	return(outfile)
	
	
	

#Obsolete function: The below function will assign annotations to the none cosmic mutations which is not necessary
'''def gene50excelLoad(infile, cosmic, outfile, report_no, barcode):
	infile = open(infile, 'r')
	cosmic_dic = cosmic_to_dic(cosmic)
	outfile_dir = re.findall('/home/scratch/.*/', outfile)[0]
	outfile = open(outfile, 'w')
	no_cosmic = []
	for line in infile:
		if line.startswith('Chrom'):
			#header
			continue
		else:
			line = line.split('\t')
			if line[4] == "Heterozygous" or line[4] == "Homozygous":
				#called
				if line[11] != "---":
					#has a cosmic ID, use the cosmic coordinates
					try:
						chr = cosmic_dic[line[11]][0]
						pos = cosmic_dic[line[11]][1]
						ref = cosmic_dic[line[11]][2]
						alt = cosmic_dic[line[11]][3]
						cDNA = cosmic_dic[line[11]][4]   # cDNA change
						protein = cosmic_dic[line[11]][5]   #codon change
						gene = cosmic_dic[line[11]][6]     #gene
						depth = line[18]      #read depth
						altCount = line[24]     #alt allele count
						frequency = line[6]    #alt allele frequency
						genotype = line[4]
						cosmic = line[11]
						outfile.write('\t'.join([chr, pos, cosmic, ref, alt, genotype, depth, altCount, frequency, gene, cDNA, protein,report_no, barcode]) + '\n')
					except KeyError:
						#cosmic ID not found in cosmic reference file, use the vcf coordinates
						chr = line[0]
						pos = line[14]
						ref = line[15]
						alt = line[16]
						depth = line[18]
						altCount = line[24]
						frequency = line[6]
						genotype = line[4]
						cosmic = line[11]
						#correct base if necessary
						(ref,alt, base_corrected) = correct_base(ref,alt)
						no_cosmic.append([chr, pos, cosmic, ref, alt, genotype, depth, altCount, frequency])
					
				else:
					#no cosmic ID, use the vcf coordinates
					chr = line[0]
					pos = line[14]
					ref = line[15]
					alt = line[16]
					depth = line[18]
					altCount = line[24]
					frequency = line[6]
					genotype = line[4]
					cosmic = line[11]
					#correct base if necessary
					(ref,alt, base_corrected) = correct_base(ref,alt)
					no_cosmic.append([chr, pos, cosmic, ref, alt, genotype, depth, altCount, frequency])
	#get effect of all mutations without cosmic ID
	outfile = effect(no_cosmic, outfile, report_no, barcode, outfile_dir)
	infile.close()
	outfile.close()'''
	

def gene50excelLoad(infile, outfile, report_no, barcode, callerID):
	infile = open(infile, 'r')
	outfile = open(outfile, 'w')
	for line in infile:
		if line.startswith('Chrom'):
			#header
			continue
		else:
			line = line.split('\t')
			if line[4] == "Heterozygous" or line[4] == "Homozygous":
				#called
				chr = line[0]
				pos = line[14]
				ref = line[15]
				alt = line[16]
				depth = line[18]
				altCount = line[24]
				frequency = line[6]
				genotype = line[4]
				cosmic = line[11]
				#correct base if necessary
				(ref,alt, base_corrected) = correct_base(ref,alt)
				outfile.write('\t'.join([chr, pos, cosmic, ref, alt, genotype, depth, altCount, frequency, report_no, barcode, callerID]) + '\n')
			if line[4] == "No Call" and re.search("Minimum coverage", line[5]):
				#nocall, check no call reason, report if real minimum coverage
				chr = line[0]
				pos = line[14]
				ref = line[15]
				alt = line[16]
				depth = line[18]
				altCount = line[24]
				frequency = line[6]
				genotype = line[4]
				cosmic = line[11]
				#correct base if necessary
				(ref,alt, base_corrected) = correct_base(ref,alt)
				outfile.write('\t'.join([chr, pos, cosmic, ref, alt, genotype, depth, altCount, frequency, report_no, barcode, callerID]) + '\n')
				
				
	infile.close()
	outfile.close()
	
try:
	
		
	if sys.argv[1] == 'gene50excelLoad':
		parser = optparse.OptionParser()
		parser.add_option('-I', '--infile',
		help = 'input excel file')
		parser.add_option('-o', '--outfile',
		help = 'output file')
		parser.add_option('-r', '--report_no',
		help = 'ion torrent report number')
		parser.add_option('-b', '--barcode',
		help = 'sequencing barcode (three digits)')
		parser.add_option('-c', '--callerID', default = 'None',
		help = 'variant caller ID')
		options,args = parser.parse_args()
		infile = options.infile
		outfile = options.outfile
		report_no = options.report_no
		barcode = options.barcode
		callerID = options.callerID
		gene50excelLoad(infile, outfile, report_no, barcode, callerID)

except ValueError:
	print 'Available functions:'
	print ' gene50excelLoad(): take excel output of ion torrent and output a file for dbloading'
	
