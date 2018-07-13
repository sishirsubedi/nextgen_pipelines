import sys
import optparse
import re
import commonTool
from reportlab.pdfgen import canvas
from reportlab.lib.units import inch
from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import Paragraph



def get_data(infile):
	dataPos = []
	dataNeg = []
	infile = open(infile, 'r')
	for line in infile:
		if line.startswith('#') or line.startswith('@'):
			continue
		else:
			data = line.split('\t')
			if data[7] == "0/1" or data[7] == "1/0" or data[7] == "0|1" or data[7] == "1|0":
				data[7] = "heterozygous"
			if data[7] == "1/1" or data[7] == "1|1":
				data[7] = "homozygous"
			if re.search('pathogenic', line, re.IGNORECASE):
				dataPos.append(data)
			else:
				dataNeg.append(data)
	infile.close()		
	return(dataPos, dataNeg)
	
def get_interpretation(data):
	data[23] = data[23].strip()
	if data[23] != "\N":
		interpretation = data[23]
	else:
		sig = data[17].split('|')
		if "Pathogenic" in sig:
			interpretation = "Pathogenic"
		elif "Likely_Pathogenic" in sig:
			interpretation = "Likely Pathogenic"
		elif "Uncertain_significance" in sig:
			interpretation = "Uncertain significance"
		else:
			interpretation = "other"
	return interpretation

def get_significance(data):
	sig = data[17].strip('|')
	sig = sig.split('|')
	disease = data[19].split('|')
	n = 0
	significance = "This mutation is "
	while n < len(sig) and n < len(disease):
		if re.search('pathogenic', sig[n], re.IGNORECASE):
			significance = significance + " %s for %s; " % (sig[n], disease[n]) 
		n += 1
	return(significance)

def check_positive(infile):
	infile = open(infile, 'r')
	positive_check = "n"
	for line in infile:
		if re.search('pathogenic', line, re.IGNORECASE):
			positive_check = "y"
			break
	return(positive_check)

def plot_id(c, id):
	c.setFillColor(colors.lightgrey)
	c.setStrokeColor(colors.white)
	c.rect(1*inch, 10*inch, 6*inch, 0.5*inch, fill = 1)
	c.setFillColor(colors.grey)
	c.drawString(1.1*inch, 10.2*inch, "ID")
	return c	
	
def draw_header_positive(c):
	c.setFillColorRGB(1,0,0)
	c.setStrokeColor(colors.white)
	c.rect(1*inch, 9.4*inch, 6*inch, 0.5*inch, fill = 1)
	c.setFillColorRGB(1,1,1)
	c.setFont("Helvetica",15)
	c.drawString(1*inch, 9.6*inch, "Result:POSITIVE -- PATHOGENIC MUTATIONS IDENTIFIED")
	return(c)
	
def draw_header_negative(c):
	c.setFillColorRGB(0.8,0.8,0.8)
	c.setStrokeColor(colors.white)
	c.rect(1*inch, 9.4*inch, 6*inch, 0.5*inch, fill = 1)
	c.setFillColorRGB(0,0,0)
	c.setFont("Helvetica",15)
	c.drawString(1*inch, 9.6*inch, "Result:NEGATIVE -- NO PATHOGENIC MUTATIONS IDENTIFIED")
	return(c)

def plot_detail(c, n, mutation, significance, property):
	if n % 2 == 0:
		#decide if it is the first or second record plotted
		ymove = 4
	else:
		ymove = 0
	#print "mutation is %s" % mutation
	#print "significance is %s" % significance
	#print "property is %s" % property
	styleSheet = getSampleStyleSheet()
	style = styleSheet['BodyText']
	detail = "mutation %s results in %s <br />\n <br />\n  %s" % (mutation, property, significance)
	p = Paragraph(detail, style)
	w,h = p.wrap(6.5*inch, 3*inch)
	p.drawOn(c, 1*inch, (7-ymove)*inch - h/2) 
	return(c)	

def plot_record_positive(data,n,c):
	#print data
	if n % 2 == 0:
		#decide if it is the first or second record plotted
		ymove = 4
	else:
		ymove = 0
	#record header boxes
	c.setFillColor(colors.lavenderblush)
	c.setStrokeColor(colors.white)
	c.rect(1*inch, (9-ymove)*inch, 1*inch, 0.3*inch, fill = 1)
	c.rect(2.1*inch, (9-ymove)*inch, 1.9*inch, 0.3*inch, fill = 1)
	c.rect(4.1*inch, (9-ymove)*inch, 2.9*inch, 0.3*inch, fill = 1)
	c.rect(1*inch, (7.5-ymove)*inch, 6*inch, 0.3*inch, fill = 1)
	#record header texts
	c.setFillColor(colors.gray)
	c.setFont("Helvetica",10)
	c.drawString(1.1*inch, (9.1-ymove)*inch, "GENE")
	c.drawString(2.2*inch, (9.1-ymove)*inch, "MUTATION")
	c.drawString(4.2*inch, (9.1-ymove)*inch, "INTERPRETATION")
	c.drawString(1.1*inch, (7.6-ymove)*inch, "Details")
	#record contenet texts
	#plot gene, mutation
	c.setFillColor(colors.black)
	c.setFont("Helvetica",10)
	c.drawString(1.1*inch, (8.6-ymove)*inch, data[16])
	c.drawString(2.2*inch, (8.8-ymove)*inch, data[14]+"(" +data[15] + ")")
	c.drawString(2.2*inch, (8.3-ymove)*inch, data[7])
	#plot interpretation
	interpretation = get_interpretation(data)
	c.drawString(4.2*inch, (8.6-ymove)*inch, interpretation)
	#plot details
	significance = get_significance(data)
	property = commonTool.hgvs_protein_parser(data[15])
	c = plot_detail(c, n, data[14], significance, property)
	return(c)
	

def plot_uncertain_header(c,n):

	#plot additional finding header
	#print "plotting uncertain header"
	if n%2 == 0:
		ymove = 4
	else:
		ymove = 0
	c.setFillColor(colors.gray)
	c.rect(1*inch, (9.4-ymove)*inch, 6*inch, 0.3*inch, fill = 1)
	c.setFillColor(colors.black)
	c.drawString(1.1*inch, (9.5-ymove)*inch, "VARIANTS OF UNCERTAIN SIGNIFICANCE")
	ymove = ymove + 0.4*inch
	return(c)
	
def plot_record_negative(data,n,c,m):
	#print data
	if n % 2 == 0:
		#decide if it is the first or second record plotted
		ymove = 4
	else:
		ymove = 0
	
	#plot uncertain variants
	#header boxes
	
	c.setFillColor(colors.gray)
	c.setStrokeColor(colors.white)
	c.rect(1*inch, (9-ymove)*inch, 1*inch, 0.3*inch, fill = 1)
	c.rect(2.1*inch, (9-ymove)*inch, 1.9*inch, 0.3*inch, fill = 1)
	c.rect(4.1*inch, (9-ymove)*inch, 2.9*inch, 0.3*inch, fill = 1)
	c.rect(1*inch, (7.5-ymove)*inch, 6*inch, 0.3*inch, fill = 1)
	#print "print uncertain records"
	#header texts
	c.setFillColor(colors.black)
	c.setFont("Helvetica",10)
	c.drawString(1.1*inch, (9.1-ymove)*inch, "GENE")
	c.drawString(2.2*inch, (9.1-ymove)*inch, "MUTATION")
	c.drawString(4.2*inch, (9.1-ymove)*inch, "INTERPRETATION")
	c.drawString(1.1*inch, (7.6-ymove)*inch, "Details")
	#record contenet texts
	#plot gene, mutation
	c.setFillColor(colors.black)
	c.setFont("Helvetica",10)
	c.drawString(1.1*inch, (8.6-ymove)*inch, data[16])
	c.drawString(2.2*inch, (8.8-ymove)*inch, data[14]+"(" +data[15] + ")")
	c.drawString(2.2*inch, (8.3-ymove)*inch, data[7])
	#plot interpretation
	interpretation = get_interpretation(data)
	c.drawString(4.2*inch, (8.6-ymove)*inch, interpretation)
	#plot details
	significance = "The significance of this mutation is not determined"
	property = commonTool.hgvs_protein_parser(data[15])
	c = plot_detail(c, n, data[14], significance, property)
	return(c)

def plot_report(infile, outFile):
	c = canvas.Canvas(outFile, pagesize=letter)
	c = plot_id(c, "ID")
	positive_check = check_positive(infile)
	if positive_check == "y":
		c = draw_header_positive(c)
	if positive_check == "n":
		c = draw_header_negative(c)
	n = 1
	(dataPos, dataNeg) = get_data(infile)
	for data in dataPos:
		c = plot_record_positive(data, n, c)
		if n % 2 == 0:
			#two records/page
			c.showPage()
			c = plot_id(c, "ID")
		n += 1
	m = 0
	for data in dataNeg:
		#print data
		#print n
		if m == 0:
			c = plot_uncertain_header(c,n)
		c = plot_record_negative(data, n, c,m)
		if n % 2 == 0:
			c.showPage()
			c = plot_id(c, "ID")
		n += 1
		m += 1
	c.save()
	
try:
	if sys.argv[1] == 'plot_report':
		parser = optparse.OptionParser()
		parser.add_option('-I', '--infile',
		help = 'input file')
		parser.add_option('-o', '--outfile',
		help = 'output file')
		options,args = parser.parse_args()
		infile = options.infile
		outfile = options.outfile
		plot_report(infile,outfile)
except ValueError:
		print 'Functions:'
		print ' plot_report(): Take input breast hereditary clinical text report and generate pdf report'