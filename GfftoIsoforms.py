#!/usr/local/bin/env python
#Version 1.1 David Clarke, 12.10.2014.
#Contact David Clarke, daveclarke.bioinfo@gmail.com

#--------------------------------
#GfftoIsoforms.py
#--------------------------------

#This script is for identifying overlapping CDS sequences in a gff3 file
#The input is a gff3 file (-i) and output name without the file extension (-o)
#There is optional input of a non standard 'Name' field (-n)
#The output is a txt file showing entries that form overlapping clusters and a modified gff3 file
#Output options are --rename, renameonly and --clusterall
#The default is to create new entries in the gff3 for CDS lines that form clusters
#default entries are 'ClusterName','Cluster' and 'ClusterNum'
#If --rename is specified then the 'ClusterName' data will instead be written into 'Name'
# --renameonly will do the same as above but also not include the 'Cluster' and 'ClusterNum' entries
#If --clusterall is specified then all CDS lines will be assigned to clusters (clusters of one)

###Import modules
import time
import math
import sys
import argparse;
from operator import itemgetter

###Function definitions.
def READLINES(lines,idfield):
	ExonList = []
	NotEXON = []
	hash = "#"
	for line in lines:
		if line[0] != hash:
			rline = line.rstrip()
			data = rline.split('\t')
	
			if len(data) == 9:
				scaff = data[0]
				source = data[1]
				feature = data[2]
				exonstart = int(data[3])
				exonend = int(data[4])
				quality = data[5]
				exondir = data[6]
				frame = data[7]
				exoninfo = data[8]
	
				Name = 'none'
				Parent = 'none'
				ID = 'none'
				INFO = {}
				try:
					if exoninfo[-1] == ';':
						exoninfo = exoninfo[:-1]
					Notes = exoninfo.split(';')
					for info in Notes:
						parts = info.split('=')
						INFO[parts[0]] = parts[1]
					if 'Name' in INFO and idfield == "":
						Name = INFO['Name']
					elif idfield != "" and idfield in INFO:
						Name = INFO[idfield]
					if 'Parent' in INFO:
						Parent = INFO['Parent']
					if 'ID' in INFO:
						ID = INFO['ID']

				except Exception, e:
					print "Error with some gff data" + str(e)
			
				if feature == "CDS":
					ExonList.append(NewGffExon(scaff,source,feature,exonstart,exonend,quality,exondir,frame,exoninfo,Name,Parent,ID,INFO))
					if Name == 'none':
						sys.exit("CDS entry did not have a 'Name' assigned, check gff3 for a uniqe naming field in the info\n"+rline+"\nUse '-n NameID' to indicate the naming field\nExiting")
				else:
					NotEXON.append(NewGffExon(scaff,source,feature,exonstart,exonend,quality,exondir,frame,exoninfo,Name,Parent,ID,INFO))

	return ExonList

def WRITELINES(lines,idfield,rename,GENEDIC,CLUSTERDIC,outgff):
#	ExonList = []
#	NotEXON = []
	intpos = len(str(len(CLUSTERDIC)))
	hash = "#"
	for line in lines:
		if line[0] == hash:
#			print "Not a data line"
			outgff.write(line)
		else:	
			rline = line.rstrip()
			data = rline.split('\t')
	
			if len(data) != 9:
#				print "Data line does not have 9 elements"
				outgff.write(line)
			else:
				scaff = data[0]
				source = data[1]
				feature = data[2]
				exonstart = int(data[3])
				exonend = int(data[4])
				quality = data[5]
				exondir = data[6]
				frame = data[7]
				exoninfo = data[8]
	
				Name = 'none'
				Parent = 'none'
				ID = 'none'
				INFO = {}
				try:
					
					if exoninfo[-1] == ';':
						exoninfo = exoninfo[:-1]
					Notes = exoninfo.split(';')
					for info in Notes:
						parts = info.split('=')
						INFO[parts[0]] = parts[1]
					if 'Name' in INFO and idfield == "":
						Name = INFO['Name']
						idfield = 'Name'
					elif idfield != "" and idfield in INFO:
						Name = INFO[idfield]
					if 'Parent' in INFO:
						Parent = INFO['Parent']
					if 'ID' in INFO:
						ID = INFO['ID']
				except Exception, e:
					print "Error with some gff data" + str(e)
			
				if feature == "CDS":
					Thisexon = NewGffExon(scaff,source,feature,exonstart,exonend,quality,exondir,frame,exoninfo,Name,Parent,ID,INFO)
					if Thisexon.HitIdent in GENEDIC:
						INFO['Cluster'] = "CLUST" + "0"*(intpos-len(str(GENEDIC[Thisexon.HitIdent])))+str(GENEDIC[Thisexon.HitIdent])
						isoform = CLUSTERDIC[GENEDIC[Thisexon.HitIdent]].index(Thisexon.HitIdent) + 1
						isopos = len(str(len(CLUSTERDIC[GENEDIC[Thisexon.HitIdent]])))
						INFO['ClusterNum'] = "0"*(isopos-len(str(isoform)))+str(isoform) + "_" + str(len(CLUSTERDIC[GENEDIC[Thisexon.HitIdent]]))
						if rename == 0:
							INFO['ClusterName'] = Name + "_" + INFO['Cluster'] + "_" + INFO['ClusterNum']
						else:
							INFO[idfield] = Name + "_" + INFO['Cluster'] + "_" + INFO['ClusterNum']
						newinfo = []
						if "ID" in INFO:
							newinfo.append("ID" + "=" + INFO["ID"])						
						if 'Name' in INFO:
							newinfo.append('Name' + "=" + INFO['Name'])						
						elif idfield in INFO:
							newinfo.append(idfield + "=" + INFO[idfield])						
						if 'Parent' in INFO:
							newinfo.append('Parent' + "=" + INFO['Parent'])						
						if 'ClusterName' in INFO:
							newinfo.append('ClusterName' + "=" + INFO['ClusterName'])						
						for meta in INFO:
							noextra = []
							if rename == 2:
								noextra = ['Cluster','ClusterNum']
							if meta not in ["ID",'Name','Parent','ClusterName',idfield] + noextra:
								newinfo.append(meta + "=" + INFO[meta])						
						outgff.write(Thisexon.noinfogff() + ";".join(newinfo) + "\n")
					else:
						outgff.write(line)
				else:
	#				NotEXON.append(NewGffExon(scaff,source,feature,exonstart,exonend,quality,exondir,frame,exoninfo,Name,Parent,ID,INFO))
					outgff.write(line)

#	return ExonList

def SORTEXONS(scaffold,ExonList):
	CSLF = []
	CSLR = []
	print "\n" + "Finding exons in " + scaffold
	for y in ExonList:
		if y.scaff == scaffold:
			if y.exondir == "+":# or "-":
				CSLF.append(y)
			elif y.exondir == "-":# or "+":
				CSLR.append(y)
			else:
				print "Read dir error\n" + str(y)

	RSCSLF = sorted(CSLF, key=lambda gffexon: gffexon.exonend, reverse=True)
	RSCSLR = sorted(CSLR, key=lambda gffexon: gffexon.exonstart, reverse=False)
			
	FSCSLF = sorted(RSCSLF, key=lambda gffexon: gffexon.exonstart, reverse=False)
	FSCSLR = sorted(RSCSLR, key=lambda gffexon: gffexon.exonend, reverse=True)
		
	return FSCSLF,FSCSLR
						
def CLUSTSCAF(HExons,dir):

	blank = (" "," ")
	print "Number of " + dir + " exons : " + str(len(HExons))

	Hnames = []
	for x in HExons:
		if x.HitIdent not in Hnames:
			Hnames.append(x.HitIdent)
			
	PPOSITION = []
	posrange = [0,0]
	for x in HExons:
		lp = len(PPOSITION)
		exonstart = int(x.exonstart)-1
		exonend = int(x.exonend)
		if posrange[1] == 0:
			posrange[0] = exonstart
			posrange[1] = exonend
			PPOSITION.append([])
			lp = len(PPOSITION)
			if x in HExons:
				PPOSITION[lp-1].append(x)
			else:
				print "Exon not matched " + str(x)
		elif exonstart < posrange[1] and dir == "+":
			if x in HExons:
				PPOSITION[lp-1].append(x)
			else:
				print "Exon not matched " + str(x)
			if posrange[1] < exonend:
				posrange[1] = exonend
		elif exonend > posrange[0] and dir == "-":
			if x in HExons:
				PPOSITION[lp-1].append(x)
			else:
				print "Exon not matched " + str(x)
			if posrange[0] > exonstart:
				posrange[0] = exonstart
		else: #Start new position
			posrange = [exonstart,exonend]
			PPOSITION.append([])
			lp = len(PPOSITION)
			if x in HExons:
				PPOSITION[lp-1].append(x)
			else:
				print "Exon not matched " + str(x)
	
	PSEQS = {}
	n = 0
	for x in PPOSITION:
		n += 1
#		print "Position " + str(n)
		for y in x:
#			print "\t" + str(y)
			if str(y.HitIdent) in PSEQS:
				tsl = PSEQS[str(y.HitIdent)]
				tsl.append(n)
				PSEQS[str(y.HitIdent)] = tsl
			else:
				PSEQS[str(y.HitIdent)] = [n]

#	for x in PSEQS:
#		print str(x) + " " + str(PSEQS[x])

		
	CLUST = {}	
	for n in Hnames:
		Overlap = []
		PeptideExons = []
#		print n
		for e in HExons:
			if e.HitIdent == n:
				PeptideExons.append(e)	
		for p in PSEQS[n]:
#			print str(p)
			for t in PPOSITION[p-1]:
				if t.HitIdent not in Overlap:
					Overlap.append(t.HitIdent)
#					print t.HitIdent
		CLUST[n] = Overlap
		
	clusted = []
	clusters = []
	for n in Hnames:
		cluster = []
		if n not in clusted:
			cluster	= CHECKDIC(n,CLUST,[])
			for c in cluster:
				clusted.append(c)
			clusters.append(cluster)
			
		
	return clusters,CLUST
		
def CHECKDIC(n,CLUST,cluster):
	cluster.append(n)
	for x in CLUST[n]:
		if x not in cluster:
			cluster = CHECKDIC(x,CLUST,cluster)
	return cluster
	
class NewGffExon:
	scaff = ""
	source = ""
	feature = ""
	exonstart = -1
	exonend = -1
	quality = "."
	exondir = ""
	frame = "."
	exoninfo = ""
	Name = ""
	Parent = ""
	ID = ""
	INFO = {}
	HitIdent = ""
	def __init__(self,scaff,source,feature,exonstart,exonend,quality,exondir,frame,exoninfo,Name,Parent,ID,INFO):
		self.scaff = scaff
		self.exondir = exondir
		self.source = source
		self.feature = feature
		self.exonstart = exonstart
		self.exonend = exonend
		self.quality = quality
		self.exoninfo = exoninfo
		self.Name = Name
		self.Parent = Parent
		self.ID = ID
		self.INFO = INFO
		self.frame = frame
		if ID not in ['none',""]:
			self.HitIdent = Name + "-" + str(ID)
		else:
			self.HitIdent = Name
#	def __repr__(self):
#		return 'GffExon %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s' % (str(self.scaff),str(self.exondir),str(self.exonstart),str(self.exonend),str(self.HitName),str(self.HitNum),str(self.peplen),str(self.pepstart),str(self.pepend),str(self.RawScore),str(self.frame),str(self.ExonPlace),str(self.ExonTotal),str(self.ExonAAS),str(self.ExonAAE))
	
	def printgff(self):
		return '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (str(self.scaff),str(self.source),str(self.feature),str(self.exonstart),str(self.exonend),str(self.quality),str(self.exondir),str(self.frame),str(self.exoninfo))

	def noinfogff(self):
		return '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (str(self.scaff),str(self.source),str(self.feature),str(self.exonstart),str(self.exonend),str(self.quality),str(self.exondir),str(self.frame))
		
#	def __eq__(self, other): 
#		return self.__dict__ == other.__dict__

def main(ingff, outfile, idfield, rename, clusterall):	

	localtime = time.localtime(time.time())
	year = localtime[0]
	month = localtime[1]
	day = localtime[2]
	date =(str(day) + "." + str(month) + "." + str(year))
	print "day : " + str(day) + ", month : " + str(month) + ", year : " + str(year)

	if outfile + ".gff3" == ingff:
		print "Output file will override input file"
		sys.exit("Exiting")

	match = open(outfile + ".txt",'w')
	outgff = open(outfile + ".gff3",'w')
	#mtable = open(outfile + "2",'w')
	#ctable = open(outfile + "3",'w')

	try:
		cdslines = []
		cdsfile = open(ingff,'r')
		cdslines = cdsfile.readlines()
		cdsfile.close()
	except:
		print "problem reading gff files"
		
	PepExons = []
		
	try:
		PepExons = READLINES(cdslines,idfield)
	except Exception, e:
		print "problem reading gff lines" + str(e)
		
	try:
		Scaffolds = set()
		for x in PepExons:
			Scaffolds.add(x.scaff)
	except:
		print "error finding Scaffolds"

	CLUSTERDIC = {}
	GENEDIC = {}

	n = 0
	for x in sorted(Scaffolds):
		PEXONS = SORTEXONS(x,PepExons)
		PEF = PEXONS[0]
		PER = PEXONS[1]
		
		ScaffSeq = ""
		
		CSF = CLUSTSCAF(PEF,"+")
		CSR = CLUSTSCAF(PER,"-")
		
		for CS in [CSF,CSR]:
			for cluster in CS[0]:
				if len(cluster) > 1 or clusterall == 1:
					n += 1
					CLUSTERDIC[n] = cluster
					for gene in cluster:
						GENEDIC[gene] = n

	if len(CLUSTERDIC) > 0:
		intpos = len(str(len(CLUSTERDIC)))
		for n in CLUSTERDIC:
			cluster = CLUSTERDIC[n]
			if len(cluster) > 1:
				clust = "CLUST" + "0"*(intpos-len(str(n)))+str(n)
				match.write(str(clust) + "\t" + ",".join(cluster) + "\n")
	else:
		match.write("No Clusters found\n")
		
	try:
		WRITELINES(cdslines,idfield,rename,GENEDIC,CLUSTERDIC,outgff)
	except Exception, e:
		print "problem writing gff lines; " + str(e)
			

	outgff.close()
	match.close()
	#mtable.close()
	#ctable.close()

if __name__== '__main__':
	###Argument handling.
	arg_parser = argparse.ArgumentParser(description='');
	arg_parser.add_argument("-i","--input", default=None, help="Gff3 input file to be clustered.");
	arg_parser.add_argument("-o","--outfile", default=None, help="Desired name of output files (two files created .gff3 and .txt)");
	arg_parser.add_argument("-n","--idfield", default='Name', help="Field in gff3 notes column to use as unique identifier for CDS annotations, defaults to 'Name'.");
	arg_parser.add_argument("--rename", action='store_true', default=False, help="If specified then the 'ClusterName' data will instead be written into 'Name' field.");
	arg_parser.add_argument("--renameonly", action='store_true', default=False, help="Same as 'rename option' but will not also include 'Cluster' and 'ClusterNum' as additional note entries.");
	arg_parser.add_argument("--clusterall", action='store_true', default=False, help="If specified then all CDS lines will be assigned to clusters (clusters of one).");
	#arg_parser.add_argument("-b","--bool", action='store_true', default=False, help="");
	#arg_parser.add_argument("-v", '--verbose', action="count", help="Counts instaces of -v, -vv = very detailed");
	
	#Parse args and assign to var 'args'
	args = arg_parser.parse_args();

	###Variable Definitions
	#verbose=args.verbose;
	ingff=args.input;
	outfile=args.outfile;
	idfield=args.idfield;
	rename=args.rename;
	if rename:
		rename=1;
	renameonly=args.renameonly;
	if renameonly:
		rename=2;
	clusterall=args.clusterall;
	if clusterall:
		clusterall=1;

	#Feed to main()
	main(ingff, outfile, idfield, rename, clusterall);
