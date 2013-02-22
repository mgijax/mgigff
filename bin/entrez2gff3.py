#
# entrez2gff.py
#
# EntrezGene gene model to GFF3 format converter.
#
# ENTREZGENE FORMAT ----------------------------------------------
#
# EntrezGene format is a TAB-delimited text format with the
# following columns:
#	0: tax_id
#	1: chromosome
#	2: chr_start
#	3: chr_stop
#	4: chr_orient
#	5: contig
#	6: ctg_start
#	7: ctg_stop
#	8: ctg_orient
#	9: feature_name
#	10: feature_id
#	11: feature_type
#	12: group_label
#	13: transcript
#	14: evidence_code
#
# Comment lines begin with a HASH (#) character.
#
# The first line of the file is a comment giving the column labels.
#
# The following feature_types are used:
#	GENE, PSEUDO, RNA, UTR, CDS
#
# All features of a gene are tied to that gene by having the same
# feature_id. The parent/child relationships of subfeatures is
# implied by the order in which they occur in the file. 
# A feature of type UTR or CDS belongs to the most recent RNA.
#
# CONVERSION TO GFF3 ---------------------------------------------
#
# In addition to being GFF3, also make it look like the output of
# the gtf2gff (Perl) converter. I.e., try to follow the same name
# conventions, etc.
#
# Map feature types to SO:
#	GENE	-> gene
#	PSEUDO	-> pseudogene
#	RNA	-> mRNA or transcript, depending on if CDS if present or not
#	CDS	-> CDS
#	UTR	-> five_prime_UTR or three_prime_UTR, depending
#		   on position relative to CDS
#
# Create exons for each transcript out of 3' and 5' utrs + cds.
# 
# ----------------------------------------------------------------

import sys
import os
import re
import gff3
import logging
import optparse

# character constants
HASH	= '#'
TAB	= '\t'
NL	= '\n'

# GFF3 empty/undefined value
EMPTY	= '.'

# constants defining the input fields
iTAX	= 0 # tax_id
iCHR	= 1 # chromosome
iSTART	= 2 # chr_start
iEND	= 3 # chr_stop
iSTRAND	= 4 # chr_orient
iCONTIG	= 5 # contig
iCSTART	= 6 # ctg_start
iCEND	= 7 # ctg_stop
iCSTRAND= 8 # ctg_orient
iFNAME	= 9 # feature_name
iFID	= 10 # feature_id
iFTYPE	= 11 # feature_type
iGROUP	= 12 # group_label
iTRANS	= 13 # transcript
iECODE	= 14 # evidence_code

# NCBI file contains annotations from multiple assemblies. We
# just want the standard black-6 genes.
ASSEMBLY = 'MGSCv37-C57BL/6J'

class Entrez2GFF3(object):
    def __init__(self):
	self.iFileName = None
	self.iFile = None
	self.oFileName = None
	self.oFile = None

	self.currLineNum = None
	self.currLine = None
	self.currTokens = None
	self.currOutLineNum = None

	self.currModelID = None
	self.currModel = []
	self.idCounters = {}

    #-----------------------------------------------------
    def parseArgs(self, argv):
	self.optParser = optparse.OptionParser()
        self.optParser.add_option(
            "-a",
            dest="assembly",
            default=None,
            metavar="ASSEMBLY",
            help="The NCBI data file includes multiple assemblies. " +
		"This specifies which assembly to select. " +
		"Example black 6 build 37: 'MGSCv37-C57BL/6J' " +
		"Example black 6 build 38: 'GRCm38-C57BL/6J' " +
		""
	    )
	(self.options,self.posArgs) = self.optParser.parse_args(argv)
	if self.options.assembly is None:
	    self.optParser.error("No assembly specified.")
	if len(self.posArgs) == 2:
	    self.iFileName = self.posArgs[1]
	else:
	    self.iFileName = '-'

    #-----------------------------------------------------
    def openFiles(self):
	if self.iFileName == '-':
	    self.iFile = sys.stdin
	else:
	    self.iFile = open(self.iFileName, 'r')
	self.currLineNum = 0
	self.currOutLineNum = 0
	self.oFile = sys.stdout

    #-----------------------------------------------------
    def getIDCounter(self, key):
        c = self.idCounters.get(key,1)
	self.idCounters[key] = c+1
	return c

    #-----------------------------------------------------
    def resetIDCounters(self):
	for k in self.idCounters.keys():
	    self.idCounters[k] = 1

    #-----------------------------------------------------
    def convertLine(self, ts):
	chrom	= ts[iCHR]
	start	= ts[iSTART]
	end	= ts[iEND]
	strand	= ts[iSTRAND]
	fname	= ts[iFNAME]
	fid	= ts[iFID]
	ftype	= ts[iFTYPE]
	trans	= ts[iTRANS]

	return gff3.Feature([
	    chrom,
	    ftype,
	    ftype,
	    start,
	    end,
	    EMPTY,
	    strand,
	    EMPTY,
	    { 'original' : ts }
	    ])

    #-----------------------------------------------------
    # This is where all the real work occurs...
    # self.currModel contains a list of features for one gene model.
    #
    def massageModel(self):
	#
	# local function to be used in sorting a gene model's featues
	def fkey( f ):
	    t = f.type
	    if t =="GENE":
		return (0,"",0,0)
	    elif t =="PSEUDO":
		return (1,"",0,0)
	    else:
		trans_id = f.attributes['original'][13]
		if t == "RNA":
		    x = 0
		else:
		    x = 1
		return (2, trans_id, x, f.start)
	#
	self.currModel.sort(None,fkey)
	gene = self.currModel[0]

	# Check for multiple genes (models) w/ same ID.
	# Report them and exclude from further processing.
	if len(self.currModel) > 1 and self.currModel[1].type == "GENE":
	    logging.warn("Duplicate GeneID. Skipping: %s\n%s\n"%(self.currModelID,str(gene)))
	    self.currModel = []
	    return
	    

	#extract original data line and reset the attributes
	orig = gene.attributes['original']
	gene.attributes = { 'ID' : self.currModelID }

	# To be consistent with gtf2gff3's output... put original type in
	# source col. The type col is always "gene".
	gene.source = gene.type 
	gene.type = "gene"

	strand   = gene.strand

	massaged = [gene]

	# loop through the current gene model. For each
	# transcript, record the start/end indexes bracketing
	# the UTRs and CDSs. 
	transcripts = []
	for i,x in enumerate(self.currModel):
	    if x.type == "RNA":
		# Each entry looks like this: [startindex, stopindex, hasCDS, hasRNA]
		#
		transcripts.append([i,i,False,True])
	    elif x.type in ["UTR","CDS"]:
		if len(transcripts) == 0:
		    # There are cases (Ig's) of CDSs w/o any RNA record.
		    transcripts.append([i,i,False,False])
		transcripts[-1][1] = i
		if x.type == "CDS":
		    transcripts[-1][2] = True
	    elif x.type == "PSEUDO":
		gene.source = "PSEUDO"

	#
	# Now massage each transcript.
	#
	for i,j,hasCDS,hasRNA in transcripts:
	    # ----------------------
	    # If there is explicit transcript feature, convert it.
	    # Otherwise, create a placeholder (presumed transcript).
	    #
	    if hasRNA:
		# normal case. 
		t = self.currModel[i]
		#id = t.attributes['original'][iFID]
		id = t.attributes['original'][iTRANS]
	    else:
		# No transcript (gene segment from NCBI). Create a placeholder.
		t  = gff3.Feature(gene)
		id = "%s:presumedTranscript_%d" % \
		    (gene.attributes['ID'],self.getIDCounter("presumed"))
		# compute start/stop from subfeatures
		start = 900000000
		end = -1
		for u in xrange(i,j+1):
		    f = self.currModel[u]
		    start = min(start, f.start)
		    end   = max(end, f.end)
		t.start = start
		t.end = end
	    #
	    t.type = "mRNA"
	    t.attributes = {}
	    t.attributes['ID'] = id
	    t.attributes['Parent'] = gene.attributes['ID']
	    massaged.append(t)

	    #
	    exons = []
	    utrcds = []
	    iadjust = 0
	    if hasRNA:
		iadjust = 1
	    if hasCDS:
		# generate/add exons
		#
		curr_e = None
		for u in xrange(i+iadjust,j+1):
		    f = self.currModel[u]
		    utrcds.append(f)
		    if curr_e and curr_e.end+1 == f.start:
			# extend current exon
			curr_e.end = f.end
		    else:
			# add a new exon
			curr_e = gff3.Feature(self.currModel[u])
			curr_e.type = "exon"
			exons.append(curr_e)
	    else:
		# transcripts w/o any CDS are re-typed as "transcript",
		# all the UTRs become exons
		if t:
		    t.type = "transcript"
		for u in xrange(i+iadjust,j+1):
		    f = self.currModel[u]
		    #utrcds.append(f)
		    e = gff3.Feature(f)
		    e.type = "exon"
		    exons.append(e)

	    # For gene models on the - strand, the GFF3 standard is to
	    # list the subfeatures in reverse order.
	    #
	    if strand=='-':
		utrcds.reverse()
		exons.reverse()

	    #
	    # Refine the UTR types (to 5' or 3').
	    # Calculate CDS phase.
	    #
	    start_codon = None
	    stop_codon = None
	    first_cds = None
	    last_cds = None
	    utr = 'five_prime_UTR'
	    cdslen = 0
	    if hasCDS:
		# Only do this if there is a CDS. (Otherwise, cannot
		# define 5' vs 3'.)
		for f in utrcds:
		    if f.type == "UTR":
			f.type = utr
		    elif f.type=="CDS":
			utr = 'three_prime_UTR'
			last_cds = f
			if not first_cds:
			    first_cds = f
			    first_cds.phase = '0'
			phase = cdslen % 3
			if phase > 0:
			    phase = 3-phase
			f.phase = str(phase)
			cdslen += f.end - f.start + 1

	    #
	    # Generate start and stop codons.
	    #
	    if t and first_cds:
		if (first_cds.end - first_cds.start) >= 2:
		    start_codon = gff3.Feature(first_cds)
		    start_codon.type = "start_codon"
		    start_codon.phase = "."
		    if strand == "+":
			start_codon.end = start_codon.start + 2
		    else:
			start_codon.start = start_codon.end - 2
		    utrcds.insert(0,start_codon)
		if (last_cds.end - last_cds.start) >= 2:
		    stop_codon = gff3.Feature(last_cds)
		    stop_codon.type = "stop_codon"
		    stop_codon.phase = "."
		    if strand == "+":
			stop_codon.start = stop_codon.end - 2
		    else:
			stop_codon.end = stop_codon.start + 2
		    utrcds.append(stop_codon)

	    #
	    # Assemble all the pieces of this transcript in the proper order
	    #
	    keyfn = lambda x: (x.strand=="+" and 1 or -1)*x.start
	    exons.sort(  None, keyfn )
	    utrcds.sort( None, keyfn )
	    tmassaged = exons + utrcds
	    parent = t
	    if t is None:
		parent = gene
	    for f in tmassaged:
		forig = f.attributes['original']
		f.attributes = {}
		f.attributes['Parent'] = parent.attributes['ID']
	    #
	    # Add massaged transcript to model
	    #
	    massaged += tmassaged
	    

	# end of for each transcript loop
	self.currModel = massaged


    #-----------------------------------------------------
    def startNewGene(self):
        if self.currModel:
	    self.finishGene()
	self.currModelID = self.currTokens[iFID][7:] # strip off "GeneID:"
	self.resetIDCounters()
	self.continueGene()

    #-----------------------------------------------------
    def continueGene(self):
	g = self.convertLine(self.currTokens)
	self.currModel.append(g)

    #-----------------------------------------------------
    def finishGene(self):
	if len(self.currModel) > 0:
	    self.massageModel()
	    for g in self.currModel:
		print gff3.format(g),
        self.currModel = []
	self.currModelID = None

    #-----------------------------------------------------
    def main(self, argv):
	self.parseArgs(argv)
        self.openFiles()
	currID = None
	for self.currLine in self.iFile:
	    self.currLineNum += 1
	    if self.currLine[0:1] == HASH:
	        continue
	    self.currTokens = self.currLine.split(TAB)
	    self.currTokens[-1] = self.currTokens[-1].strip()

	    # only interested in Black 6 at this time
	    if self.currTokens[iGROUP] != self.options.assembly:
	        continue

	    # new gene?
	    if self.currTokens[iFID] != currID:
		self.startNewGene()
		currID = self.currTokens[iFID]
	    else:
		self.continueGene()
	#
	self.finishGene()

#-----------------------------------------------------
#-----------------------------------------------------
if __name__=="__main__":
    Entrez2GFF3().main(sys.argv)
