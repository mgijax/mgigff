"""
#
# mgigff.py
#
# Top level script for the unified MGI GFF File pipeline.
#
# To run this script, edit the default config file (config.cfg
# in the parent directory), then:
#	% python mgigff.py
#
# You can also create/read additional config file(s):
#	% python mgigff.py file1.cfg file2.cfg ...
# Later config files can override settings from earlier config files.
#
#
# HISTORY
#   jer - April 2011
#   Created. Phase 1 (genes with models) only.
#
#   jer - August 2011
#   Added phase 2 (genes without models, but with sequences).
#
"""

#----------------------------------------------------------------------------
import os
import types
import time
import sys
import ConfigParser
import re
import mgiadhoc as db

import gff3
import psl
import logging
import getopt

#----------------------------------------------------------------------------
COMMA	= ','
TAB	= '\t'
NL	= '\n'
SP	= ' '
PIPE	= '|'
DOT	= '.'

#----------------------------------------------------------------------------
class MGIGFFMaker(object):
    def __init__(self,args):
	self.myfile = None	# this script file
	self.home = None	# where I'm installed
	#
	self.defaultConfigFile = None # config file if none specified
	self.configParser = None# config file parser
	self.config = None	# config data extracted into a dict
	#
	self.mgiDate = "?"
	self.providers = []	# list of model providers
	self.k2provider = {}	# MGI logicalDb key-to-provider index
	#
	self.mid2models = {}  # { mgiid -> { ldbname-> [ id1, id2, ... ] } }
	self.mid2marker = {}  # { mgiid -> marker result record }
	self.model2mids = {}  # { ldbname->{ id  -> [ mgiid1, mgiid2, ... ] } }
	self.allModels = set()
	#
	self.count = 0	      # output line counter
	#
	#
	self.doPhase1 = True
	self.doPhase2 = True
	#
	self.init(args)

    #---------------------------------------------------------------
    #---------------------------------------------------------------

    def doCmd(self, cmd):
	"""
	Executes a shell command.
	"""
	logging.info(cmd)
	os.system(cmd)

    def showConfig(self, level=logging.INFO):
	"""
	Logs the current configuration.
	"""
	import pprint
	s = pprint.pformat(self.config)
	logging.log(level, s)
	sys.stdout.write(s)
	sys.stdout.write('\n')

    def go(self):
	"""
	Main.
	"""
	if self.doPhase1:
	    self.phase1()
	if self.doPhase2:
	    self.phase2()

	self.processSorted()
	self.generateExomeFile()
	self.moveToArchive()
	logging.info("mgigff finished.\n\n")

    #---------------------------------------------------------------
    #---------------------------------------------------------------

    def init(self, args):
	self.initContext()
	self.initConfig(args)
	self.initLogging()
	self.initSortPipe()
	self.initProviders()

    def initContext(self):
	"""
	Find where I am installed, find my config file, etc...
	"""
	self.myfile = os.path.abspath(__file__)
	(dir,modfile) = os.path.split(self.myfile)
	(self.home,bin) = os.path.split(dir)
	self.defaultConfigFile = os.path.join(self.home, "config.cfg")

    def initConfig(self, args):
	# Config file(s) may be specified on command line.
	opts, posargs = getopt.getopt( args, "v12d:c:s:", ["version", "define=", "config=", "section="] )
	posargs = [self.defaultConfigFile] + posargs
	# Read all config files 
	self.cp = ConfigParser.ConfigParser()
	for cfname in posargs:
	    if not os.path.exists(cfname):
		raise IOError("Config file not found: " + cfname)
	    self.cp.read(cfname)
	self.versionOnly = False
	#
	# Now look through options specified.
	dfltsection = "MGIGFF"
	for n,v in opts:
	    if n in ["-d","--define"]:
		# cmd line definitions
		# v has the form "name=value"
		m = re.match( r'\s*([\w.]+:)?([\w.]+)\s*=(.*)', v)
		if not m:
		    raise RuntimeError("Bad format in definition: %s" % v)
		section = m.group(1) and m.group(1)[:-1] or dfltsection
		if not self.cp.has_section(section):
		    self.cp.add_section(section)
		self.cp.set(section, m.group(2), m.group(3))
	    elif n in ["-s","--section"]:
		dfltsection = v
	    elif n in ["-c","--config"]:
		if not os.path.exists(v):
		    raise IOError("Config file not found: " + v)
		self.cp.read(v)
	    elif n == "-v":
		# print version/config info and exit
		self.versionOnly = True
	    elif n == "-1":
		self.doPhase2 = False
	    elif n == "-2":
		self.doPhase1 = False

	if not (self.doPhase1 or self.doPhase2):
	    self.doPhase1 = True
	    self.doPhase2 = True

	#
	# Load all parms into a plain dict. Parse dotted
	# section names and interpret/translate to dict hierarchy.
	# E.g., section [Foo.Bar] ===>  { "Foo" : { "Bar" : {} } }
	tld = dict()
	ks = os.environ.keys()
	ks.sort()
	for k in ks:
	    #print '[', k, '=', os.environ[k]
	    pass
	for s in self.cp.sections():
	    cd = tld
	    for sn in s.split("."):
		cd = cd.setdefault(sn,{})
	    for o in self.cp.options(s):
		v = self.cp.get(s,o,False,os.environ)
		cd[o] = v
	self.config = tld
	#
	cfg = self.config['MGIGFF']
	#
	self.dataDirectory = cfg['datadirectory']
	if not (os.path.exists(self.dataDirectory) and os.path.isdir(self.dataDirectory)):
	    raise RuntimeError("Data directory not found or is not a directory: %s"%self.dataDirectory)
	self.workingDirectory = cfg['workingdirectory']
	if not os.path.exists(self.workingDirectory):
	    os.makedirs(self.workingDirectory)
	elif not os.path.isdir(self.workingDirectory):
	    raise RuntimeError("Working directory is not a directory: %s"%self.workingDirectory)
	self.distributionDirectory = cfg.get('distributiondirectory','')
	if self.distributionDirectory and not (os.path.exists(self.distributionDirectory) and os.path.isdir(self.distributionDirectory)):
	    raise RuntimeError("Distribution directory not found or is not a directory: %s"%self.distributionDirectory)

	tstamp = time.strftime("%Y%m%d") # e.g. "20110503" for May 3, 2011
	self.mgigffFileC = cfg['mgiprefix'] + ".gff3"			# MGI.gff3
	self.mgigffFileV = cfg['mgiprefix'] + ".%s.gff3" % tstamp	# MGI.20110503.gff3
	self.exomeFileC  = cfg['exomeprefix'] + ".gff3"			# MGI.exome.gff3
	self.exomeFileV  = cfg['exomeprefix'] + ".%s.gff3" % tstamp	# MGI.exome.20110503.gff3
	#
	self.readmeFileC = cfg['mgireadme']
	self.readmeFileV = cfg['mgireadme']

	self.furlC = os.path.join(cfg['baseurl'], self.mgigffFileC)	# http://ftp.../pub/mgigff/MGI.gff3
	self.furlV = os.path.join(cfg['baseurl'], self.mgigffFileV)	# http://ftp.../pub/mgigff/MGI.20110503.gff3
	self.rurlC = os.path.join(cfg['baseurl'], self.readmeFileC)
	self.rurlV = os.path.join(cfg['baseurl'], self.readmeFileV)

	# thresholds
	self.phase2minidentity = int(cfg.get('phase2minidentity',0))
	self.phase2minlength   = int(cfg.get('phase2minlength',0))
	self.phase2maxdistance = int(cfg.get('phase2maxdistance',0))

	self.logFile = cfg['logfile']
	self.orphanMgiFile = os.path.join(self.workingDirectory, "LOG_orphan_MGI.txt")
	self.orphanModelFile = os.path.join(self.workingDirectory, "LOG_orphan_Model.txt")
	self.notFoundFile = os.path.join(self.workingDirectory, "LOG_not_found.txt")

	# file of seqids (+assoc. gene info) for "orphan" mgi genes.
	self.phase2seqids = os.path.join(self.workingDirectory, "phase2seqids.txt" )

	# the sequences (in FASTA format) for the seqids in self.phase2seqids
	self.phase2sequences = os.path.join(self.workingDirectory, "phase2sequences.fa")

	# the Blat results (in psl format) for the sequences in self.phase2sequences
	self.phase2blat = os.path.join(self.workingDirectory, "phase2blat.psl")

	# the pslReps results (in psl format) for the alignments in self.phase2blat
	self.phase2psl = os.path.join(self.workingDirectory, "phase2pslreps.psl")
	# the other file that pslReps outputs - not used so far...
	self.phase2pslr = os.path.join(self.workingDirectory, "phase2pslreps.psr")

	# the alignments from self.phase2psl, decorated with associated gene info (from self.phase2seqids)
	self.phase2merged = os.path.join(self.workingDirectory, "phase2merged.psl")

	# log file of alignments/genes rejected from the result. Possible reasons:
	#	1. multiple alignments for one sequence
	#	2. insufficient quality (%identity<threshold)
	#	3. chromosome disagreements
	self.phase2rejects = os.path.join(self.workingDirectory, "phase2rejects.txt")
	self.phase2rejectsgff = os.path.join(self.workingDirectory, "phase2rejects.gff3")

	# log file of genes added by phase 2.
	self.phase2added = os.path.join(self.workingDirectory, "phase2added.txt")

    def initLogging(self):
	"""
	Initialize logging.
	"""
	#
	ll = self.config['MGIGFF'].get('loglevel','DEBUG')
	try:
	    self.logLevel = int(ll)
	except ValueError:
	    try: 
		self.logLevel = vars(logging)[ll]
	    except KeyError:
		self.logLevel = logging.DEBUG
	#
	self.logFile = os.path.join(self.workingDirectory,self.logFile)
	logging.basicConfig(
	    level = self.logLevel,
	    format='%(asctime)s %(levelname)s %(message)s',
	    filename=self.logFile,
	    filemode='a')
	logging.info("\n\nMGIGFF started...")

	self.showConfig()
	if self.versionOnly:
	    sys.exit(0)

    def initSortPipe(self):
	"""
	Opens a R/W pipe to a sort command. We write modified GFF lines
	to the sort's input (we add certain tags), we sort on the desired
	fields, then we read the sort's output and strip off the added tags.

	To achieve the desired sort order, each GFF3 line is tagged
	with the associated gene's start position and a global, incremented 
	counter. We sort by (1) chromosome (i.e. col 1 of the GFF), (2) gene
	start position, and (3) counter.
	"""
	sortcmd = '%s -k3,3 -k1,1n -k2,2n' % self.config['MGIGFF']['gnusort']
	self.sortpipein,self.sortpipeout = os.popen2(sortcmd)

    def initProviders(self):
	"""
	Creates a Provider object for each configured provider. 
	(Each provider has a section in the config file named provider.XXX,
	where XXX is the provider's name. See config file.)
	"""
        for (pname,provider) in self.config.get("provider", {}).iteritems():
	    p = Provider(pname,provider,self)
	    self.providers.append(p)
	    self.k2provider[p.key] = p

    def writeToSortPipe(self, feature, currmgistart):
	# Write MGI feature to the sort's input
	self.count += 1
	self.sortpipein.write("%s\t%d\t%s" % \
	    (currmgistart,self.count,gff3.format(feature)) )

    #---------------------------------------------------------------
    #---------------------------------------------------------------

    def phase1(self):
	logging.info("Starting phase 1...")
	self.initPhase1Logs()
	self.getMGIGeneInfo()
	self.convertProviderFiles()
	self.mergeProviderFiles()
	self.closePhase1Logs()
	logging.info("End of phase 1.")

    def initPhase1Logs(self):
	self.orphanMgiFd = open(self.orphanMgiFile, 'w')
	self.orphanModelFd = open(self.orphanModelFile, 'w')
	self.notFoundFd = open(self.notFoundFile, 'w')
	
    def closePhase1Logs(self):
	self.orphanMgiFd.close()
	self.orphanModelFd.close()
	self.notFoundFd.close()
	
    def convertProviderFiles(self):
	"""
	Converts each provider's raw input file to standard GFF3 format.
	"""
	logging.info("Converting provider files to GFF3.")
        for p in self.providers:
	    p.convert()

    def mergeProviderFiles(self):
	"""
	Reads each provider's GFF3 (converted) file, and writes it to 
	the sort pipe's input. Along the way it:
	    - attaches any associated MGI IDs
	    - prepends two columns at the beginning of each line which are used
	      to drive the sorting. The columns are: (1) the startposition of the 
	      associated mgi gene and (2) a global counter incremented for each
	      line. Sorting on col 1 then col 2 has the effect of grouping all
	      lines associated with a gene together, ordering those groups by
	      startposition of the associated gene, preserving the input order
	      within each group.
	    - moves original type into col 9 (for genes/models), and puts
	      the provider's name into col 2.
	    - keeps running tally of min/max coordinates for features of a gene
	    - writes to the sort pipe's input
	    - detects/reports models that are referenced in MGI but not found
	      in provider file
	"""
	logging.info("Merging GFF3 provider files with MGI data.")
	mids = None
	currmgistart = 0
	# foreach provider
	for provider in self.providers:
	    # get the dict mapping model ids to MGI ids
	    id2mgi = self.model2mids[provider.name]
	    # open the GFF3 file for this provider
	    file = open(provider.outfile, 'r')
	    currModel = None
	    isPseudo = False
	    for line in file:
		if line.startswith(gff3.COMMENT_CHAR):
		    continue
		feature = gff3.Feature(line)
	        if feature is None:
		    continue
		# Remove Name attribute from all provider features
		feature.attributes.pop('Name',None)
		#
		fid=feature.attributes.get('ID',None)
		if fid:
		    feature.attributes['ID'] = "%s:%s" % (provider.name,fid)
		pid=feature.attributes.get('Parent',None)
		if pid:
		    feature.attributes['Parent'] = "%s:%s"%(provider.name,pid)
		#
		if feature.type in ['gene','pseudogene']:
		    #
		    currModel = feature
		    # FIXME. The following weirdness compensates for gtf2gff3's
		    # usage of cols 2 and 3. (This odd usage has been mimicked in
		    # entrez2gff3.py so this code can treat them uniformly.)
		    if feature.type == 'gene':
			s = feature.source
			feature.attributes['bioType'] = s
			isPseudo = ("pseudo" in s.lower())
			if isPseudo:
			    feature.type = 'pseudogene'
		    # For each top-level feature (ie. for each model)
		    # find its corresponding MGI gene(s). If none found,
		    # write entry in orphan model file, but still output
		    # the model to the GFF file.
		    self.allModels.add( feature.attributes['ID'] )
		    #
		    mids = id2mgi.pop(fid,None)
		    if mids:
			feature.attributes['Dbxref'] = map(lambda x:'MGI:%s'%x,mids)
			try:
			    mkr = self.mid2marker[mids[0]]
			except KeyError:
			    logging.error("Internal error: marker lookup failed. " + \
				"key=%s" % mids[0])
			currmgistart = mkr.attributes['sortKey'] # for sorting
		    else:
			self.orphanModelFd.write(fid+NL)
			currmgistart = feature.start # for sorting
		else:
		    # For bottom-level features, remove their (unneeded) ID attribute
		    if feature.type not in ["transcript","mRNA"]:
			feature.attributes.pop('ID',None)
		    # Propagate MGI ids to all subfeatures
		    gid = currModel.attributes['ID']
		    if mids:
			feature.attributes['Dbxref'] = \
				map(lambda x:'MGI:%s'%x,mids) + [gid]
		    else:
			feature.attributes['Dbxref'] = [gid]

		    #
		    if isPseudo:
			feature.type = "pseudogenic_"+feature.type

		# Accumulate min/max coords for the genes
		if mids:
		    for mid in mids:
			mkr = self.mid2marker.get(mid,None)
			mkr.start = min(mkr.start,feature.start)
			mkr.end   = max(mkr.end,  feature.end)

		#
		feature.source = provider.name

		# Write feature to the sort's input
		self.writeToSortPipe(feature, currmgistart)

	    # end for line in provider file

	    # anything left in id2mgi here represents model IDs associated with
	    # genes in MGI, but not found in the provider file.
	    ks = id2mgi.keys()
	    ks.sort()
	    for k in ks:
		self.notFoundFd.write("%s (%s)\n" % (k, COMMA.join(id2mgi[k])))

	# end for each provider

    #---------------------------------------------------------------
    #---------------------------------------------------------------

    def phase2(self):
	logging.info("Starting phase 2...")

	self.initPhase2Logs()

	# 1. Query MGI for seqids of genes w/out gene models
	logging.info("Querying MGI for seqids for orphan genes...")
	c = Converter('phase2query', self.config['converter']['PHASE2QUERY'], self)
	c.convert( self.orphanMgiFile, self.phase2seqids )

	# 2. Fetch the sequences from Entrez into a fasta file.
	logging.info("Retrieving sequences from Entrez...")
	c = Converter('seq2fa', self.config['converter']['SEQID2FA'], self)
	c.convert( self.phase2seqids, self.phase2sequences )

	# 3. Blat the sequences against the assembly.
	logging.info("Blat'ing sequences against genome...")
	c = Converter( 'blat', self.config['converter']['GFCLIENT'] , self)
	c.convert( self.phase2sequences, self.phase2blat )

	# 4. Run pslReps to filter the alignments to single best hits.
	logging.info("Filtering alignments...")
	c = Converter( 'pslReps', self.config['converter']['PSLREPS'] , self)
	c.convert( self.phase2blat, self.phase2psl+" "+self.phase2pslr )

	# 5. Attach MGI gene info to each alignment in psl file
	logging.info("Merging MGI gene info into .psl file...")
	c = Converter( 'phase2merge', self.config['converter']['PHASE2MERGE'] , self)
	c.convert( self.phase2seqids + ' ' + self.phase2psl, self.phase2merged )

	# 6. Read the psl file, do filtering/sanity checking, and convert to GFF features.
	self.phase2convert()

	self.closePhase2Logs()
	logging.info("End of phase 2.")

    def initPhase2Logs(self):
	self.phase2rejectsFd = open(self.phase2rejects, 'w')
	self.phase2rejectsgffFd = open(self.phase2rejectsgff, 'w')
	self.phase2addedFd = open(self.phase2added, 'w')

    def closePhase2Logs(self):
	self.phase2rejectsFd.close()
	self.phase2rejectsgffFd.close()
	self.phase2addedFd.close()

    def phase2convert(self):
	logging.info("Converting alignments to features...")
	currMgiId = None
	currAlignments = {}
	for a in psl.iterate(self.phase2merged):
	    # The alignments in this file have been sorted by gene.
	    mgiid,seqid,xxx = a.qName.split(PIPE,2)
	    if mgiid != currMgiId:
		self.processAlignments(currMgiId,currAlignments)
		currMgiId = mgiid
		currAlignments = {}
	    currAlignments.setdefault(seqid,[]).append(a)
	self.processAlignments(currMgiId,currAlignments)

    def prependDatabase(self, id):
	if id.startswith("MGI:"):
	    return "MGI:"+id
	elif id[2:3] == "_":
	    return "RefSeq:"+id
	else:
	    return "GenBank:"+id

    def processAlignments(self, mgiid, aligns):
	if mgiid is None or len(aligns) == 0:
	    return

	qName = None
	symbol = None
	name = None
	chr = None
	bioType = None

	minCoord = 999999999
	maxCoord = 0

	badChr = False

	dbxrefs = []

	def makeGffFeature(a, type, score="."):
	    return gff3.Feature([
		    a.tName[3:],
		    "phase2",
		    type,
		    a.tStart,
		    a.tEnd,
		    score,
		    a.strand,
		    '.',
		    {
			'seqid' : seqid,
			'qname' : a.qName
		    } ])

	for seqid,saligns in aligns.items():
	    if qName is None:
		qName = saligns[0].qName
		xxx,xxx,bioType,symbol,name,chr,xxx = qName.split(PIPE,6)
	    if len(saligns) > 1:
		# Report and remove seqids with multiple alignments
		self.phase2rejectsFd.write(
		  "IGNORING SEQUENCE (%s) FOR GENE (%s) - MULTIPLE ALIGNMENTS\n"%(seqid,mgiid))
		for a in saligns:
		    self.phase2rejectsFd.write(str(a))
		    f = makeGffFeature(a, "multipleAlignments")
		    self.phase2rejectsgffFd.write(str(f))
		del aligns[seqid]
		continue

	    # replace list with single item
	    a = saligns[0]
	    aligns[seqid] = a
	    #
	    hspLen = a.qEnd-a.qStart
	    pident = 100 * float(a.matches)/hspLen
	    plen = 100 * float(hspLen)/a.qSize
	    if (self.phase2minidentity>0 and pident < self.phase2minidentity) \
	    or (self.phase2minlength>0 and plen < self.phase2minlength):
		self.phase2rejectsFd.write(
		    "\nREJECTING ALIGNMENT - SEQUENCE (%s) GENE (%s) PCT identity/length(%1.2f/%1.2f)\n" \
		    % (seqid,mgiid,pident,plen))
		self.phase2rejectsFd.write(str(a))
		f = makeGffFeature(a, "pctIdentity", "%1.2f/%1.2f"%(pident,plen))
		self.phase2rejectsgffFd.write(str(f))
		del aligns[seqid]
		continue
	    #
	    badChr = badChr or (chr != a.tName)
	    #
	    minCoord = min(minCoord, a.tStart)
	    maxCoord = max(maxCoord, a.tEnd-1)
	    #
	    dbxrefs.append(self.prependDatabase(seqid))

	# if chromosomes disagree, reject the lot
	if badChr:
	    self.phase2rejectsFd.write(
	        "\nREJECTING GENE - CHROMOSOMES DISAGREE (%s)\n"%mgiid)
	    for s,a in aligns.iteritems():
		self.phase2rejectsFd.write(str(a))
		f = makeGffFeature(a, "chrMismatch")
		self.phase2rejectsgffFd.write(str(f))
	    return

	mid = self.prependDatabase(mgiid)

	# create feature for MGI gene
	gf = gff3.Feature(
	    [ chr[3:], 'MGI', self.convertMGIType(bioType), minCoord, maxCoord, '.', '.', '.', {} ])
	gf.ID = mgiid
	gf.Name = symbol
	gf.mgiName = name
	gf.bioType = bioType
	gf.Dbxref = dbxrefs
	# register it
	self.mid2marker[mgiid] = gf

	if len(aligns) == 0:
	    self.phase2rejectsFd.write(
	        "\nNO MODELS REMAIN AFTER FILTERING - GENE (%s)\n"%mgiid)
	    self.phase2rejectsFd.write(str(gf))
	    return

	# write the feature representing the MGI gene
	self.writeToSortPipe( gf, gf.start )
	sids = []
	for s,a in aligns.iteritems():
	    ss = self.prependDatabase(s)
	    score = a.matches/float(a.qSize)
	    # write the "match" feature (this alignment)
	    mg = gff3.Feature(
	        [a.tName[3:], 'Blat', 'match', a.tStart, a.tEnd-1, "%1.2f"%score, a.strand, '.', {}])
	    mg.ID = ss
	    mg.Dbxref = mid
	    self.writeToSortPipe( mg, gf.start )
	    self.allModels.add(ss)
	    sids.append(ss)
	    for i in range(a.blockCount):
		sz = a.blockSizes[i]
		st = a.tStarts[i]+1
		se = st + sz - 1
		# write "match-part" feature (alignment fragment)
		mpf = gff3.Feature([a.tName[3:], 'Blat', 'match-part', st, se, '.', a.strand, '.', {}])
		mpf.Parent = ss
		mpf.Dbxref = mid
		self.writeToSortPipe( mpf, gf.start )
	self.phase2addedFd.write( mgiid + TAB + COMMA.join(sids) + NL )

    #---------------------------------------------------------------
    #---------------------------------------------------------------

    def generateHeader(self):

	MGIGFF_HEADER = '''#
	# File: %s.gz
	# Description: MGI Unified Mouse Gene Catalog. See: %s
	# URL: %s.gz
	# Latest-version: %s.gz
	# Generated: %s
	# Genome build: %s
	# Contact: %s
	#
	# MGI database:
	#   Last updated: %s
	''' % (
	    self.mgigffFileV, 
	    self.rurlV,
	    self.furlV,
	    self.furlC,
	    time.asctime(time.localtime(time.time())), 
	    self.config['MGIGFF']['genomebuild'],
	    self.config['MGIGFF']['contactemail'],
	    self.mgiDate)
	
	for p in self.providers:
	    MGIGFF_HEADER += '''#
	    # Provider: %s
	    #   File: %s
	    #   Modified: %s
	    ''' % (p.name, os.path.basename(p.file), 
		   time.asctime(time.localtime(p.fstat.st_mtime)))
	
	p = re.compile(r'^\s+', re.M)
	MGIGFF_HEADER = p.sub('',MGIGFF_HEADER)

	return gff3.HEADER + MGIGFF_HEADER

    def processSorted(self):
	"""
	Final pass over the sorted features. Strips off the tags
	added for sort step, and replaces MGI marker start/end
	with computed values. Decides when/where to output "###" dividers.
	"""
	logging.info("Merging/sorting phase 1 and 2 results.")
	gene_re= re.compile(r'^[^\t]*\t([^\t]*)\t(gene|pseudogene|sequence_feature)')
	self.sortpipein.close()
	fd = open(os.path.join(self.workingDirectory,self.mgigffFileV), 'w')
	fd.write(self.generateHeader())
	needDivider = False
	lastTop = None
	for line in self.sortpipeout:
	    #
	    # Split off and discard the sort columns, then
	    # instantiate feature from remainder.
	    #
	    tokens = line.split(TAB, 2)
	    line = tokens[2]
	    f = gff3.Feature(line)
	    #
	    # Is it a "top level" feature? Here we define top level as
	    # 1. an MGI gene with at least one model, or
	    # 2. a gene/pseudogene feature from a provider, with no MGI id
	    #
	    isTop = f.type in ("gene","pseudogene","sequence_feature") and \
		    (f.source == "MGI" or "MGI" not in f.attributes.get("Dbxref",""))
	    
	    writeIt = True
	    if isTop and f.source == "MGI":
		# calculate top level start/stop positions.
		mkr = f
		mkr.attributes.pop('sortKey', None)
		mgiid = mkr.ID
		mkr.ID = "%s:%s"%(mkr.source,mkr.ID)

		# set start/end coords to min/max values
		minmax = self.mid2marker[mgiid]
		mkr.start = minmax.start
		mkr.end   = minmax.end
		# before writing out this gene, make sure at least
		# one of its named models was actually seen
		writeIt = False
		xrs = mkr.attributes.get('Dbxref',[])
		if type(xrs) is types.StringType:
		    xrs = [xrs]
		for id in xrs:
		    if id in self.allModels:
			writeIt = True
			break
		if not writeIt:
		    logging.warn("Models not found: %s\t%s\t%s\t%s" % \
			(mgiid, mkr.attributes['Name'], mkr.attributes['bioType'], 
			 str(mkr.attributes.get('Dbxref',[]))))
		    continue
	    #
	    # If f is a top-level feature (TLF), and the previous TLF doesn't overlap this one,
	    # write a separator before f.
	    # (Big "Hmmmmm...." here. Are these the right criteria?)
	    if isTop:
		if lastTop and (lastTop.seqid != f.seqid or lastTop.end < f.start):
		    fd.write('###\n')
		lastTop = f
	    fd.write(gff3.format(f))
	#
	# end of loop
	#
	fd.close()

    def generateExomeFile(self):
	import mgiexome
	mgiexome.main(os.path.join(self.workingDirectory,self.mgigffFileV), 
		      os.path.join(self.workingDirectory,self.exomeFileV))

    def moveToArchive(self):
	if not self.distributionDirectory:
	    return
	logging.info("Copying files to distribution directory.")
	sf = os.path.join(self.workingDirectory,self.mgigffFileV)
	df = os.path.join(self.distributionDirectory,self.mgigffFileV)
	dlf= os.path.join(self.distributionDirectory,self.mgigffFileC)
	#
	# remove destination files from distrib directory, if they exist
	self.doCmd( "rm -f %s.gz %s.gz" % (df,dlf) )
	#
	# compress the GFF file and write to distrib directory
	self.doCmd("gzip -c %s > %s.gz" % (sf,df))
	#
	# Create symbolic link
	self.doCmd( "cd %s; rm -f %s.gz; ln -s %s.gz %s.gz" % \
	   (self.distributionDirectory,self.mgigffFileC,self.mgigffFileV,self.mgigffFileC))

    #-------------------------------------------------
    # FIXME: The remainder is all for getting gene data from MGI.
    # Should be split out...
    #-------------------------------------------------

    def parseMGIdate(self, r):
	self.mgiDate = r['lastdump_date'].strftime('%b %d, %Y')

    def parseIDs(self, r):
	"""
	Processes each MGI/provider ID pair. Creates bi-directional, 
	multi-valued mappings.
	"""
	ldbname = self.k2provider[r['_logicaldb_key']].name
	mid = r['mgiid']
	id = r['accid']
	self.mid2models.setdefault(mid, {}).setdefault(ldbname,[]).append(id)
	self.model2mids.setdefault(ldbname, {}).setdefault(id, []).append(mid)

    def convertMGIType(self, t):
	if "segment" in t:
	    return "sequence_feature"
	elif "pseudo" in t:
	    return "pseudogene"
	elif "gene" in t:
	    return "gene"
	else:
	    return "sequence_feature"

    def parseMarker(self, m):
	"""
	Processes each mouse marker. Turns it into a Feature, registers it,
	and writes it to the sort pipe's input.
	"""
	mid = m['mgiid']
	ftype = m['featuretype']
	ftype2 = self.convertMGIType(ftype)
	sortKey = str(int(m['startcoordinate']))
	attrs = {
	    "ID" : mid,
	    "Name" :  m['symbol'],
	    "bioType" : ftype,
	    "mgiName" : m['name'],
	    "sortKey" : sortKey, # temporary, will be removed later
	    }

	noModels = True
	for n,ids in self.mid2models.get(mid,{}).iteritems():
	    noModels = False
	    #attrs["%s_ID" % n.lower()] = ids
	    attrs.setdefault('Dbxref',[]).extend( map(lambda x:'%s:%s'%(n,x), ids) )
	
	if noModels:
	    logging.error("No models found?? %s"%mid)

	mgiMarker = gff3.Feature([
	    m['chromosome'],
	    "MGI",
	    ftype2,
	    int(m['startcoordinate']),
	    int(m['endcoordinate']),
	    ".",
	    m['strand'],
	    ".",
	    attrs
	    ])
	self.mid2marker[mid] = mgiMarker
	self.writeToSortPipe(mgiMarker, int(m['startcoordinate']))

    def parseOrphanMarker(self, m):
	self.orphanMgiFd.write( TAB.join( [
	    m['mgiid'],
	    m['symbol'],
	    m['featuretype'],
	    ]) + NL)

    def getMGIGeneInfo(self):
	"""
	Queries MGI for MGI genes and gene model associations.
	"""
	logging.info("Getting data from MGI...")
	cfg = self.config["MGIGFF"]
	#
	db.HOST     = cfg.get("dbserver", db.HOST)
	db.DATABASE = cfg.get("dbdatabase", db.DATABASE)
	db.USER     = cfg.get("dbuser", db.USER)
	db.PASSWORD = cfg.get("dbpassword", db.PASSWORD)

	logging.info("Database connection: Host(%s) Database(%s) User(%s) " % \
	   (db.HOST,db.DATABASE,db.USER))

	conn = db.connect()

	#
	# Retrieve MGI update date
	#
	db.sql('''
	    select lastdump_date
	    from MGI_dbinfo
	    ''', self.parseMGIdate, conn)

	#
	# Generate temp table of  MGI_id/model_id pairs
	#
	db.sql('''
	    create temporary table _ids (
		accid 		varchar(30),
		_LogicalDB_key	int,
		mgiID		varchar(30),
		_Marker_key	int
	    )
	    ''', None, conn)

	db.sql('''
	    insert into _ids
	    select distinct 
	       a.accid, a._LogicalDB_key, a2.accID, a._Object_key
	    from ACC_Accession a, ACC_Accession a2, MRK_Marker mm
	    where a._LogicalDB_key in (%s)
	    and a._Object_key = a2._Object_key
	    and a._MGIType_key = 2
	    and a._MGIType_key = a2._MGIType_key
	    and a2._LogicalDB_key = 1
	    and a2._MGIType_key = 2
	    and a2.preferred = 1
	    and a2._Object_key = mm._Marker_key
	    and mm._marker_status_key in (1,3)  /* include interim symbols */
	    ''' % COMMA.join(map(str,self.k2provider.keys())), None, conn)

	#
	db.sql('''
	    create index ix0 on _ids(_Marker_key)
	    ''',None, conn)

	# retrieve the MGI/provider id pairs (to build mappings)
	db.sql('''
	    select * from _ids
	    ''', self.parseIDs, conn)

	#
	# Retrieve additional information about the mouse genes. From
	# the first query. (I.e., all returned genes have at least one
	# model association.)
	#
	db.sql('''
	    select 
		a.accID as "mgiid", 
		v._Marker_key, 
		v.symbol, 
		v.name, 
		c.genomicchromosome as chromosome, 
		c.startCoordinate, 
		c.endCoordinate, 
		c.strand,
		mcv.term as featuretype
	    from 
		ACC_Accession a,
		MRK_Marker v, 
		MRK_Location_Cache c, 
		MRK_MCV_Cache mcv
	    where v._Marker_key in (
		select _Marker_key
		from _ids)
	    and v._Marker_key = a._Object_key
	    and a._LogicalDB_key = 1
	    and a._MGIType_key = 2
	    and a.preferred = 1
	    and a.private = 0
	    and v._Marker_key = c._Marker_key
	    and c.startCoordinate is not null
	    and v._Marker_key = mcv._Marker_key
	    and mcv.qualifier='D'
	    ''', self.parseMarker, conn)
	
	#
	# Retrieve MGI genes and pseudogenes that do not have models.
	#
	db.sql('''
	    select 
		a.accID as "mgiid", 
		v.symbol, 
		c.genomicchromosome as chromosome, 
		c.startCoordinate, 
		c.endCoordinate, 
		c.strand,
		mcv.term as "featuretype"
	    from 
		ACC_Accession a,
		MRK_Marker v, 
		MRK_Location_Cache c, 
		MRK_MCV_Cache mcv
	    where v._Marker_key in (
		select distinct _Marker_key
		from MRK_MCV_Cache
		where term in ('gene','pseudogene'))
	    and v._Marker_key not in (
		select _Marker_key
		from _ids)
	    and v._Marker_key = a._Object_key
	    and a._LogicalDB_key = 1
	    and a._MGIType_key = 2
	    and a.preferred = 1
	    and a.private = 0
	    and v._Marker_key = c._Marker_key
	    and v._Marker_key = mcv._Marker_key
	    and mcv.qualifier='D'
	    ''', self.parseOrphanMarker, conn)
	
	#
	conn.close()

#----------------------------------------------------------------------------
class Provider(object):
    """
    A Provider encapsulates a source of gene models; it must have
    a corresponding logical database in MGI. In addition to its
    logical db key, a provider has a name, a directory where its
    raw input file is found, and a regular expression for finding 
    that file within that directory. It also has a Converter, that
    will convert its raw data file into a standard GFF3 file.
    """
    def __init__(self, name, parms, gffMaker):
        self.name = name
	self.key = int(parms['mgildbkey'])
	self.dir = parms.get('datadirectory', gffMaker.dataDirectory)
	self.file_re = re.compile(parms['fileregex'])
	self.file = None
	self.outfile = os.path.join(gffMaker.workingDirectory, name+".gff3")
	for fn in os.listdir(self.dir):
	    if self.file_re.match(fn):
	        self.file = os.path.join(self.dir,fn)
		break
	if self.file is None:
	    raise IOError("No input file for provider " + self.name)
	self.fstat = os.stat(self.file)
	n = parms.get('converter')
	self.converter = Converter(n, gffMaker.config['converter'][n], gffMaker)

    def convert(self):
	self.converter.convert(self.file, self.outfile)


#----------------------------------------------------------------------------
class Converter(object):
    """
    A converter is something you can run to translate a gene model file
    in some input format into a standardized GFF3 output.
    """
    def __init__(self, name, parms, gffMaker):
        self.name = name
	self.executable = parms['executable']
	self.cmdTmplt = parms["command"]
	self.logFile = gffMaker.logFile

    def convert(self, infile, outfile):
	cmd = self.cmdTmplt.replace('INPUT',infile).replace('OUTPUT', outfile)
	logging.info("%s: Command=%s" % (self.name,cmd))
	if newer(infile.split(), outfile.split()):
	    status=os.system(cmd + ' 2>> ' + self.logFile)
	    if status != 0:
		raise RuntimeError("Converter failed. status=%s cmd=%s"%(status,cmd))
	else:
	    logging.info("Command was short circuited.")


#----------------------------------------------------------------------------
def newer( sourceFile, targetFile ):
    """
    # Utility function that returns True iff any sourceFile has been modified more
    # recently than some targetFile (or if some targetFile does not exist). If any 
    # sourceFile does not exist, throws an IOException.
    # Args:
    #    sourceFile (string or list) Path to source file, or list of paths.
    #    targetFile (string or list) Path to target file, or list of paths.
    # Returns:
    #    boolean
    """
    if type(sourceFile) is types.StringType:
	sourceFile = [sourceFile]
    if type(targetFile) is types.StringType:
	targetFile = [targetFile]
    ss = max(map( lambda x:os.stat(x).st_mtime, sourceFile))
    try:
	ts = min(map(lambda x:os.stat(x).st_mtime, targetFile))
    except:
	return True

    return ss > ts
    
#----------------------------------------------------------------------------

if __name__ == "__main__":
    m=MGIGFFMaker(sys.argv[1:])
    m.go()
