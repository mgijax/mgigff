#
# convertNCBI.py
#
# Filters NCBI file for features on one of the chromosomes in the chromosome file.
# Usage:
#  python convertNCBI.py /path/to/ncbidatafile.gff3 > output.gff3
# 
# 
# 

import sys
import os
import types
import gff3
import optparse

TAB	= '\t'
HASH	= '#'

class ConvertNCBI:
    def __init__(self, args):
	self.ncbiGffFile = args[1]
	self.id2chr = {} 	# maps chr seq id to chr, e.g., NC_000079.6 to chr 13.
	self.iid2xid= {}	# maps internal id to external, e.g., gene47 -> GeneID:10012341
	self.currentRegionId = None
	self.currentRegion = None
	
    def getGeneID(self, f):
	ids = f.attributes.get('Dbxref',[])
	if type(ids) is types.StringType:
	    ids = [ids]
	for x in ids:
	    if x.startswith('GeneID:'):
		return x[7:]
	#raise RuntimeError("No GeneID found in " + str(f))
	return None

    # The default sort of the file from NCBI is such that the features of one gene can be
    # interspersed with features from another. The mgigff process requires that all features for
    # a gene be together in an unbroken sequence. This iterator iterates over features sorted the
    # NBCI way and yields them sorted the way mgigff needs it.
    #
    def regroupingIterator(self, features):
	lastSeq = None
	groups = []
	for f in features:
	    #
	    if f.seqid != lastSeq:
		for glist, gset in groups:
		   for f2 in glist:
		       yield f2
		groups = []
	    lastSeq = f.seqid
	    #
	    p = f.attributes.get('Parent', None)
	    if not p:
		# feature is top-level
		groups.append( ([f], set([f.ID])) )
	    else:
		# feature is part of something else. 
		pid = p
		for i in range(len(groups)-1, -1, -1):
		  glist, gset = groups[i]
		  if pid in gset:
		    glist.append(f)
		    gset.add(f.ID)
		    break
		else:
		  raise RuntimeError("Orphan feature: " + str(f))
	for glist, gset in groups:
	    for f in glist:
	        yield f

    def process(self, f):
	# NCBI file has multi-level sort, with first level being by region 
	# (e.g. a chromosome or a contig)
	# A feature with 'region' in col 3 introduces a new region, and its features
	# will follow. See if it's one we want.
	if f[2] == 'region':
	  if f[0] == self.currentRegionId:
	    # A region feature within the current region. Skip.
	    return None
	  # Region with a different ID. See if it's one we care about.
	  self.currentRegionId = f[0]
	  chr = f.attributes.get('chromosome',None)
	  map = f.attributes.get('map', None)
	  genome = f.attributes.get('genome', None)
	  if f.attributes.get('strain',None) != 'C57BL/6J':
	    # skip anything not on B6
	    self.currentRegion = None
	  elif chr == 'Unknown':
	    # unplaced contig
	    self.currentRegion = f[0]
	  elif genome == 'genomic':
	    # unlocalized contig
	    self.currentRegion = chr + '|' + f[0]
	  elif genome == 'chromosome':
	    # regular ol' chromosome
	    self.currentRegion = chr
	  elif genome == 'mitochondrion':
	    # mitochondrion
	    self.currentRegion = 'MT'
	  else:
	    # something else
	    self.currentRegion = None
	  return None
	#
	# A feature in the current region
	if f[0] != self.currentRegionId:
	  raise RuntimeError('Internal error. Region id mismatch detected.' + str(f))
	if self.currentRegion is None:
	  # don't care about this region
	  return None
	#  
	if f[2] in ['match','cDNA_match']:
	  return None
	#
	f[0] = self.currentRegion
	isPseudo = f.attributes.get('pseudo',False)
	f[1] = isPseudo and 'PSEUDO' or f[2]
	fgid = self.getGeneID(f)
	fxid = f.attributes.get('transcript_id',None)
	fp  = f.attributes.get('Parent',None)
	fid = f.ID
	f.attributes.clear()
	if not fp:
	    if fgid:
		f.ID = fgid
		self.iid2xid[fid]=fgid
	    else:
		f.ID = fid
	else:
	    f.Parent = self.iid2xid.get(fp,fp)
	    if fxid and ('RNA' in f[2] or 'transcript' in f[2]):
	        f.ID = fxid
		self.iid2xid[fid]=fxid
	    else:
	        f.ID = fid

	return f

    def mainIterator(self):
	for f in gff3.iterate(self.ncbiGffFile):
	    if self.process(f):
		yield f

    def main(self):
        for f in self.regroupingIterator(self.mainIterator()):
	    sys.stdout.write(str(f))

#
if __name__ == "__main__":
    ConvertNCBI(sys.argv).main()
