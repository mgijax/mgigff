#
# convertNCBI.py
#
# Filters NCBI file for features on one of the chromosomes in the chromosome file.
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
	self.chrConvFile = os.path.join(os.path.dirname(self.ncbiGffFile), 'chr_accessions_GRCm38.p2')
	self.id2chr = {} 	# maps chr seq id to chr, e.g., NC_000079.6 to chr 13.
	self.iid2xid= {}	# maps internal id to external, e.g., gene47 -> GeneID:10012341
	
    def loadChrConvFile(self, fname):
	'''
	Loads the chromosome lookup table. This is a tab delimited
	'''
	ix = {}
	fd = open(fname, 'r')
	for l in fd:
	    if l.startswith(HASH):
		continue
	    flds = l.split(TAB)
	    ix[flds[1]] = flds[0]
	fd.close()
	return ix
    
    def getGeneID(self, f):
	ids = f.attributes.get('Dbxref',[])
	if type(ids) is types.StringType:
	    ids = [ids]
	for x in ids:
	    if x.startswith('GeneID:'):
		return x[7:]
	#raise RuntimeError("No GeneID found in " + str(f))
	return None

    def process(self, f):
	chr = self.id2chr.get(f[0], None)
	if chr is None or f[2] in ['region','cDNA_match']:
	    return None
	f[0] = chr
	isPseudo = f.attributes.get('pseudo',False)
	f[1] = isPseudo and 'PSEUDO' or f[2]
	fgid = self.getGeneID(f)
	fp  = f.attributes.get('Parent',None)
	f.attributes.clear()
	if not fp:
	    if fgid:
		f.ID = fgid
		self.iid2xid[fid]=fgid
	    else:
		f.ID = fid
	else:
	    f.ID = fid
	    f.Parent = self.iid2xid.get(fp,fp)

	return f

    def main(self):
	self.id2chr = self.loadChrConvFile(self.chrConvFile)
	for f in gff3.iterate(self.ncbiGffFile):
	    if self.process(f):
		sys.stdout.write(str(f))

#
if __name__ == "__main__":
    ConvertNCBI(sys.argv).main()
