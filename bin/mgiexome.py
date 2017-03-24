#
# mgiexome.py
#
# Post process mgi gff file to produce mgi exome file.
#
# The exome file contains only exon features. These exons 
# are derived by "unioning" (per gene) the exons in the MGI 
# GFF3 file. Each exon is given an ID and is tagged with the
# mgi id and symbol of its gene. 
#
# The previous statement is almost correct. This program also 
# outputs pseudogenic_exons and match_parts.
#
# Genes associated with multiple models union the exons across
# all models. Models associated with multiple genes contribute
# to the exons for all genes. 
#
# 2013-03-18 jer - Handle pseudogenes (look for "pseudogenic_exons").
# Special handling for some NCBI pseudogenes, which have no internal structure, just start/end.
# The search for "pseudogenic_exon" doesn't find anything for these. 
# 

import sys
import gff3
import types
import logging

#
def compress(id,exons):
    """
    Given a set of possibly overlapping exons associated with an id,
    produces and returns a set of nonoverlapping exons by "unioning"
    the input, i.e., by merging overlapping exons into larger
    exons that span them. The output exons are given ID attributes
    constructed from the input id plus an appended exon counter.
    E.g., the exons for "MGI:12345" will have IDs "MGI:12345_001",
    "MGI:12345_002", ...
    """
    exons.sort( None, lambda x : x.start )
    oexons = []
    laste = None
    for i,e in enumerate(exons):
	if i==0 or e.start > laste.end:
	    # Start a new exon
	    oexons.append(e)
	    e.source = "MGI"
	    e.attributes = {"ID":"%s_%03d"%(id,len(oexons))}
	    laste = e
	else:
	    # New exon overlaps current. Merge by extending coords.
	    # Dont have to worry about start coord because of sort order.
	    # So just extend the end coord.
	    laste.end = max(laste.end, e.end)
    return oexons

#
def flush( mgi2exons, mgi2symbol, ofile=sys.stdout ):
    """
    For each gene id, compress and output its exons.
    Tag exons with gene's id and symbol.
    Then clear the buffer.
    """
    for id,exons in mgi2exons.iteritems():
	cexons = compress(id,exons)
	for e in cexons:
	    e.mgiId = id
	    e.mgiSymbol = mgi2symbol.get(id, '???') # FIXME
	    ofile.write( str(e) )
    mgi2exons.clear()

#
def main(infile,outfile=None):
    logging.info("Generating exome file.")

    #
    # Set up output and write GFF3 header.
    #
    if type(outfile) is types.StringType:
	ofile = open(outfile, 'w')
    else:
	ofile = sys.stdout
    ofile.write(gff3.HEADER)

    mgi2symbol = {}

    #
    # Iterate by groups, separated by "###". 
    #
    for grp in gff3.iterate(infile,True): 
	mgi2exons = {}
	currMgiId = None
	currMgiSymbol = None
	mgipseudos = []
	for f in grp:
	    if f.source == "MGI":
		mgiid = f.ID[4:]
		mgi2symbol[mgiid] = f.Name
		if f.type == "pseudogene":
		    mgipseudos.append(f)
	    elif f.type in ["exon","pseudogenic_exon","match-part"] and f.attributes.has_key('Dbxref'):
		#
		# Exon with a Dbxref attribute. Add this exon to all MGI ids
		#
		dbx = f.Dbxref
		if type(dbx) is types.StringType:
		    dbx = [dbx]
		for d in dbx:
		    #
		    # Accept both "MGI:MGI:12345" and "MGI:12345"
		    # Accumulate exon features for each MGI id seen.
		    # 
		    if d.startswith("MGI:"):
			mgiid = d
			if mgiid.startswith("MGI:MGI:"):
			    mgiid = mgiid[4:]
			g = gff3.Feature(f)
			mgi2exons.setdefault(mgiid,[]).append(g)
	# end: for f 

	for f in mgipseudos:
	    mgiid = f.ID[4:]
	    if mgiid not in mgi2exons:
		e = gff3.Feature(f)
		mgi2exons[mgiid] = [e]

	flush(mgi2exons, mgi2symbol, ofile)

    # end: for grp 
    ofile.close()

#
if __name__ == "__main__":
    if len(sys.argv) == 2:
	main(sys.argv[1])
    else:
	main("-")
