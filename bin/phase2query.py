#
# phase2query.py
#
# First step of phase 2: take the MGI genes that don't have
# models (from phase 1), and get their sequence ids.
#
# Usage:
#    % python phase2query.py < INFILE > OUTFILE
#
# Input: file of MGI gene ids (e.g., the orphan MGI report). Any tab delimited
# file with MGI id in the first column will do.
#
# Processing: The MGI ids are first uploaded into a temp table.
# Then MGI is queried to retrieve sequence ids associated with
# the MGI ids. The sequences are restricted as follows:
# - from Genbank/EMBL/DDBJ or RefSeq only
# - must be associated with only one gene.
# - only DNA and RNA (no protein)
# - if DNA, must be < 10 kb in length
#
# Output: Tab delimited file having columns:
#	seqId, mgiId, symbol, name, type, chromosome
#
# Example:
#	M57401        MGI:102792      Mcptl   Gene    14
#	NM_008573     MGI:102792      Mcptl   Gene    14
# 

import sys
import mgiadhoc as db

def uploadIds( ids, conn, tname="_ids", colname="id", index="_ids_ix_0", chunksize=500):
    db.sql(('create temporary table %s (%s varchar(255))' % (tname,colname)), None, conn)
    for i in xrange(0,len(ids),chunksize):
	cmds = []
	for id in ids[i:i+chunksize]:
	    cmds.append("insert into %s values ('%s')" % (tname,id))
	db.sql(cmds,'ignore',conn)
    if index:
	db.sql(('create unique index %s on %s(%s)' % (index, tname, colname)), None, conn)
    return tname


def getSeqData(ids):
    # create db connection
    conn = db.connect()

    #
    # Create temp table of MGI ids
    idTbl = uploadIds(ids, conn, tname="_ids")

    #
    # Create temp table containing basic info for markers
    db.sql('''
	create temporary table _markers (
	    _marker_key int,
	    mgiid varchar(30),
	    symbol varchar(255), 
	    name varchar(255),
	    chromosome varchar(4),
	    mgitype varchar(50),
	    mcvtype varchar(255)
	    )
	''',None,conn)
    db.sql('''
	insert into _markers
	select 
	    m._marker_key, 
	    a.accid as "mgiid", 
	    m.symbol,
	    m.name,
	    m.chromosome,
	    t.name as "mgitype",
	    v.term as "mcvtype"
	from MRK_Marker m, MRK_Types t, ACC_Accession a, MRK_MCV_Cache v
	where m._Organism_key = 1
	and m._Marker_Type_key = t._Marker_Type_key
	and m._Marker_key = a._Object_key
	and a._MGIType_key = 2
	and a._LogicalDB_key = 1
	and a.prefixPart = 'MGI:'
	and a.preferred = 1
	and a.accID in (select id from %s)
	and m._Marker_key = v._Marker_key
	and v.qualifier = 'D'
	'''%idTbl, None, conn)
    db.sql('''
	create index idx0 on _markers(_marker_key)
	''', None, conn)
    #
    # Temp table of sequences (keys) that are associated with only 1 marker
    db.sql('''
	create temporary table _seq1marker (_sequence_key int)
	''', None, conn)
    db.sql('''
	insert into _seq1marker
	select _sequence_key
	from SEQ_Marker_Cache
	group by _Sequence_key
	having count(*) = 1
	''',None, conn)
    db.sql('''
	create index s1m_idx1 on _seq1marker(_sequence_key)
	''', None, conn)

    #
    #
    results = db.sql('''
	select
	    m._marker_key, m.mgiid, m.symbol, m.name, m.chromosome, m.mgitype, m.mcvtype,
	    c.accid as "seqid", s.length, t.term as "seqtype"
	from
	    _markers m,
	    SEQ_Marker_Cache c, 
	    SEQ_Sequence s,
	    VOC_Term t
	where
	    c._Marker_key = m._marker_key
	and c._Sequence_key = s._Sequence_key
	and c._LogicalDB_key in (9,27) /* (9,27,35,36,53) ??*/
	and s._Sequence_key in (select _Sequence_key from _seq1marker)
	and s._SequenceType_key = t._Term_key
	order by
	    m.mgiid, c.accid
	''', None, conn)

    #
    conn.close()
    #
    return results

def getSeqsForOrphanMarkers( infile, outfile ):
    if infile == "-":
	fd = sys.stdin
    else:
	fd = open(infile, 'r')
    lines = fd.readlines()
    fd.close()
    ids = map(lambda x:x.split()[0], lines)
    seqs = getSeqData(ids)
    if outfile == "-":
	ofd = sys.stdout
    else:
	ofd = open(outfile, 'w')
    TAB = '\t'
    NL = '\n'
    for s in seqs:
	writeit = False
	# restrict size if its a DNA seq. Also, no protein seqs.
	if s['seqtype'] == "DNA":
	    writeit = (s['length'] <= 10000)
	else:
	    writeit = (s['seqtype'] != "Polypeptide")

	if writeit:
	    ofd.write(TAB.join([
		s['seqid'],
		s['mgiid'],
		s['symbol'],
		s['name'],
		s['mcvtype'],
		s['chromosome'],
		]) + NL)
    ofd.close()


if __name__ == "__main__":
    if len(sys.argv) == 1:
	getSeqsForOrphanMarkers( '-','-' )
    elif len(sys.argv) == 2:
	getSeqsForOrphanMarkers( sys.argv[1], '-' )
    else:
	getSeqsForOrphanMarkers( sys.argv[1], sys.argv[2] )
