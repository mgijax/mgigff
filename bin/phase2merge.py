#
# phase2merge.py
#
# Attaches MGI gene information to Blat alignments.
# 
# Usage:
#    % python phase2merge.py SEQIDFILE PSLFILE > OUTPUTPSLFILE
#

import sys
import psl

TAB	= '\t'
NL	= '\n'
PIPE	= '|'
DOT	= '.'

def readSeqIdFile(fname):
    '''
    Reads the file containing seqids and MGI gene info. Returns
    dict from seqid to the info.
    '''
    idx = {}
    fd = open(fname, 'r')
    for line in fd:
	seqid, mgiid, symbol, name, type, chr = line.strip().split(TAB)
	idx[seqid] = [mgiid,seqid,type,symbol,name,"chr%s"%chr]
    return idx

def merge(fname, idx):
    for a in psl.iterate(fname):
	seqid = a.qName.rsplit(DOT,1)[0]
	mgistuff = idx.get(seqid,None)
	if mgistuff:
	    a.qName = PIPE.join(mgistuff+[seqid])
	    sys.stdout.write(str(a))
	else:
	    sys.stderr.write('Skipping id %s. Data = %s\n' % (seqid, str(a)))

if len(sys.argv) != 3:
    sys.stderr.write('Usage:\n\tpython %s SEQIDFILE PSLFILE > OUTPUTPSLFILE\n'%sys.argv[0])
    sys.exit(-1)

idx = readSeqIdFile(sys.argv[1])
merge(sys.argv[2], idx)
