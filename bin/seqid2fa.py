#
# seqid2fa.py
#
# Fetches sequences from entrez by id.
#

import sys
import types
import time
import urllib
import xml.dom.minidom

FETCHURL  = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
BATCHSIZE =1000
SLEEPTIME =1
TOOL      = "MGI"
EMAIL     = "Joel.Richardson@jax.org"

def openSequenceFetch(ids, tool, email, db='nucleotide', retmode='text', rettype='fasta',batchsize=BATCHSIZE,sleeptime=SLEEPTIME):
    lasttime = 0
    for i in xrange(0, len(ids), batchsize):
	params = urllib.urlencode(
	   {'db': db,
	    'retmode' : retmode,
	    'rettype' : rettype,
	    'id' : ",".join(ids[i:i+batchsize]),
	    'tool' : tool,
	    'email' : email
	    })
	#
	t = time.time()
	st = max(0,sleeptime - (t - lasttime))
	if st > 0:
	    time.sleep(st)
	lasttime = t+st
	#
	f = urllib.urlopen(FETCHURL, params)
	for line in f:
	    yield line
	f.close()

def fetchSequences(
	infile,
	outfile,
	tool,
	email,
	db='nucleotide',
	retmode='text',
	rettype='fasta',
        batchsize=BATCHSIZE,
	sleeptime=SLEEPTIME):

    if infile == "-" or infile is sys.stdin:
	fd = sys.stdin
	infile = "<stdin>"
    else:
	fd = open(infile, 'r')
    lines = fd.readlines()
    fd.close()
    #
    ids = map(lambda l:l.split()[0], lines)
    if outfile == "-" or outfile is sys.stdout:
	fd = sys.stdout
    else:
	fd = open(outfile, 'w')
    for line in openSequenceFetch(ids,tool,email,db,retmode,rettype,batchsize,sleeptime):
	fd.write(line)
    if fd is not sys.stdout:
	fd.close()

if __name__ == "__main__":
    if len(sys.argv) == 1:
	fetchSequences('-', '-', TOOL, EMAIL)
    elif len(sys.argv) == 2:
	fetchSequences(sys.argv[1], '-', TOOL, EMAIL)
    elif len(sys.argv) == 3:
	fetchSequences(sys.argv[1], sys.argv[2], TOOL, EMAIL)
