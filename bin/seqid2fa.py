#
# seqid2fa.py
#
# Fetches sequences from entrez by id. Input is a single-column file of ids.
# Output is a fasta-formatted file of sequences from Genbank (technically, 
# from the nucleotide database accessible via NCBI eutils).
# Usage:
#	% python seqid2fa.py INFILE OUTFILE
#
#

import sys
import types
import time
import urllib
import xml.dom.minidom

FETCHURL  = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
BATCHSIZE = 500
SLEEPTIME = 3
TOOL      = "MGI"
EMAIL     = "Joel.Richardson@jax.org"
NTRIES	  = 3

def openSequenceFetch(ids, tool, email, db='nucleotide', retmode='text', rettype='fasta',batchsize=BATCHSIZE,sleeptime=SLEEPTIME):
    for i in xrange(0, len(ids), batchsize):

	# Get the next batch of ids.
	params = urllib.urlencode(
	   {'db': db,
	    'retmode' : retmode,
	    'rettype' : rettype,
	    'id' : ",".join(ids[i:i+batchsize]),
	    'tool' : tool,
	    'email' : email
	    })

	# For each batch, try up to NTRIES time to get the sequences. Provides
	# some protection from intermittant errors from the eUtils server.
	for ntry in range(NTRIES):

	    try:
		time.sleep(sleeptime)
		fd = urllib.urlopen(FETCHURL, params)
	    except:
		continue

	    line = fd.readline()
	    if not line:
		break
	    elif not line.startswith(">"):
		if ntry == NTRIES-1:
		    # Last try. Write the error to stderr.
		    sys.stderr.write(line)
		    sys.stderr.write(fd.read())
		fd.close()
		continue

	    # success!
	    yield line
	    for line in fd:
		yield line
	    fd.close()
	    break
	else:
	    # if we exhaust the loop, there was an error
	    raise RuntimeError("Error from eutils.")

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
