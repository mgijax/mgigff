#
# eutils.py
#
# Library for using NCBI's eUtils tools to query/fetch
# objects from Entrez.
#

import sys
import time
import urllib
import xml.dom.minidom

BASEURL   = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
FETCHURL  = BASEURL + 'efetch.fcgi'
SEARCHURL = BASEURL + 'esearch.fcgi'
POSTURL   = BASEURL + 'epost.fcgi'
INFOURL   = BASEURL + 'einfo.fcgi'
SUMMARYURL= BASEURL + 'esummary.fcgi'

BATCHSIZE=1000
SLEEPTIME=1

#
# Iterator for fetching sequences (in fasta format) from NCBI.
# Uses eUtils.
#
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
	sys.stderr.write(FETCHURL + '\n' + str(params) + '\n')
	f = urllib.urlopen(FETCHURL, params)
	for line in f:
	    yield line
	f.close()

if __name__ == "__main__":
    if len(sys.argv) == 1:
	lines = sys.stdin.readlines()
	ids=map(lambda l:l.split()[0], lines)
    else:
	ids = sys.argv[1:]
    for line in openSequenceFetch(ids, 'MGI', 'Joel.Richardson@jax.org'):
	print line,

