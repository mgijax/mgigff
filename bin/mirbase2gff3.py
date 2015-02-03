
#
# mirbase2gff3.py
#
import sys
import gff3

for f in gff3.iterate( sys.stdin ):
    if f[2] == "miRNA_primary_transcript":
	if f[0].startswith('chr'):
	    f[0] = f[0][3:]
	f[1] = f[2]
	f[2] = 'gene'
	sys.stdout.write(str(f))
