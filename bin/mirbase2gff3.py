
#
# mirbase2gff3.py
#
# The GFF3 file from mirBASE has two feature types: miRNA_primary_transcript and miRNA.
# miRNAs are mature products derived from primary transcripts. They have a "Derives_from" 
# attribute that points to their transcript.
# 
# The mapping to MGI's GFF3 representation is as follows:
#	- Each primary transcript generates a gene feature and a transcript feature.
#	  The transcript points to the gene via the Parent attribute.
#	  The gene gets the ID. The transcript gets its own ID (created from the original).
#	- Each mature miRNA generates an exon.
#	  The Parent of the exon is the transcript designated in the Derives_from attribute.
#
import sys
import gff3

seenIds = set()
tCounts = {}
tMap = {}
for f in gff3.iterate( sys.stdin ):
    i = f.ID
    if i in seenIds:
	sys.stderr.write('DUPLICATE ID (skipped): ')
	sys.stderr.write(str(f))
	sys.stderr.write('\n')
	continue
    else:
	seenIds.add(i)

    if f[0].startswith('chr'):
	f[0] = f[0][3:]
    f[1] = "miRNA"

    if f[2] == "miRNA_primary_transcript":
	#
	# Gene feature
	#
	f[2] = 'gene'	# 
	sys.stdout.write(str(f))
	#
	# Transcript feature
	#
	f[2] = 'transcript'
	f.Parent = i
	nt = tCounts.setdefault(i,0) + 1
	tCounts[i] = nt
	f.ID = i + "_T"
	tMap[i] = f.ID
	sys.stdout.write(str(f))
    elif f[2] == "miRNA":
	#
	# Exon feature
	#
	f[2] = 'exon'
	df = f.Derives_from
	del f.attributes['Derives_from']
	f.ID = i
	f.Parent = tMap[df]
	sys.stdout.write(str(f))

