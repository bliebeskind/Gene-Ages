#! /usr/bin/env python

## For download taxonomic information from NCBI for the species in eggnog
## To produce the file eggnogv4.taxonomies.txt, I used the following:
## grep -v "#" eggnogv4.species.txt | python fetch_taxonomies.py
##
## Can also be used with a smaller file either piped, as above, or as the first argument
## as long as the file has the species name in the first tab-delimited column and the
## NCBI taxon id in the second.

from Bio import Entrez
import sys

Entrez.email = "bliebeskind@austin.utexas.edu"

try:
	infile = sys.argv[1]
	with open(infile) as f:
		for line in f:
			line = line.strip().split("\t")
			taxon,taxID = line[0],line[1]
			sys.stderr.write("Fetching %s\n" % taxon)
			handle = Entrez.efetch(db="taxonomy",id=taxID)
			record = Entrez.read(handle)
			print ''.join([taxon,"\t",taxID,"\t",str(record[0]),'\n'])
except IndexError:
	infile = sys.stdin
	for line in infile:
		line = line.strip().split("\t")
		taxon,taxID = line[0],line[1]
		sys.stderr.write("Fetching %s\n" % taxon)
		handle = Entrez.efetch(db="taxonomy",id=taxID)
		record = Entrez.read(handle)
		print ''.join([taxon,"\t",taxID,"\t",str(record[0]),'\n'])

