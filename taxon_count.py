#! /usr/bin/env python

import sys
import pickle
from LECA import parsers,functions
from collections import Counter

## To get the get the number of each group in each NOG.
## Groups are: Archaea,Bacteria,Bikont,Bilateria,Mammalia,Metazoa,Unikont,Vertebrata 
## Obviously these groups are nested, but "Vertebrate" means a non-mammal 
## Vertebrate, Bilateria means a non-vertebrate Bilaterian, etc.
##

def get_taxonD(path_to_taxa):
	'''Load taxonD'''
	with open(path_to_taxa) as p:
		return pickle.load(p)
		
def filter_eggnog(src):
	'''Read from eggnog_parser and filter out just "nog" and "taxonID"
	Return generator of tuples'''
	### Could take advantage of other info from parser.
	return ((i["nog"],i["taxonID"]) for i in src)

def filter_phylomeDB(src):
	'''Read from phylome_parser and filter out just "gene" and "ortholog"
	Return generator of tuples'''
	return ((i["gene"],i["ortholog"]) for i in src)

def get_stream(infile,db,phylomeDB_taxonD=None):
	'''Take database file and return (category,taxonid) tuples from parsers 
	to feed to flatten()'''
	db = db.lower()
	funcsD = {'eggnog':(filter_eggnog, 
						parsers.eggnog_parser,
						{"infile":infile}),
			'phylomedb':(filter_phylomeDB, 
						parsers.phylome_parser,
						{"infile":infile,"as_taxid":True,"taxonD":phylomeDB_taxonD})}
	assert db in funcsD, 'Database %s not found' % db
	funcs = funcsD[db]
	return funcs[0](funcs[1](**funcs[2]))
	
def taxon_stream(infile,db,taxonD):
	'''Convert taxonids to taxonomic names using a the dictionary taxonD.
	Return (group,taxon) generator.'''
	src = get_stream(infile,db)
	for i,j in src:
		yield i,taxonD[j] # could add error handling?
		
def flat_stream(infile,db,taxonD):
	'''Flatten groups. Return generator of tuples (group, (taxon,taxon...))'''
	unique_categories = []
	src = taxon_stream(infile,db,taxonD)
	return (i for i in functions.flatten(src))
	
def counter_stream(infile,db,taxonD):
	'''Collect tuples from flat_stream and count occurrences of each taxon id 
	in the second element. Return generator of tuples (group, {taxon1: 2, taxon2: 39,...})'''
	for i,j in flat_stream(infile,db,taxonD):
		yield i, Counter(j)
		
def print_taxon_count(infile,db,taxonD):
	line_count = 0
	taxon_list = sorted(list(set([i for i in taxonD.itervalues()]))) # unique taxa from taxonD
	print "\t".join(["group"] + taxon_list)
	for i,j in counter_stream(infile,db,taxonD):
		counts_list = [str(j[taxon]) for taxon in taxon_list]
		assert len(counts_list) == len(taxon_list)
		print i,"\t","\t".join(counts_list)
		line_count +=1
		if line_count % 100 == 0:
			sys.stderr.write(str(line_count)+"\n")
	
		
if __name__ == '__main__':
	try:
		path_to_taxa = sys.argv[3]
	except IndexError:
		path_to_taxa = "/project/LECA/info_files/taxonomyD.p"
	infile,dbtype = sys.argv[1],sys.argv[2]
	taxonD = get_taxonD(path_to_taxa)
	print_taxon_count(infile,dbtype,taxonD)