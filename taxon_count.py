#! /usr/bin/env python

import sys
import pickle

### BROKEN - taxD and groupset need to come from a pickle file made by taxonomies.py

## To get the get the number of each group in each NOG associated with human diseases
## Groups are: Archaea,Bacteria,Bikont,Bilateria,Mammalia,Metazoa,Unikont,Vertebrata 
## Obviously these groups are nested, but "Vertebrate" means a non-mammal 
## Vertebrate, Bilateria means a non-vertebrate Bilaterian, etc.
##
## awk 'NR==FNR{a[$2];next}($1 in a)' HumanNOGsDiseases.tsv NOG.members.txt | python taxon_count.py

## Minor problems:
## 	TaxIDs 240176,931890 are not in eggnogv4.taxonomies.tsv
##	No taxon level for Opisthokonts, i.e. Fungi are grouped with Dicty
		
	
def taxon_count_gen(handle,taxD,group_set):
	'''
	Returns generator yielding a tuple: (KOG, group_counts), where
	group counts is a dictionary holding the counts for each taxonomic
	group.
	
	handle: an opened file like NOG.members.txt, or similar
	taxa_file: should be eggnogv4.taxonomies.tsv, or similar
	'''
	kog = None
	group_counts = {t:0 for t in group_set}
	curr_taxIDs = []
	for line in handle:
		if line.startswith("#"): # skip header
			continue
		line = line.strip().split("\t")
		curr_kog = line[0].strip()
		if curr_kog != kog: # new kog, yield and purge
			if not kog == None: # if not first iteration
				yield kog, group_counts
			kog = curr_kog
			group_counts = {g:0 for g in group_set}
			curr_taxIDs = []
		taxID = line[1].split(".")[0]
		if taxID in curr_taxIDs: # ignore paralogs
			continue # taxon already found
		curr_taxIDs.append(taxID)
		try:
			group = taxD[taxID]
			group_counts[group] += 1
		except KeyError: # TaxIDs 240176,931890 are not in eggnogv4.taxonomies.tsv
			continue
	yield kog, group_counts
			
def print_taxon_count(handle,taxD,group_set):
	'''
	Print output from taxon_count_gen as a tab-delimited file
	
	handle: an opened file like NOG.members.txt, or similar
	taxD:
	group_set:
	'''
	is_first = True
	for i,j in taxon_count_gen(handle,taxD,group_set):
		if is_first:
			header = "\t".join(["NOG"] + sorted(j.keys()))
			yield header
			is_first = False
		yield "\t".join([i] + [str(j[taxon]) for taxon in sorted(j.keys())])
		
if __name__ == '__main__':
	path_to_taxa = "/project/LECA/eggNOG/info_files/eggnogv4.taxonomies.p"
	with open(path_to_taxa) as pickleFile:
		taxonD,group_set = pickle.load(pickleFile)
	try:
		infile = sys.argv[1] # an opened file like NOG.members.txt
		with open(infile) as f:
			for i in print_taxon_count(f,taxonD,group_set):
				print i
	except IndexError: # reading from stdin
		infile = sys.stdin
		for i in print_taxon_count(infile,taxonD,group_set):
			print i
