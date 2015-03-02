#! /usr/bin/env python

import sys,os

## To get the get the number of each group in each NOG associated with human diseases
## Groups are: Archaea,Bacteria,Bikont,Bilateria,Mammalia,Metazoa,Unikont,Vertebrata 
## Obviously these groups are nested, but "Vertebrate" means a non-mammal 
## Vertebrate, Bilateria means a non-vertebrate Bilaterian, etc.
##
## awk 'NR==FNR{a[$2];next}($1 in a)' HumanNOGsDiseases.tsv NOG.members.txt | python taxon_count.py

## Minor problems:
## 	TaxIDs 240176,931890 are not in eggnogv4.taxonomies.tsv
##	No taxon level for Opisthokonts, i.e. Fungi are grouped with Dicty


def taxon_id_gen(infile):
	'''
	Parses eggnogv4.taxonomies.tsv. Returns generator of tuples: (taxonid,class),
	where class can be one of the following:
	Mammal, Vertebrate, Bilateria,Metazoa,Unikont,Bikont,Bacteria,Archaea.
	Obviously these groups are nested, but "Vertebrate" means a non-mammal 
	Vertebrate, Bilateria means a non-vertebrate Bilaterian, etc.
	'''
	L = [(19,"Mammalia"),(11,"Vertebrata"),(7,"Bilateria"),(5,"Metazoa")]
	
	Unikonts = ["Opisthokonta","Amoebozoa"]
	
	Bikonts = ["Alveolata","Euglenozoa","Fornicata","Heterolobosea","Parabasalia",
	"Rhodophyta","Stramenopiles","Viridiplantae"]

	with open(infile) as f:
		for line in f:
			line = line.strip().split("\t")
			taxonid = line[1]
			try:
				if line[3] == 'Bacteria':
					yield taxonid, "Bacteria"
				elif line[3] == 'Archaea':
					yield taxonid,"Archaea"
				elif line[4] in Bikonts:
					yield taxonid, "Bikont"
				elif line[4] in Unikonts:
					for i in L:
						try:
							if line[i[0]] == i[1]:
								yield taxonid,i[1]
								break
						except IndexError:
							continue
					else:
						yield taxonid,"Unikont"
				else:
					print line
					break
			except IndexError:
				print line
				break
				
def taxonD(infile):
	'''Create and return a dictionary from taxon_id_gen and
	the set of groups found (non-redundant list of values
	in the dictionary)'''
	#return {i:j for i,j in taxon_id_gen(infile)}
	group_set = []
	taxD = {}
	for i,j in taxon_id_gen(infile):
		taxD[i] = j
		if j in group_set:
			continue
		else:
			group_set.append(j)
	return taxD, group_set
		
	
def taxon_count_gen(handle,taxa_file):
	'''
	Returns generator yielding a tuple: (KOG, group_counts), where
	group counts is a dictionary holding the counts for each taxonomic
	group.
	
	handle: an opened file like NOG.members.txt, or similar
	taxa_file: should be eggnogv4.taxonomies.tsv, or similar
	'''
	kog = None
	taxD,group_set = taxonD(taxa_file) # taxonIDs mapped to groups
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
			
def print_taxon_count(handle,taxa_file):
	'''
	Print output from taxon_count_gen as a tab-delimited file
	
	handle: an opened file like NOG.members.txt, or similar
	taxa_file: should be eggnogv4.taxonomies.tsv, or similar
	'''
	is_first = True
	for i,j in taxon_count_gen(handle,taxa_file):
		if is_first:
			header = "\t".join(["NOG"] + sorted(j.keys()))
			yield header
			is_first = False
		yield "\t".join([i] + [str(j[taxon]) for taxon in sorted(j.keys())])
		
if __name__ == '__main__':
	path_to_taxa = "/project/LECA/eggNOG/info_files/eggnogv4.taxonomies.tsv"
	try:
		infile = sys.argv[1]
		with open(infile) as f:
			for i in print_taxon_count(f,path_to_taxa):
				print i
	except IndexError: # reading from stdin
		infile = sys.stdin
		for i in print_taxon_count(infile,path_to_taxa):
			print i
