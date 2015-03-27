#! /usr/bin/env python

import sys,os
import pickle

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
	
def pickle_taxonD(infile,pickle_file):
	D,group_set = taxonD(infile)
	with open(pickle_file,'w') as f:
		pickle.dump(D,f)
	
if __name__ == '__main__':
	path_to_taxa = "/project/LECA/eggNOG/info_files/eggnogv4.taxonomies.tsv"
	pickle_file = sys.argv[1]
	pickle_taxonD(path_to_taxa,pickle_file)