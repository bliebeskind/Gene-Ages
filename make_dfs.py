### Two functions used to compare the distributions over functional categories that are
### present in each of the files created by funccat_counts. 
### 
### These distributions can be compared to each other or to the distribution present in
### the database as a whole. This latter gives functional enrichment for type of 
### taxonomic distribution, e.g. NOGs found across eukaryotes but not in prokaryotes 
### (LECA-specific NOGs).

import pandas as pd

def load_counts(infile_list):
	'''
	Takes a list of files like those produced by "funccat_counts.sh". This will
	return a pandas DataFrame object where each column corresponds to a certain
	taxonomic distribution, and each row is an annotation type. These distributions
	over annotation types are normalized to the total number of NOGs in each original
	file.
	'''
	D = {}
	for f in infile_list:
		group = f.split(".")[0]
		assert group not in D, "Group %s found twice" % group
		D[group] = {}
		with open(f) as f:
			for i in f:
				line = i.strip().split()
				category = line[1].strip()
				if category =='S' or category == 'R': # skip uncharacterized NOGs
					continue
				D[group][category] = line[0]
	df = pd.DataFrame(D).astype("float").fillna(value=0)
	return df/df.sum()
	
def load_background_dist(infile):
	'''Load a file like one created by "funccat_counts.sh". Return a normalized
	pandas Series.'''
	D = {}
	with open(infile) as f:
		for i in f:
			line = i.strip().split()
			category = line[1].strip()
			if category =='S' or category == 'R': # skip uncharacterized NOGs
				continue
			D[category] = line[0]
	S = pd.Series(D).astype("float")
	return S/S.sum()