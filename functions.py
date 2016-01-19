#! /usr/bin/env python

import sys
import glob
import os.path
import cPickle as pickle
from collections import Counter

def csv_parser(infile):
	'''Parse a csv file and return index and values dictionary
	File must have a header with no first column
	The rest of the file has an index value in first column
	After index there are the same number of fields as header
	
	E.g.
	,col1,col2
	1,val1,val2
	2,val3,val4
	
	Returns a generator:
	(1,{col1:val1,col2:val2})
	(2,{col1:val3,col2:val4})
	'''
	with open(infile) as f:
		header = f.readline().strip().split(",")[1:]
		for line in f:
			line = line.strip().split(",")
			index, values = line[0], line[1:]
			assert len(values) == len(header), "Header dimensions don't match body: %s" % prot
			valuesD = dict(zip(header,values))
			yield index, valuesD

def flatten(tuple_gen,group_col=0,value_col=1):
	'''
	(a,1)(a,2)(b,1) ---> (a,(1,2))(b,1)
	
	Infile must be sorted on group column.
	'''
	group = None
	values = []
	for tup in tuple_gen:
		assert len(tup) == 2, "Tuples must be length 2"
		assert group_col != value_col, "Group and Value columns must be different"
		curr_group = tup[group_col].strip()
		curr_value = tup[value_col].strip()
		if group != curr_group:
			if not group == None: # if not first iteration
				yield group, tuple(values) # yield last group
			group = curr_group
			values = []
		values.append(curr_value)
	yield group, tuple(values)
	
def id_convert(prot_stream,mapping=None):
	'''
	For converting protein ids. Takes a stream of ids and a path to a pickled diciontary
	
	Open mapping (a pickle file).'''
	if mapping == None:
		mapping = "/project/LECA/info_files/all_prot2gene.p"
		assert os.path.exists(mapping), "No file %s" % mapping
	with open(mapping) as f:
		D = pickle.load(f)
		assert type(D) is dict, "Mapping must be pickled dictionary"
	for prot in prot_stream:
		yield D[prot]
		
def stream_2cols(infile,delim=None,header=False):
	'''
	Open a two column file and return a generator of tuples corresponding to the two columns.
	delim: optional delimiter string, default is to split line on whitespace
	header: default is False - if True, skip first line of infile
	'''
	with open(infile) as f:
		if header:
			f.next()
		for line in f:
			if delim:
				line = line.strip().split(delim)
			else:
				line = line.strip().split()
			yield line[0],line[1]

def pickle_2cols(col_stream,pickle_file):
	'''Given an imput stream of length 2 tuples, make dictionary mapping first to second
	element, and pickle the dictionary.'''
	D = {}
	count = 0
	for tup in col_stream:
		D[tup[0]] = tup[1]
		count +=1
	with open(pickle_file,'w') as f:
		pickle.dump(D,f)
	print "Pickled %i lines" % count
	
def taxon_lookup(DB='eggnog'):
	'''Return taxon lookup dictionary'''
	db_taxonPath = {'eggnog':"/project/LECA/eggNOG/info_files/eggnogv4.taxonomies.p",
	'phylomedb':"/project/LECA/eggNOG/info_files/eggnogv4.taxonomies.p",
	'orthomcl':"/project/LECA/eggNOG/info_files/eggnogv4.taxonomies.p"}
	db = DB.lower()
	assert db in db_taxonPath, "Database %s not found" % DB
	taxon_path = db_taxonD[db]
	assert os.path.exists(taxon_path), "No file: %s" % taxon_path
	with open(taxon_path) as f:
		taxonD = pickle.load(f)
	return taxonD

	
def _extract_prot(file):
	split_file = file.split("-")
	if split_file[0] == 'nan':
		return split_file[2]
	else:
		return split_file[0]
	
def _get_redundant_prots():
	file_iter = (i for i in glob.iglob("*-*-*-*.csv"))
	return [i for i,j in 
		Counter((_extract_prot(file) for file in file_iter)).iteritems() 
			if j>1]
		
def _count_sum(infile,protein):
	for prot, lineD in csv_parser(infile):
		if prot == protein:
			return float(lineD["COUNT"])
	else:
		raise Exception("Protein %s not found in %s" % (protein,infile))
	
def nonRedundant_filestream():
	'''
	Get a file stream of Claire's tables where no protein ids are repeated. There are
	repeats due to one-to-many mapping between Uniprot and Ensembl, e.g.:
		'O00241-SIRB1_HUMAN-ENSP00000371016-HUMAN.csv',
		'O00241-SIRB1_HUMAN-nan-HUMAN.csv'
	This function takes only the table with the highest "COUNT" column for the 
	target protein, i.e. the one with the fewest missing algorithms.
	'''
	file_iter = glob.iglob("*-*-*-*.csv")
	redundant_prots = _get_redundant_prots()
	for file in file_iter:
		prot = _extract_prot(file)
		if prot not in redundant_prots:
			yield file
	for p in redundant_prots:
		files = glob.iglob("*%s*.csv" % p)
		yield max(((f,_count_sum(f,p)) for f in files),key=lambda x: x[1])[0]
		
	
if __name__ == '__main__':
	infile = sys.argv[1]
	with open(infile) as f:
		line_gen = (line.strip().split('\t') for line in f)
		for l in flatten(line_gen):
			print l
			