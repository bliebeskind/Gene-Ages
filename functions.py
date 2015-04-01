#! /usr/bin/env python

import sys,pickle

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
	with open(mapping) as f:
		D = pickle.load(f)
		assert type(D) is dict, "Mapping must be pickled dictionary"
	for prot in prot_stream:
		yield D[prot]
		
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
	
	
if __name__ == '__main__':
	infile = sys.argv[1]
	with open(infile) as f:
		line_gen = (line.strip().split('\t') for line in f)
		for l in flatten(line_gen):
			print l
			