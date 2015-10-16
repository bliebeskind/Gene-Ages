import pandas as pd
import cPickle as pickle
from functions import csv_parser
from itertools import combinations
from collections import OrderedDict


def load_pickle(infile):
	'''Read in pickled files'''
	with open(infile) as f:
		return pickle.load(f)
	
def all_by_all_dists(infile,node_distsD):
	'''Read a csv file of node age calls and return a generator of OrderedDicts
	that give the distances between databases for each gene. These will be tallied
	by sum_dists.'''
	parsed = csv_parser(infile)
	for prot, dbAgeD in parsed:
		dbs = sorted(dbAgeD.keys())
		dbDists = {db: {} for db in dbs}
		for db1, db2 in combinations(dbs,r=2): # corner matrix, off diagonal
			node1,node2 = dbAgeD[db1],dbAgeD[db2] # get inferred ages
			if node1 == 'None' or node2 == 'None': # If node == None make distance largest possible
				dbDists[db1][db2] = 20
				continue
			dbDists[db1][db2] = node_distsD[node1][node2] # get distance between nodes 
		yield OrderedDict(sorted(dbDists.iteritems(), key=lambda x: x[0]))
	
def add_dicts(d1,d2): # should be done with counter addition
	'''
	Given two nested dictionaries, each of depth 2, add the sums held at the second level, e.g.:
	
	d1 = {'a':{'b':1,'c':3}}
	d2 = {'a':{'b':3,'c':2}}
	add_dicts(d1,d2)
		{'a':{'b':4,'c':5}}
	'''
	out = {}
	for col in d1:
		assert col not in out, "Repeat column names: %s" % col
		out[col] = {}
		for row in d1[col]:
			assert row not in out[col], "Repeat row names: %s" % row
			try:
				d2_value = d2[col][row]
			except KeyError, e:
				raise Exception("Couldn't find value %s" % (e))
			out[col][row] = d1[col][row] + d2_value
	return out
		
def sum_dist(infile,nodeDistsFile):
	'''Sum of patristic distances between database age calls'''
	dists = load_pickle(nodeDistsFile)
	is_first = True
	count = 0
	for distD in all_by_all_dists(infile,dists):
		if is_first:
			D = distD
			is_first = False
		else:
			D = add_dicts(D,distD)
		count += 1
		if count % 100 == 0:
			print count
	return pd.DataFrame(D), count
	
def avg_dist(infile,nodeDistsFile):
	'''Return a pandas DataFrame holding average patristic distances between 
	database age calls'''
	sumDistDF,count = sum_dist(infile,nodeDistsFile)
	return sumDistDF.applymap(lambda x: x/float(count))
