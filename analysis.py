from itertools import combinations
from collections import OrderedDict, Counter
import cPickle as pickle
import pandas as pd

def load_pickle(infile):
	'''Read in pickled files'''
	with open(infile) as f:
		return pickle.load(f)
		
def calc_dists(infile,node_distsD):
	'''Read a database age file created by infer_age.serialize_dbAgeNodes and return a dictionary holding
	distances (in number of branches) between every non-redundant, non-self pairs of databases.'''
	dbAgeD = load_pickle(infile)
	dbs = sorted(dbAgeD.keys())
	dbDists = {db: {} for db in dbs}
	for db1, db2 in combinations(dbs,r=2): # corner matrix, off diagonal
		node1,node2 = dbAgeD[db1],dbAgeD[db2] # get inferred ages
		if node1 == None or node2 == None: # NOT IDEAL - if node == None make distance largest possible
			dbDists[db1][db2] = 20
			continue
		dbDists[db1][db2] = node_distsD[node1][node2] # get distance between nodes 
	return OrderedDict(sorted(dbDists.iteritems(), key=lambda x: x[0]))
	
def add_dicts(d1,d2):
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
		
def sum_dist(infile_stream,nodeDistsFile):
	'''Sum of distances'''
	dists = load_pickle(nodeDistsFile)
	is_first = True
	count = 0
	for f in infile_stream:
		if is_first:
			D = calc_dists(f,dists)
			is_first = False
		else:
			D = add_dicts(D,calc_dists(f,dists))
		count +=1
	return D, count

def node2taxonAge(infile,conversion_dictionary):
	'''Convert database node label mapping to a dictionary containing taxonomic names. Input conversion
	dictionary must map nodes to taxon labels.'''
	nodeAgeD = load_pickle(infile)
	taxonAgeD = {}
	for db in nodeAgeD:
		assert db not in taxonAgeD, "Repeated database name: %s" % db
		age = nodeAgeD[db]
		if age == None: # skip databases with nodes == None
			continue
		taxonAgeD[db] = conversion_dictionary[age]
	return taxonAgeD

def taxonAgeCount(infile_stream,conversion_file):
	'''Create distribution of taxon-based ages for each database'''
	conversion_dictionary = load_pickle(conversion_file)
	is_first = True
	for f in infile_stream:
		taxonAgeD = node2taxonAge(f,conversion_dictionary)
		if is_first:
			taxAgeCounter = {i:Counter([j]) for i,j in taxonAgeD.iteritems()}
			is_first = False
		else:
			for db,age in taxonAgeD.iteritems():
				if db in taxAgeCounter:
					taxAgeCounter[db].update([age])
				else:
					taxAgeCounter[db] = Counter([age])
	return taxAgeCounter
		