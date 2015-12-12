import numpy as np
import pandas as pd
import cPickle as pickle
from functions import csv_parser
from itertools import combinations_with_replacement
from collections import OrderedDict,Counter

def load_pickle(infile):
	'''Read in pickled files'''
	with open(infile) as f:
		return pickle.load(f)
	
def all_by_all_dists(infile,node_distsD,nonePenalty=None):
	'''Read a csv file of node age calls and return a generator of OrderedDicts
	that give the distances between databases for each gene. These will be tallied
	by sum_dists.'''
	parsed = csv_parser(infile)
	for prot, dbAgeD in parsed:
		dbs = sorted(dbAgeD.keys())
		dbDists = Counter({db: Counter({}) for db in dbs})
		for db1, db2 in combinations_with_replacement(dbs,r=2): # corner matrix, off diagonal
			node1,node2 = dbAgeD[db1],dbAgeD[db2] # get inferred ages
			if node1 == 'None' or node2 == 'None': # If node == None make distance largest possible
				if nonePenalty == None:
					continue
				else:
					dbDists[db1][db2] = nonePenalty
					continue
			dbDists[db1][db2] = node_distsD[node1][node2] # get distance between nodes
		yield dbDists
		
def sum_dist(infile,nodeDistsFile,nonePenalty=None):
	'''Sum of patristic distances between database age calls'''
	dists = load_pickle(nodeDistsFile)
	is_first = True
	count = 0
	for distD in all_by_all_dists(infile,dists,nonePenalty):
		if is_first:
			D = distD
			dbCounts = {db:0 for db in distD.keys()}
			is_first = False
		else:
			D = D + distD # add Counter objects
		for db in distD:
			if distD[db] != Counter(): # skip dbs missing gene
				dbCounts[db] += 1
		count += 1
		if count % 100 == 0:
			print count
	return pd.DataFrame(D).fillna(np.nan), dbCounts # replace None with NaN
	
def avg_dist(infile,nodeDistsFile,nonePenalty=None):
	'''
	Return a pandas DataFrame holding average patristic distances between 
	database age calls. 
	
	Missing data is not penalized by default, but can be by setting nonePenalty
	to the number of branches (patristic distance) to penalize a database for 
	missing data. When nonePenalty is None, the average for each pair of databases
	is calculated only on genes found in both databases.
	'''
	sumDistDF,counts = sum_dist(infile,nodeDistsFile)
	return pd.DataFrame(
		{col:sumDistDF[col].map(lambda x: x/float(counts[col]))
		for col in sumDistDF.columns})
