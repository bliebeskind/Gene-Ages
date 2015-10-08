import sys
import numpy as np
import pandas as pd
import cPickle as pickle
import multiprocessing as mult
from itertools import combinations
from collections import OrderedDict, Counter

# shared funcs

def load_pickle(infile):
	'''Read in pickled files'''
	with open(infile) as f:
		return pickle.load(f)

def _paraStream(func,instream):
	chunk = mult.cpu_count() - 1
	pool = mult.Pool(processes=chunk)
	return pool.imap_unordered(func,instream,chunk)

# Sum distances
		
def all_by_all_dists(infile,node_distsD):
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
	
def one_by_all_dists(infile,node_distsD,database):
	'''Read in database age file created by infer_age.serialize_dbAgeNodes and calculate the average 
	distance between an input database and all other databases. Return a tuple of the protein 
	(from infile) and the distance.'''
	dbAgeD = load_pickle(infile)
	prot = infile.split(".")[0]
	assert database in dbAgeD, "Database %s not found" % database
	otherDBs = dbAgeD.keys()
	otherDBs.remove(database)
	totalDist = 0
	totalDBs = 0.0
	focalDBAge = dbAgeD[database]
	if focalDBAge == None:
		return prot, "None"
	for i in otherDBs:
		node = dbAgeD[i]
		if node == None:
			continue
		else:
			dist = node_distsD[node][focalDBAge]
			totalDist += dist
			totalDBs += 1
	try:
		avgDist = totalDist/totalDBs
	except ZeroDivisionError:
		return prot, 0
	return prot, avgDist
	
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
		
def sum_dist(infile_stream,nodeDistsFile):
	'''Sum of distances'''
	dists = load_pickle(nodeDistsFile)
	is_first = True
	count = 0
	for f in infile_stream:
		if is_first:
			D = all_by_all_dists(f,dists)
			is_first = False
		else:
			D = add_dicts(D,all_by_all_dists(f,dists))
		count +=1
	return D, count

# Binned ages	

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

def _deTuppler(args): # crappy workaround to pass pickleable function with only 1 argument to Pool
	return node2taxonAge(*args)

def taxonAgeCount(infile_stream,conversion_file):
	'''Create distribution of taxon-based ages for each database'''
	conversion_dictionary = load_pickle(conversion_file)
	is_first = True
	ingen = ((f,conversion_dictionary) for f in infile_stream)
	outgen = _paraStream(_deTuppler,ingen) # create parallelized generator
	for taxonAgeD in outgen:
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
	
## Get ages

def ages_to_tsv(infile_stream,as_taxa=False,conversion_dictionary=None):
	'''Print pickled age dictionaries to a csv'''
	is_first = True
	if conversion_dictionary and None not in conversion_dictionary:
		conversion_dictionary[None] = "None"
	for f in infile_stream:
		protein = f.split(".")[0]
		ages = load_pickle(f)
		if is_first: # print header with leading tab character
			dbs = sorted([''] + ages.keys())
			print "\t".join(dbs)
			is_first = False
		else:
			assert sorted(ages.keys()) == dbs, "Databases don't match: %s" % protein
			if as_taxa:
				try:
					print "\t".join([protein] + [conversion_dictionary[ages[i]] for i in dbs])
				except TypeError:
					raise Exception("conversion dictionary must be supplied if as_taxa is True")
			else:
				print "\t".join([protein] + [str(ages[i]) for i in dbs])
				
## Polarization stat

def _within(L1,L2,node_distsD,dbAgeD):
	'''Calculate the average within group dists for the two groups L1 and L2'''
	dists = []
	for l in [L1,L2]:
		for i,j in combinations(l,r=2):
			n1,n2 = dbAgeD[i],dbAgeD[j] # convert dbs to nodes
			if n1 == None or n2 == None or n1 == 'None' or n2 == 'None':
				continue
			dists.append(node_distsD[n1][n2]) # get dist between nodes
	return sum(dists)/float(len(dists))

	
def _between(L1,L2,node_distsD,dbAgeD):
	'''Calculate the between groups dists for L1 and L2'''
	dists = []
	for i in L1:
		n1 = dbAgeD[i] # convert dbs to nodes
		if n1 == None or n1 == 'None':
			continue
		for j in L2:
			n2 = dbAgeD[j] # convert dbs to nodes
			if n2 == None or n2 == 'None':
				continue
			dists.append(node_distsD[n1][n2]) # get dist between nodes
	return sum(dists)/float(len(dists))

def _checkNames(L,D):
	for name in L:
		try:
			assert name in D
		except AssertionError:
			raise Exception("%s not found" % name)

def polarization(dbAgeD,gene,node_distsD,class1,class2,polMetric=False):
	'''
	Calculate the polarization score statistic or metricfor a give set of age calls. The
	score is a ratio of the within/between patristic distances for the two classes. The
	metric is the absolute value of the difference (between - within) and can be used
	as a penalty.
	
	Input should
	be a dictionary mapping each database to its inferred age node. These nodes must match the
	entries in the node_distsD, which gives the patristic distance between each pairs of nodes.
	
	Suggested classes to compare are:
	class1 = ["InParanoid","InParanoidCore","OMA_Groups","OMA_Pairs","PANTHER8_LDO","RSD"]
	class2 = ["Orthoinspector","Hieranoid_2","EnsemblCompara_v2","PANTHER8_all","Metaphors","PhylomeDB"]
	'''
	_checkNames(class1,dbAgeD)
	_checkNames(class2,dbAgeD)
	try:
		wInDists = _within(class1,class2,node_distsD,dbAgeD)
		betweenDists = _between(class1,class2,node_distsD,dbAgeD)
	except ZeroDivisionError:
		sys.stderr.write("%s: too many None's to calculate\n" % gene)
		return
	if polMetric:
		score = abs(betweenDists - wInDists)
	else:
		score = wInDists/betweenDists # won't raise zerodivision error because these are type numpy.float64
		if np.isnan(score): # but will return nan
			if wInDists == 0.0:
				return "\t".join([gene,"1"])
			else:
				raise Exception("%s: between dist is zero, but within is >0\n" % gene)
	return ",".join([gene, str(score)])
	
	
### LDO analysis - see whether InParanoid and OMA are collapsing co-orthologous groups

def _findLDObreak(orthoAges,odb,ydb):
	for gene in orthoAges.index:
		youngAge,oldAge = orthoAges.ix[gene,ydb], orthoAges.ix[gene,odb]
		if np.isnan(youngAge) or np.isnan(oldAge): # skip missing values
			continue
		elif youngAge >= oldAge: # skip genes where 'young' db calls older value
			continue # watch for cases where these all fail
		else:
			try:
				ldos = len(orthoAges.loc[orthoAges[ydb] == oldAge]) # num rows where young db found the older age
			except ValueError:
				print ydb,oldAge,gene
				raise
			yield gene, bool(ldos) # 1/0 --> True/False - i.e. LDOs were found or not
				
	
## Would like this to spit out a tuple of gene, whether an LDO split was detected (T/F), for which database pair
def _LDOcomp(orthoAges,oldGroup,youngGroup,binnedConversion):
	'''Do the analysis for a single orthogroup'''
	for odb in oldGroup:
		for ydb in youngGroup:
			orthoAgesTrimmed = orthoAges[[odb,ydb]] # trim DF - now just orthos and two DBs.
			if binnedConversion:
				func = lambda x: binnedConversion[x]
				orthoAgesTrimmed = orthoAgesTrimmed.applymap(func)
			for gene, value in _findLDObreak(orthoAgesTrimmed,odb,ydb):
				yield gene, value, odb, ydb


def run_LDOcomp(coOrthoFile,ageFile,oldGroup,youngGroup,binnedConversion=None):
	'''
	coOrthos is a file like coOrthoGroups.txt
	ageFile is a csv file holding either node ages or categorical. If categorical, set binnedConversion = True
	oldGroup and youngGroup are lists of databases for comparison. Must match headers in ageFile
	'''
	if binnedConversion:
		binnedConversion = {'Cellular_organisms':7,'Euk_Archaea':6,'Eukaryota':5,'Opisthokonta':4,'Eumetazoa':3,'Vertebrata':2,'Mammalia':1,np.nan:np.nan}
	ages = pd.read_csv(ageFile,index_col=0,na_values=["None"])
	outD = {} # {gene : {[oldDB, younDB]:True/False,...}}
	comps = 0
	with open(coOrthoFile) as f:
		for line in f:
			orthos = [o for o in line.strip().split(",") if o in ages.index]
			orthoAges = ages.loc[orthos] # trim DF
			orthoAges.dropna(how='all',inplace=True)
			if len(orthoAges.index) <= 1: # so must check that more than one gene found after drop
				continue
			for gene,value,odb,ydb in _LDOcomp(orthoAges,oldGroup,youngGroup,binnedConversion):
				if gene in outD:
					dbs = (odb,ydb)
					if dbs in outD[gene]: # only care that it's True once
						if value == True and outD[gene][dbs] == False:
							outD[gene][dbs] = value
					else:
						outD[gene][dbs] = value
				else:
					outD[gene] = {(odb,ydb):value}
				comps +=1
				if comps % 100 == 0:
					print comps
	return outD
	
def percLDOs(resultD):
	'''Return a Dataframe giving, for each gene the percent of pairwise 
	DB comparisons for which an LDO break was found.'''
	calc = lambda d: float(len([i for i in d.itervalues() if i == True]))/len(d)
	D = {}
	for gene in resultD:
		assert gene not in D, "Gene: %s found twice" % gene
		D[gene] = calc(resultD[gene])
	return D
