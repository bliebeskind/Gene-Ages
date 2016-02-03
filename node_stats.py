import sys
import numpy as np
import pandas as pd
from itertools import combinations

## Consistency

def nodeError(ages,node_distsD):
	'''Calculate the node error score for a single parsed row of age calls
	ages: a dictionary mapping each database to it's node age call
	node_distsD: a dictionary giving the distances between nodes'''
	length = 0.0
	totalDists = 0.0
	ageCalls = [a for a in ages.itervalues()]
	for i,j in combinations(ageCalls,r=2):
		if i != 'None' and j != 'None':
			assert i in node_distsD, "Node %s not found in supplied dictionary" % i
			assert j in node_distsD, "Node %s not found in supplied dictionary" % j
			totalDists += int(node_distsD[i][j])
			length += 1
	if length == 0:
		return None
	else:
		return totalDists/length

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

def bimodality(dbAgeD,gene,node_distsD,class1,class2,polMetric=True):
	'''
	Calculate the bimodality score statistic or metric for a give set of age calls. The
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
		return ""
	if polMetric:
		score = betweenDists - wInDists
	else:
		score = wInDists/betweenDists # won't raise zerodivision error because these are type numpy.float64
		if np.isnan(score): # but will return nan
			if wInDists == 0.0:
				return 1
			else:
				raise Exception("%s: between dist is zero, but within is >0\n" % gene)
	return score
	
	
### LDO analysis - see whether DBs are over-splitting co-orthologous groups

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
				
def _get_falsePos(infile,dbs):
	'''Read a false positive (lossStats) file, and return list of genes where
	reference algorithms used for LDO analysis are flagged as being a false positive'''
	genes = []
	with open(infile) as f:
		for line in f:
			line = line.strip().split(",")
			outliers = line[3].split()
			if len(set(outliers).intersection(set(dbs))) > 0: # reference LDO algorithm has a false positive flag
				genes.append(line[0])
	return genes

def run_LDOcomp(coOrthoFile,ageFile,oldGroup,youngGroup,falsePos,binnedConversion=False):
	'''
	coOrthos is a file like coOrthoGroups.txt
	ageFile is a csv file holding either node ages or categorical. If categorical, set binnedConversion = True
	oldGroup and youngGroup are lists of databases for comparison. Must match headers in ageFile
	falsePos is a file output by loss_stats
	'''
	if binnedConversion:
		binnedConversion = {'Cellular_organisms':7,'Euk_Archaea':6,'Eukaryota':5,'Opisthokonta':4,'Eumetazoa':3,'Vertebrata':2,'Mammalia':1,np.nan:np.nan}
	ages = pd.read_csv(ageFile,index_col=0,na_values=["None"])
	outD = {} # {gene : {[oldDB, younDB]:True/False,...}}
	comps = 0
	falsePosGenes = _get_falsePos(falsePos,oldGroup)
	with open(coOrthoFile) as f:
		for line in f:
			orthos = [o for o in line.strip().split(",") if o in ages.index]
			orthoAges = ages.loc[orthos] # trim DF
			orthoAges.dropna(how='all',inplace=True)
			if len(orthoAges.index) <= 1: # so must check that more than one gene found after drop
				continue
			for gene,value,odb,ydb in _LDOcomp(orthoAges,oldGroup,youngGroup,binnedConversion):
				if gene in falsePosGenes:
					continue
				if gene in outD:
					if ydb in outD[gene]: # only care that it's True once
						if value == True and outD[gene][ydb] == False:
							outD[gene][ydb] = value
					else:
						outD[gene][ydb] = value
				else:
					outD[gene] = {ydb:value}
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
	
## Loss Taxa

def loss_stats(infile,stdDevs=2,var_filter=3):
	'''From a file of taxon losses, return a generator of csv lines holding:
	gene, mean, variance, and outliers.
	
	stdDevs: Number of stdevs above the mean loss taxa at which to call and algorithm an outlier
	var_filter: Amount of variance at which to consider the possibility of outliers.
	'''
	df = pd.read_csv(infile,index_col=0)
	count = 0
	yield ",".join(['']+["mean","variance","outliers"])
	for ind, row in df.iterrows():
		mean,var = row.mean(), row.var()
		if var <= var_filter:
			outliers = ''
		else:
			outliers = ' '.join((i for i in row[(row-row.mean()) > row.std()*stdDevs].index))
		yield ",".join([ind, str(mean), str(var), outliers])
		count +=1
		if count % 100 == 0:
			sys.stderr.write(str(count)+"\n")