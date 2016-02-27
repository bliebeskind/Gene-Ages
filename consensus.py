import sys
import math
import pandas as pd
import cPickle as pickle
from LECA import csv_parser
from collections import Counter
	
### Functions to calculate LDO-corrected consensus


def _load_pickle(infile):
	'''Load dictionary mapping genes to a dictionary of Databases:True/False, i.e.
	whether each database was detected as having an LDO for that gene.'''
	with open(infile) as f:
		return pickle.load(f)
		
def _ageDist_gen(infile,LDO_dict=None,lossTaxa_dict=None):
	'''
	Loop over lines in infile and return a generator of (gene, age distribution)
	Filters out databases that have evidence of an LDO break and those that miss that gene
	'''
	if LDO_dict:
		ldos = _load_pickle(LDO_dict)
	if lossTaxa_dict:
		lossTaxa = _load_pickle(lossTaxa_dict)
	for gene,line in csv_parser(infile):
		ageVec = []
		num_ldos = 0
		num_lossTaxa = 0
		for db,age in line.iteritems():
			if LDO_dict:
				if gene in ldos and db in ldos[gene] and ldos[gene][db] == True:
					num_ldos += 1
					continue
			if lossTaxa_dict:
				if gene in lossTaxa and db in lossTaxa[gene]:
					num_lossTaxa += 1
					continue
			if age == 'None':
				continue
			else:
				ageVec.append(age)
		ageCounts = Counter(ageVec) # count age occurrences
		numDBsContributing = len(ageVec)
		try:
			assert numDBsContributing > 0
		except AssertionError:
			sys.stderr.write("No databases have age information for %s\n" % gene)
			continue
		normCounts = [(i,float(j)/numDBsContributing) for i,j in ageCounts.iteritems()] # normalize
		yield gene, sorted(normCounts, key=lambda x:x[1]), numDBsContributing, num_ldos+num_lossTaxa # sort
		
def _get_ages(infile):
	df = pd.read_csv(infile,index_col=0,na_values=["None"])
	ages = set()
	for i in df:
		ages = ages.union(set(df[i].value_counts().index))
	return ages
	
def consensus_ages(infile,ages,LDO_dict=None,lossTaxa_dict=None):
	'''
	Create csv file holding the distribution over age calls and the score columns "modeAge", 
	"NumDBsContributing", "NumDBsFiltered", "entropy"
	'''
	yield ",".join([''] + [age for age in ages] + ["modeAge"] + ["NumDBsContributing"] + ["NumDBsFiltered"] + ["entropy"]) # header
	func = lambda x: str(ageD[x]) if x in ageD else "0"
	for gene, ageProbs, numDBs, num_filtered in _ageDist_gen(infile,LDO_dict,lossTaxa_dict):
		try:
			modeAge = ageProbs[-1][0] # because sorted in _ageDist_gen
		except IndexError:
			print ageProbs
		ageD = dict(ageProbs)
		entropy = -(sum(x * math.log(x) for x in ageD.itervalues()))
		yield ",".join([gene] + [func(i) for i in ages] + [modeAge] + [str(numDBs)] + [str(num_filtered)] + [str(entropy)])