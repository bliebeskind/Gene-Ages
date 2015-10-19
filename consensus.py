import math
import cPickle as pickle
from LECA import csv_parser
from collections import Counter
	
### Functions to calculate LDO-corrected consensus


def _load_LDO_dict(infile):
	'''Load dictionary mapping genes to a dictionary of Databases:True/False, i.e.
	whether each database was detected as having an LDO for that gene.'''
	with open(infile) as f:
		return pickle.load(f)
		
def _ageDist_gen(infile,LDO_dict=None,filterLDOs=False):
	'''
	Loop over lines in infile and return a generator of (gene, age distribution)
	Filters out databases that have evidence of an LDO break and those that miss that gene
	'''
	if filterLDOs:
		ldos = _load_LDO_dict(LDO_dict)
	for gene,line in csv_parser(infile):
		ageVec = []
		num_ldos = 0
		for db,age in line.iteritems():
			if filterLDOs:
				if gene in ldos and db in ldos[gene] and ldos[gene][db] == True:
					num_ldos +=1
					continue
			if age == 'None':
				continue
			else:
				ageVec.append(age)
		ageCounts = Counter(ageVec) # count age occurrences
		numDBsContributing = len(ageVec)
		assert numDBsContributing > 0, "No databases have age information for %s" % gene
		normCounts = [(i,float(j)/numDBsContributing) for i,j in ageCounts.iteritems()] # normalize
		yield gene, sorted(normCounts, key=lambda x:x[1]), numDBsContributing, num_ldos # sort
	
def consensus_ages(infile,ages,LDO_dict=None,filterLDOs=True):
	'''
	Create csv file holding the distribution over age calls and the score columns "modeAge", 
	"NumDBsContributing", "NumDBsFiltered", "entropy"
	'''
	yield ",".join([''] + [age for age in ages] + ["modeAge"] + ["NumDBsContributing"] + ["NumDBsFiltered"] + ["entropy"]) # header
	func = lambda x: str(ageD[x]) if x in ageD else "0"
	for gene, ageProbs, numDBs, num_ldos in _ageDist_gen(infile,LDO_dict,filterLDOs):
		try:
			modeAge = ageProbs[-1][0] # because sorted in _ageDist_gen
		except IndexError:
			print ageProbs
		ageD = dict(ageProbs)
		entropy = -(sum(x * math.log(x) for x in ageD.itervalues()))
		yield ",".join([gene] + [func(i) for i in ages] + [modeAge] + [str(numDBs)] + [str(num_ldos)] + [str(entropy)])