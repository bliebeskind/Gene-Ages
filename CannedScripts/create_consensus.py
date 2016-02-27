from LECA.consensus import consensus_ages
import cPickle as pickle
import sys, os

### This program will create the consensus (mode) age calls 
### by trimming databases that oversplit co-orthologous groups.
###
### **Note: if this script does not find the file LDORESULTS, it will
### silently calculate a consensus without it, so make sure this path is
### correct.


############# User input #######################

INFILE = "binAges_<SPECIES>.csv"
LDORESULTS = "../Errors/Oversplitting/<SPECIES>_LDO_results.p"
FALSEPOSITIVES = "../Errors/Losses/FalsePos_<SPECIES>.p"
TAXON = "<SPECIES>"

## Set one or other to None to create a consensus output without
## filtering algorithms by oversplitting or false positive criteria

#FALSEPOSITIVES = None
#LDORESULTS = None

############ Don't change #######################

with open("../OtherInput/ageLists.p") as f:
	ageLists = pickle.load(f)
assert TAXON in ageLists, "Taxon %s not found in age order file" % TAXON
AGES = ageLists[TAXON]

if os.path.exists(LDORESULTS):
	for line in consensus_ages(INFILE,AGES,LDORESULTS,FALSEPOSITIVES):
		print line
else:
	for line in consensus_ages(INFILE,AGES,LDO_dict=None,lossTaxa_dict=FALSEPOSITIVES):
		print line
