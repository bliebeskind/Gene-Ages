from LECA.node_stats import nodeError,bimodality
from LECA import csv_parser
import cPickle as pickle

### This program will calculate the consistency score for each gene
### It prints a comma separated stream (gene,conScore)

############# User input #######################

INFILE = "nodeAges_<SPECIES>.csv"
CLASS1 = ["InParanoid","InParanoidCore","OMA_Groups","OMA_Pairs","PANTHER8_LDO","RSD"]
CLASS2 = ["EggNOG","Orthoinspector","Hieranoid_2","EnsemblCompara_v2","PANTHER8_all","Metaphors","PhylomeDB"]

############ Don't change #######################

with open("../nodeDists/<SPECIES>_nodeDists.p") as f:
    nodeDistsD = pickle.load(f)

print ",".join(['',"NodeError","Bimodality"])
	
for gene,ageD in csv_parser(INFILE):
	err = str(nodeError(ageD,nodeDistsD))
	bi = str(bimodality(ageD,gene,nodeDistsD,CLASS1,CLASS2))
	print ",".join([gene,err,bi])
