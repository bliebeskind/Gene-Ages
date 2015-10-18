from LECA.node_stats import ageConsistency
from LECA import csv_parser
import cPickle as pickle

### This program will calculate the consistency score for each gene
### It prints a comma separated stream (gene,conScore)

############# User input #######################

INFILE = "nodeAges_HUMAN.csv"

############ Don't change #######################

with open("NodeDists.p") as f:
    nodeDistsD = pickle.load(f)

print ",".join(['',"Consistency"])
	
for gene,ageD in csv_parser(INFILE):
    print ",".join([gene,str(ageConsistency(ageD,nodeDistsD))])