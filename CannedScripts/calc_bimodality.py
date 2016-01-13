from LECA.node_stats import bimodality
from LECA import csv_parser
import cPickle as pickle

### This program will calculate the bimodality statistic or metric.
### It prints a comma separated stream (gene,pol)

############# User input #######################

INFILE = "nodeAges_HUMAN.csv"
CLASS1 = ["InParanoid","InParanoidCore","OMA_Groups","OMA_Pairs","PANTHER8_LDO","RSD"]
CLASS2 = ["Orthoinspector","Hieranoid_2","EnsemblCompara_v2","PANTHER8_all","Metaphors","PhylomeDB"]

############ Don't change #######################

with open("NodeDists.p") as f:
    nodeDistsD = pickle.load(f)
	
print ",".join(['',"Bimodality"])

for gene,ageD in csv_parser(INFILE):
    print ",".join([gene,str(bimodality(ageD,gene,nodeDistsD,CLASS1,CLASS2))])
