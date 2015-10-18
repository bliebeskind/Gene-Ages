from LECA.node_stats import polarization
from LECA import csv_parser
import cPickle as pickle
import sys

INFILE = "nodeAges_HUMAN.csv"
CLASS1 = ["InParanoid","InParanoidCore","OMA_Groups","OMA_Pairs","PANTHER8_LDO","RSD"]
CLASS2 = ["Orthoinspector","Hieranoid_2","EnsemblCompara_v2","PANTHER8_all","Metaphors","PhylomeDB"]


with open("NodeDists.p") as f:
    nodeDistsD = pickle.load(f)

for gene,ageD in csv_parser(INFILE):
    print polarization(ageD,gene,nodeDistsD,CLASS1,CLASS2)
