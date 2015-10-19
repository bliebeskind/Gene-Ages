from LECA.node_stats import run_LDOcomp, percLDOs
from LECA import csv_parser
import cPickle as pickle
import pandas as pd

############# User input #######################

ORTHOFILE = "coOrthoGroups.txt"
AGEFILE = "nodeAges_HUMAN.csv"
YOUNGGROUP = ["InParanoid","InParanoidCore","OMA_Groups","OMA_Pairs","PANTHER8_LDO","RSD"]
OLDGROUP = ["Metaphors"]
BINNED=False


############ Don't change #######################

print "Running LDO analysis"
LDO_results = run_LDOcomp(ORTHOFILE,AGEFILE,OLDGROUP,YOUNGGROUP,BINNED)

print "Pickling results to 'LDO_results.p'"
with open("LDO_results.p",'w') as out:
	pickle.dump(LDO_results,out)

print "Writing results to 'LDO_results.csv'"
LDO_df = pd.DataFrame.from_dict(LDO_results,orient='index')
LDO_df.to_csv("LDO_results.csv")

print "Writing summary stats to 'LDO_summary.csv'"
LDOstats = percLDOs(LDO_results)
LDOstatsDF = pd.DataFrame.from_dict(LDOstats,orient='index')
LDOstatsDF.columns = ["FractionLDOs"]
LDOstatsDF.to_csv("LDO_summary.csv")