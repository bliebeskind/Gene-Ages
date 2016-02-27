from LECA.node_stats import run_LDOcomp, percLDOs
from LECA import csv_parser
import cPickle as pickle
import pandas as pd

### This program runs the oversplitting analysis, used for false negative analysis
### Outputs three files:
###		_results.p		# pickle file used by consensus scripts
###		_results.csv	# the results of the analysis in csv format
###		_summary.csv	# summary statistics (not that useful)

############# User input #######################

ORTHOFILE = "<SPECIES>_coOrthologs.txt"
AGEFILE = "../../NodeAges/nodeAges_<SPECIES>.csv"
YOUNGGROUP = ["InParanoid","InParanoidCore","PANTHER8_LDO","OMA_Groups",
	"OMA_Pairs","RSD","EggNOG","Orthoinspector","Hieranoid_2",
	"EnsemblCompara_v2","Metaphors","PANTHER8_all"]
OLDGROUP = ["PhylomeDB"]
FALSEPOS = "../Losses/lossStats_<SPECIES>.csv"
BINNED=False # Can also calculate LDOs on binned ages, but I wouldn't recommend it. Not conservative enough.


############ Don't change #######################

print "Running LDO analysis"
LDO_results = run_LDOcomp(ORTHOFILE,AGEFILE,OLDGROUP,YOUNGGROUP,FALSEPOS,BINNED)

print "Pickling results to '<SPECIES>_LDO_results.p'"
with open("<SPECIES>_LDO_results.p",'w') as out:
	pickle.dump(LDO_results,out)

print "Writing results to '<SPECIES>_LDO_results.csv'"
LDO_df = pd.DataFrame.from_dict(LDO_results,orient='index')
LDO_df.to_csv("<SPECIES>_LDO_results.csv")

print "Writing summary stats to '<SPECIES>_LDO_summary.csv'"
LDOstats = percLDOs(LDO_results)
LDOstatsDF = pd.DataFrame.from_dict(LDOstats,orient='index')
LDOstatsDF.columns = ["FractionLDOs"]
LDOstatsDF.to_csv("<SPECIES>_LDO_summary.csv")
