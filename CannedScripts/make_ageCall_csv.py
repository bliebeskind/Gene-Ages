from LECA.infer_age import ages_from_tables
from glob import iglob
import sys

#### This script will write csv files of age calls for each database and gene
#### from age tables, which can be found here: http://clairemcwhite.github.io/orthoblender/
####
#### You will also need: 1.) a tree with all nodes labeled in nexus format
####			 2.) a pickled dictionary mapping node labels to clade names (for binned age calls)
####
#### Output will be written to "binAges.csv" if you set BINNED to True, or to "nodeAges.csv" if you
#### set it to False
####
#### By default, the following databases will be in the output:
#### 	"InParanoid","InParanoidCore","OMA_Groups","OMA_Pairs","PANTHER8_LDO","RSD","EggNOG",
#### 	"Orthoinspector","Hieranoid_2","EnsemblCompara_v2","PANTHER8_all","Metaphors","PhylomeDB"

#### User Input ####

BINNED=False # whether you want node ages or binned ages
TREEFILE = "../OtherInput/RefSetSpeciesTree2014_pruned.nex"
MAPPING = "../OtherInput/nodes2taxa_bos.p"   # pickled dictionary mapping node labels to clade names


#### Don't change #####


temp_iter = iglob("*.csv")
file_iter = (i for i in temp_iter if i != "binAges.csv" and i!= "nodeAges.csv")
if BINNED:
    lineGen = ages_from_tables(file_iter,TREEFILE,True,MAPPING)
    with open("binAges.csv",'w') as out:
        for line in lineGen:
	    prot = line.split(",")[0]
	    sys.stderr.write(prot+"\n")
	    out.write(line+"\n")
else:
    lineGen = ages_from_tables(file_iter,TREEFILE)
    with open("nodeAges.csv",'w') as out:
        for line in lineGen:
            prot = line.split(",")[0]
            sys.stderr.write(prot+"\n")
            out.write(line+"\n")

