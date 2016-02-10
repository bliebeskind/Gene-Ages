from LECA.infer_age import count_lossTaxa
from LECA.functions import nonRedundant_filestream as stream
from glob import iglob
import sys

### Write number of loss taxa from Claire's files ###

INFILES=stream() 
DBS=["InParanoid","InParanoidCore","OMA_Groups","OMA_Pairs","PANTHER8_LDO",\
	"RSD","EggNOG","Orthoinspector","Hieranoid_2","EnsemblCompara_v2",\
	"PANTHER8_all","Metaphors","PhylomeDB"]
TREE="../OtherInput/RefSetSpeciesTree2014_pruned.nex"

count = 0
with open("lossTaxa_<SPECIES>.csv",'w') as out:
	for line in count_lossTaxa(INFILES,TREE,DBS):
		out.write(line + "\n")
		count +=1
		if count % 100 == 0:
			sys.stderr.write(str(count) + "\n")
