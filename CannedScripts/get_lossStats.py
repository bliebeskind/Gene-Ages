from LECA.node_stats import loss_stats
import cPickle as pickle

### Calculates statistics for false positive analysis.
### Writes two files, one .csv file and one pickle file (binary) used for
### generating a consensus.
###
### The user can tweak the sensitivity of this analysis by changing how many std. devs.
### above the mean number of loss taxa an algorithm must be before being flagged as a
### false positive.

###### USER INPUT #######

INFILE="lossTaxa_<SPECIES>.csv"
STATSOUTFILE="lossStats_<SPECIES>.csv"
PICKLEOUTFILE="FalsePos_<SPECIES>.p"
PICKLEPROTOCOL=pickle.HIGHEST_PROTOCOL
STDEVS=2 # throw out algorithm if lossTaxa are this many stdevs above the mean

##### DON'T Change ######

with open(STATSOUTFILE,'w') as out:
	for line in loss_stats(INFILE,STDEVS):
		out.write(line+"\n")
		
with open(STATSOUTFILE) as stats:
	stats.readline() # skip header
	statD = {}
	for line in stats:
		line = line.strip().split(",")
		statD[line[0]] = line[3].split()
	with open(PICKLEOUTFILE,'w') as out:
		pickle.dump(statD,out,protocol=PICKLEPROTOCOL)
