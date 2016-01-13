from LECA.node_stats import loss_stats
import cPickle as pickle

###### USER INPUT #######

INFILE="lossTaxa.csv"
STATSOUTFILE="lossStats.csv"
PICKLEOUTFILE="FalsePos.p"
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