import pandas as pd
from LECA.db_comp import avg_dist

INFILE = "nodeAges_HUMAN.csv"
NODEDISTS = "NodeDists.p"

#####

avgDistsDF = avg_dist(INFILE,NODEDISTS)
with open("patristic_distances.nex",'w') as out:
	# Write nexus headers
	out.write("#NEXUS\n\n\n")
	out.write("Begin distances;\nDimensions ntax=%i;\nformat nodiagonal;\nmatrix\n" % len(avgDistsDF.columns))
	
	# Write matrix
	out.write(avgDistsDF.columns[0] + "\n")
	for i,j in avgDistsDF.iterrows():
		out.write(" ".join([i]+[str(x) for x in j if str(x) != 'nan']) + "\n")
	out.write(";\nend;\n\n")
	
	# Write PAUP commands
	out.write("begin paup;\nlog file=paup.log;\ndset distance=user;\nhsearch;\nsavetrees file=patristicDistances.tre brlens=yes;\nend;")


