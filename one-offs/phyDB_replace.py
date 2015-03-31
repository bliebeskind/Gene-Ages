#! /usr/bin/env python

### Replace PhylomeDB names in orthologs file with Uniprot names

import sys

D = {}
with open("phylome2Ensembl_human.txt") as f:
	f.readline() # skip header
	for line in f:
		line = line.strip().split('\t')
		D[line[0]] = line[2]
		
count = 0
with open("sorted_orthologs.txt") as f:
	print f.readline()
	for line in f:
		line = line.strip().split("\t")
		line[0] = D[line[0]]
		print "\t".join(line)
		count +=1
		if count%100 == 0:
			sys.stderr.write(str(count)+"\n")


#with open("orthologs_sorted.txt") as f:
#	print f.readline()
#	for line in f:
#		line = line.strip().split("\t")
#		newline = []
#		newline.append(D[line[0]])
#		newline.append(D[line[1]])
#		newline += line[2:5]
#		if len(line) > 5:
#			coLogs = line[5].split()
#			newCoLogs = [D[log] for log in coLogs]
#			newline.append(' '.join(newCoLogs))
#		print "\t".join(newline)
#		count +=1
#		if count%100 == 0:
#			sys.stderr.write(str(count)+"\n")
	