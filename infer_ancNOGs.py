#! /usr/bin/env python

### For classifying NOGs by their taxonomic distribution. The executable portion 
### below is hard-coded to separate NOGs into four groups: those that are in both 
### Unikonts and Bikonts(LECA-specific NOGs); Unikonts, Bikonts, and Archaea; Unikonts,
### Bikonts, and Eu-Bacteria; and all four. It writes the NOGs to four files
### corresponding to these distributions.

import sys

def opener(infile):
	with open(infile) as f:
		firstline = f.readline() # skip header
		euks = firstline.strip().split()[3:]
		return f
				
def infer_from_counts(src,euk_order=["Bikont","Unikont","Metazoa","Bilateria","Vertebrata","Mammalia"]):
	'''
	For inferring from a file similar to HumanDiseaseNOGs_taxonCount.tsv
	'''
	for line in src:
		line = line.strip().split()
		nog = line[0]
		numArc,numBac = eval(line[1]),eval(line[2])
		eukD = dict(zip(euks,line[3:]))
		for level in euk_order: # move up the tree
			try:
				count = eval(eukD[level])
			except IndexError: # eukD and euks don't match
				raise Exception("%s not found in columns" % level)
			if count == 0:
				continue
			elif level == "Bikont" and count >0: # LECA
				if numBac >0:
					if numArc >0:
						yield "\t".join([nog,"LUCA"])
					else:
						yield "\t".join([nog,"LECA+Bac"])
				elif numArc >0:
					yield "\t".join([nog,"LECA+Arc"])
				else:
					yield "\t".join([nog,"LECA"])
				break
			else:
				assert count >0
				yield "\t".join([nog,level])
				break
		else:# no breaks reached in for loop - NOT WORKING ??
			raise Exception("all levels are zero: %n" % nog)
				
if __name__ == '__main__':
	try:
		infile = sys.argv[1]
		gen = infer_from_counts(opener(infile))
	except IndexError:
		infile = sys.stdin
		gen = infer_from_counts(infile)
	for i in gen:
		print i

