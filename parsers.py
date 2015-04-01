#! /usr/bin/env python

def phylome_parser(infile):
	'''For parsing a phylomeDB ortholog file'''
	values = ("prot","ortholog","type","CS","trees","co-orthologs")
	with open(infile) as f:
		for line in f:
			line = line.strip().split("\t")
			line[-1] = line[-1].split(' ')
			assert len(line) == len(values)
			yield dict(zip(values,line))