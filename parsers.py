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
			
def eggnog_parser(infile,header_string=None):
	'''Parses an eggnog .members.txt file. Returns generator of 
	dicationaries for each line.'''
	values = ("nog","taxonID","prot_name","start","end")
	with open(infile) as f:
		if header_string != None: 	# if header string
			first = f.readline()	# skip header
			assert first.startswith(header_string),"No header starting with: %s" % header_string
		for line in f:
			line = line.strip().split("\t")
			taxon,prot = line[1][:line[1].find(".")],line[1][line[1].find(".")+1:]
			line_list = [line[0],taxon,prot] + line[2:]
			assert len(line_list) == len(values)
			yield dict(zip(values,line_list))