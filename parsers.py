#! /usr/bin/env python

def phylome_parser(infile,as_taxid=False,taxonD=None,type_filter=None):
	'''
	Parses a phylomeDB ortholog file
	type_filter is a list of ortholog relation types to consider. Must be one of:
		"many-to-many","many-to-one","one-to-many", or "one-to-one"
		All types not in type_filter are returned by this function
	'''
	values = ("gene","ortholog","type","CS","trees","co-orthologs")
	with open(infile) as f:
		for line in f:
			line = line.strip().split("\t")
			if len(line) > 5:
				line[-1] = line[-1].split(' ')
			if type_filter:
				types = ["many-to-many","many-to-one","one-to-many","one-to-one"]
				type = line[2]
				assert type in types, "type %s not recognized" % type
				if type not in type_filter:
					continue
			lineD = dict(map(None,values,line))
			if as_taxid:
				assert taxonD != None, "Must specify a taxon dictionary"
				lineD["ortholog"] = taxonD[lineD["ortholog"].split("_")[1]]
			yield lineD
		
			
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