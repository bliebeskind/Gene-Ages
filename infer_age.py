#! /usr/bin/env python

import dendropy
import sys, pickle

alternate_names = {'ASPFC': 'ASPFU',
	'BACT4': 'BACTN',
	'CAEBR': 'CAEEL',
	'CANAW': 'CANAL',
	'CANAX': 'CANAL',
	'CHLTA': 'CHLTR',
	'CRYNE': 'CRYNJ',
	'GEOSN': 'GEOSL',
	'HALS3': 'HALSA',
	'LEPIR': 'LEPIN',
	'MYCTO': 'MYCTX',
	'MYCTU': 'MYCTX',
	'NEUCS': 'NEUCR',
	'PHAND': 'PHANO',
	'PSEA7': 'PSEAE',
	'STRCH': 'STRCO',
	'SULSF': 'SULSO',
	'THEMT': 'THEMA',
	'USTMD': 'USTMA',
	'YARLL': 'YARLI',
	'YEASX': 'YEAST'}
	
Archaea = ["THECO","METAC","METJA","HALSA","SULSO"]

# From taxon count file

def counts_gen(src,tree):
	'''
	Read from taxon_count file and yield two-element tuples:
	(NOG, {Metazoa:10, Bilateria:90...}
	'''
	tree_taxa = [str(i) for i in tree.taxon_set]
	infile_taxa = src.next().strip().split()[1:] # taxa should be columns 2-end
	assert sorted(tree_taxa) == sorted(infile_taxa), \
	"Taxa from infile and tree don't match:\n %s\n%s" % (str(tree_taxa),str(infile_taxa))
	for line in src:
		line = line.strip().split()
		count_cols = [int(i) for i in line[1:]]
		assert len(count_cols) == len(infile_taxa)
		yield line[0], {i:j for i,j in zip(infile_taxa,count_cols) if j > 0} # skip zero counts
		
def age_generator(src,tree_source,format='nexus',source_type='file'):
	'''Return generator of tuples where i,j are group,age'''
	tree = get_dendropy_tree(tree_source,format,source_type)
	for group, count in counts_gen(src,tree):
		if len(count.keys()) > 1:
			count_node = tree.mrca(taxon_labels=count.keys())
			age = count_node.annotations.get_value("age")
			if age == "LUCA": # This block reduces the generality, might think about breaking it out
				if "Bacteria" not in count:# into another function
					assert "Archaea" in count, "Neither Archea nor Bacteria found in LUCA group: %s" % group
					yield group, "Euk+A"
				elif "Archaea" not in count:
					assert "Bacteria" in count
					yield group, "Euk+B"
				else:
					yield group, age
			else:
				yield group, age
		else: # just one branch
			yield group, count.keys()[0] # should add error checking if all zeros
		
# from DataBase comparison files (Claire's format)

def get_dendropy_tree(tree_source,format='nexus',source_type='file'):
	'''Read in a tree using dendropy. Tree should have node labels for ancestral nodes that will be
	used for aging. They should be of format: [&age: "Metazoa"]'''
	format,source = format.lower(),source_type.lower()
	tree_funcD = {'string':dendropy.Tree.get_from_string,
				'file':dendropy.Tree.get_from_path,
				'stream':dendropy.Tree.get_from_stream}
	tree = tree_funcD[source_type](tree_source,format,extract_comment_metadata=True)
	return tree
	
def read_dbComp(infile):
	'''Read in one of Claire's ortholog X database files, and return the name of the human protein and
	a dictionary mapping each database to a list of the species with identified orthologs.'''
	protein = infile.split(".")[0]
	with open(infile) as f:
		header = f.readline().strip().split(",")
		dbList = header[3:]
		dbDict = {dbname:[] for dbname in dbList}
		for line in f:
			line = line.strip().split(",")
			try:
				species = line[1].split("_")[1]
			except IndexError:
				raise Exception("bad format: %s, species %s" % (infile, line[1]))
			if species in alternate_names:
				species = alternate_names[species]
			for db,value in zip(dbList,line[3:]):
				if value == "1":
					dbDict[db].append(species)
	return protein, dbDict

		
def get_db_age_nodes(infile,tree,as_clades=False,conversion_dictionary=None):
	'''Open one of Claire's DBcomp files and infer the age of the protein for each database therein.
	Uses a species tree with annotated interior nodes and returns a dictionary mapping each database
	name to the ancestral node representing the inferred age.
	
	Calls read_dbComp'''
	prot, dbD = read_dbComp(infile)
	ageD = {}
	for db in dbD:
		species_set = set(dbD[db])
		if len(species_set) == 0: # shouldn't be, maybe some error handling is in order here.
			ageD[db] = None
			continue
		try:
			ageNode = tree.mrca(taxon_labels=species_set).label
			if ageNode == None:
				ageD[db] = None
				continue
		except KeyError:
			taxon_labels = [i.label for i in tree.taxon_set]
			es = [i for i in species_set if i not in taxon_labels]
			print "Couldn't find taxon %s in protein %s" % (str(es),prot)
			return prot, None
		if as_clades: # use clade names instead of numerical internal node labels
			ageNode = conversion_dictionary[ageNode]
			if ageNode == "Cellular_organisms":
				for i in species_set:
					if i in Archaea:
						break
				else:
					ageNode = "Euk+Bac"
		ageD[db] = ageNode
	return prot, ageD
	
def serialize_dbAgeNodes(infile_stream,tree_source,as_clades=False,conversion_dictionary=None):
	'''Pickle the output of get_db_age_nodes to a file <prot>.p'''
	tree = get_dendropy_tree(tree_source)
	convD = None
	if as_clades:
			assert conversion_dictionary != None, \
			"Must supply a dictionary converting age nodes to taxa"
			with open(conversion_dictionary) as f:
				convD = pickle.load(f)
	for f in infile_stream:
		prot, ageD = get_db_age_nodes(f,tree,as_clades,convD)
		if ageD == None:
			continue
		with open(prot+".p",'w') as out:
			pickle.dump(ageD,out)
		print prot

def ages_from_tables(infile_stream,tree_source,as_clades=False,conversion_dictionary=None,dbs=None):
	'''
	Create a csv file of ages for each gene from age call tables.
	Takes a stream of infiles (e.g. from iglob) and a species tree as input
	By default outputs node labels as ages, but if you want to bin these nodes into clades,
	must supply a conversion dictionary mapping each node to its clade and set as_clades=True
	'''
	if dbs == None:
		columns = ["InParanoid","InParanoidCore","OMA_Groups","OMA_Pairs","PANTHER8_LDO","RSD",
				"Orthoinspector","Hieranoid_2","EnsemblCompara_v2","PANTHER8_all","Metaphors","PhylomeDB"]
	tree = get_dendropy_tree(tree_source)
	convD = None
	if as_clades:
			assert conversion_dictionary != None, \
			"Must supply a dictionary converting age nodes to taxa"
			with open(conversion_dictionary) as f:
				convD = pickle.load(f)
	yield ",".join(['']+columns)
	for f in infile_stream:
		prot, ageD = get_db_age_nodes(f,tree,as_clades,convD)
		if ageD == None:
			continue
		yield ",".join([prot]+[str(ageD[i]) for i in columns])

