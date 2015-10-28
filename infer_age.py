#! /usr/bin/env python

import dendropy
import sys, pickle
from LECA.functions import flatten
from LECA.parsers import phylome_parser

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

# Shared Functions

def get_dendropy_tree(tree_source,format='nexus',source_type='file'):
	'''Read in a tree using dendropy. Tree should have node labels for ancestral nodes that will be
	used for aging. They should be of format: [&age: "Metazoa"]'''
	format,source = format.lower(),source_type.lower()
	tree_funcD = {'string':dendropy.Tree.get_from_string,
				'file':dendropy.Tree.get_from_path,
				'stream':dendropy.Tree.get_from_stream}
	tree = tree_funcD[source_type](tree_source,format,extract_comment_metadata=True)
	return tree
	
# From a single database

def _stream_from_phylome(infile,type_filter=["many-to-many","many-to-one","one-to-many","one-to-one"]):
	'''Return stream of (protein,(species1,species2)) from a phylomeDB orthologs file.
	
	species output: if call is Phy003II51_MACMU , MACMU will be returned.
	type_filter: list of what kinds of orthology relationship to consider.'''
	return flatten(((i["gene"],i["ortholog"].split("_")[1]) for i in phylome_parser(infile,type_filter=type_filter)))
		
def age_generator(src,tree_source,conversion_dictionary=None,as_clades=False,tree_format='nexus'):
	'''Takes in a stream of (prot,[taxon1,taxon2]) tuples, and returns a stream 
	of (protein,age) tuples calculated by phylostratigraphy on a supplied tree.'''
	tree = get_dendropy_tree(tree_source,tree_format)
	strTaxonSet = [i.label for i in tree.taxon_namespace] # get strings of taxa in tree
	count = 0
	for prot, taxa in src:
		callSet = set([i for i in taxa if i in strTaxonSet]+["HUMAN"]) # ignore taxa not in tree
		if len(callSet) <=1:
			print "nothing found for %s" % prot
			continue
		ageNode = tree.mrca(taxon_labels=callSet).label # can be None
		if as_clades: # use clade names instead of numerical internal node labels
			ageNode = conversion_dictionary[ageNode]
			if ageNode == "Cellular_organisms":
				for i in callSet:
					if i in Archaea:
						break
				else:
					ageNode = "Euk+Bac"
		yield prot, ageNode
		count +=1
		if count % 100 == 0:
			sys.stderr.write(str(count)+"\n")
		
# from DataBase comparison files (Claire's format)

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
		dbs = ["InParanoid","InParanoidCore","OMA_Groups","OMA_Pairs","PANTHER8_LDO","RSD",
				"Orthoinspector","Hieranoid_2","EnsemblCompara_v2","PANTHER8_all","Metaphors","PhylomeDB"]
	tree = get_dendropy_tree(tree_source)
	convD = None
	if as_clades:
			assert conversion_dictionary != None, \
			"Must supply a dictionary converting age nodes to taxa"
			with open(conversion_dictionary) as f:
				convD = pickle.load(f)
	yield ",".join(['']+dbs)
	for f in infile_stream:
		prot, ageD = get_db_age_nodes(f,tree,as_clades,convD)
		if ageD == None:
			continue
		yield ",".join([prot]+[str(ageD[i]) for i in dbs])

def count_losses(infile_stream,tree_source,dbs=None):
	'''
	Create csv file of number of inferred gene losses for each database. The losses are calculated by 
	subtracting the number of species with orthologs from the total number of descendants of the ancestral
	node.
	'''
	if dbs == None:
		dbs = ["InParanoid","InParanoidCore","OMA_Groups","OMA_Pairs","PANTHER8_LDO","RSD",
				"Orthoinspector","Hieranoid_2","EnsemblCompara_v2","PANTHER8_all","Metaphors","PhylomeDB"]
	tree = get_dendropy_tree(tree_source)
	yield ",".join(['']+dbs)
	for f in infile_stream:
		lossD = {}
		prot, parsed = read_dbComp(f)
		prot, ageD = get_db_age_nodes(f,tree)
		if ageD == None:
			continue
		for db in dbs:
			assert db in ageD, "Database %s not found in %s" % (db,f)
			node = ageD[db]
			nodeObj = tree.find_node_with_label(node)
			assert nodeObj != None, "Internal label %s not found" % node
			numDescendants = len(set((i.taxon.label for i in nodeObj.leaf_nodes())))
			numOrthologs = len(set(parsed[db]))
			lossD[db] = numDescendants - numOrthologs
		yield ",".join([prot]+[str(lossD[i]) for i in dbs])