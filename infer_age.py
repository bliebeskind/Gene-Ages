#! /usr/bin/env python

import dendropy
import sys

def get_dendropy_tree(tree_source,format='nexus',source_type='file'):
	'''Read in a tree using dendropy. Tree should have node labels for ancestral nodes that will be
	used for aging. They should be of format: [&age: "Metazoa"]'''
	format,source = format.lower(),source_type.lower()
	tree_funcD = {'string':dendropy.Tree.get_from_string,
				'file':dendropy.Tree.get_from_path,
				'stream':dendropy.Tree.get_from_stream}
	tree = tree_funcD[source_type](tree_source,format,extract_comment_metadata=True)
	return tree
				
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
		
				
if __name__ == '__main__':
	try:
		infile,tree = sys.argv[1],sys.argv[2]
		with open(infile) as f:
			for i,j in age_generator(f,tree,'nexus','file'):
				print "\t".join(i,j)
	except IndexError:
		src = sys.stdin
		for i,j in age_generator(src,tree,'nexus','file'):
			print "\t".join(i,j)

