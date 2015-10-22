import sys
from collections import Counter

def parse(infile):
	'''Parse a CORUM database file yielding dictionaries for each line'''
	with open(infile) as f:
		header = f.readline().strip().split("\t")
		length = len(header)
		for line in f:
			line = line.strip().split("\t")
			if len(line) != length:
				print "Line %s is not complete" % line[0]
				continue
			outD = dict(zip(header,line))
			try:
				outD["subunits (UniProt IDs)"] = \
				outD["subunits (UniProt IDs)"].strip().replace("(","").replace(")","").split(",")
			except KeyError:
				raise Exception("Must have column 'subunits (UniProt IDs)'")
			yield outD
			
def renameProts(infile,convD):
	'''Rename proteins using conversion dictionary'''
	with open(infile) as f:
		header = f.readline().strip()
		yield header
		numChanged = 0
		for line in f:
			line = line.strip().split("\t")
			if len(line) != len(header.strip().split("\t")):
				print "Line %s is not complete" % line[0]
				continue
			try:
				prots = line[2].strip().strip().replace("(","").replace(")","").split(",")
			except IndexError, e:
				raise Exception("%s: line %s" % (e, line[0]))
			for p in prots:
				if p in convD:
					prots.remove(p)
					prots.append(p)
					numChanged += 1
			outline = "\t".join(line[:2] + [",".join(prots)] + line[5:])
			yield outline
		print numChanged
		
def mapProt2Complex(infile):
	'''Parse a CORUM data file and return a dictionary
	mapping each protein to the complex(es) it's in.'''
	outD = {}
	for line in parse(infile):
		prots = line["subunits (UniProt IDs)"]
		complexName = line["Complex name"]
		for p in prots:
			if p in outD:
				outD[p] += [complexName]
			else:
				outD[p] = [complexName]
	return outD

def age_dist(infile,ageD):
	'''Get distribution of ages from whole file like Human_3cols'''
	## lots of overlapping code with complex_age
	allCounts = Counter()
	with open(infile) as f:
		for line in f:
			line = line.strip().split("\t")
			if len(line) < 3: # skip ones without proteins
				continue
			subunits = line[2].replace("(","").replace(")","").split(",")
			age_list = []
			for i in subunits:
				try:
					age_list.append(ageD[i])
				except KeyError:
					sys.stderr.write("Couldn't find %s\n" % i)
					continue # proteins that aren't found in D are ignored as missing data
			ages = Counter(age_list)
			allCounts = allCounts + ages
	return allCounts

def funcat_count(infile,funcatD):
	funcat_list = sorted(funcatD.keys())
	yield "\t".join([''] + funcat_list)
	with open(infile) as f:
		for line in f:
			line = line.strip().split("\t")
			assert len(line) == 2, "%s" % line
			complex,cats = line[0],line[1]
			cats_list = cats.split(",")
			# should have something checking items in cats_list in funcat_set
			cats_presence = map(lambda x: "1" if x in cats_list else "0", funcat_list)
			yield "\t".join([complex] + cats_presence)
	
def name_funcats(infile,funcatD):
	with open(infile) as f:
		for line in f:
			line = line.strip().split("\t")
			complex,cats = line[0],line[1].split(",")
			try:
				named_cats = ",".join([funcatD[cat] for cat in cats])
			except KeyError, e:
				raise Exception("Couldn't find category %s in %s" % (e, line[1]))
			yield "\t".join([complex] + [named_cats])

#Should use FunCat dictionary instead of this set	
def get_funcat_set(infile):
	catSet = set([])
	with open(infile) as f:
		for line in f:
			line = line.strip().split("\t")
			try:
				cats = set(line[1].split(","))
			except IndexError:
				print "Couldn't find annotations for %s" % line[0]
			catSet = catSet.union(cats)
	return catSet
			

def funcat_level(src,level=1,sep="\t",header=True):
	'''
	src should be generator of (complex,funcat annotation).
	Highest level is "1".
	FunCat annotations should be joined within the column by comma, so done't use comma for "sep"
	'''
	for i,j in src:
		inlevels_list = j.split(",")
		outlevels_list = [".".join([] + l.split(".")[:level]) for l in inlevels_list]
		outlevels = ",".join(outlevels_list)
		yield sep.join((i,outlevels))

def filter_zeros(infile,outfile,header=True,threshold=1):
	'''Remove lines with zero counts from an age counts file.'''
	skipped = []
	inf,outf = open(infile),open(outfile,'w')
	if header:
		outf.write(inf.readline())
	for line in inf:
		line = line.strip().split("\t")
		rowSum = sum([int(i) for i in line[1:]])
		if rowSum >= threshold:
			outf.write("\t".join(line) + "\n")
		else:
			skipped.append(line[0])
	return skipped

def complex_age(infile,D):
	'''From a file like Human_3cols return a dictionary of counts for each category in each complex
	D: a dictionary mapping proteins to ages.'''
	with open(infile) as f:
		f.next()
		for line in f:
			line = line.strip().split("\t")
			if len(line) < 3: # skip ones without proteins
				continue
			subunits = line[2].replace("(","").replace(")","").split(",")
			age_list = []
			for i in subunits:
				try:
					age_list.append(D[i])
				except KeyError:
					sys.stderr.write("Couldn't find %s\n" % i)
					continue # proteins that aren't found in D are ignored as missing data
			ages = Counter(age_list)
			yield line[0],line[1],ages

def print_complex_age(infile,D,age_list):
	'''Print a tab delimited file of counts of each age category output by complex age.'''
	age_list = sorted(age_list)
	gen = complex_age(infile,D)
	yield "\t".join(["complex"] + age_list)
	for index,complex,ages in gen:
		complex_ages = [str(ages[i]) for i in age_list]
		assert len(complex_ages) == len(age_list)
		yield "\t".join([complex]+complex_ages)