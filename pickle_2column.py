#! /usr/bin/env python

### Takes a two column file and pickles a dictionary mapping the first to
### the second column.
###
### The first argument is the infile, the second is the intended name of the
### pickle file.


import sys,pickle

infile,pickle_file = sys.argv[1],sys.argv[2]

D = {}
count = 0
with open(infile) as f:
    for line in f:
	line = line.strip().split()
	D[line[0]] = line[1]
	count +=1

with open(pickle_file,'w') as f:
    pickle.dump(D,f)

print "Pickled %i lines" % count
