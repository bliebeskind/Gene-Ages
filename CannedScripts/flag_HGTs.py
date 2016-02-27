#! /usr/bin/env python

import pandas as pd

## This program flags genes as possible HGT events.
## How exactly this is done is covered in a ipython notebook.

############# User Input ##############

LOSSES="lossTaxa_HUMAN.csv"
STATS="lossStats_HUMAN.csv"
OUTFILE="HGTFlag_HUMAN.csv"

####### Don't change below here #######


def make_outlierD(infile):
	outlierD = {}
	with open(infile) as f:
		f.readline() # skip header
		for line in f:
			line = line.strip().split(",")
			if line[3] == '':
				continue
			outlierD[line[0]] = line[3].split()
	return outlierD
	
def create_averages(df,outlierD):
	count = 0
	avgs = pd.Series()
	for index,row in df.iterrows():
		assert index not in avgs, "Label found twice: %s" % index
		count +=1
		if count % 100 == 0:
			print count
		if index in outlierD: # gene had algorithm outliers
			dbs = [i for i in df.columns if i not in outlierD[index]]
		else:
			dbs = df.columns
		avgs[index] = row[dbs].mean() # only include algorithm that have not been trimmed
	return avgs
	
def _floatRange(start,stop,step):
    i = start
    while i <= stop:
        yield i
        i += step
	
if __name__ == '__main__':
	lossTaxa = pd.read_csv(LOSSES,index_col=0)
	
	print "Reading outliers"
	outlierD = make_outlierD(STATS)
	
	print "Calculating averages"
	averages = create_averages(lossTaxa,outlierD)
	lossTaxa["Avg"] = averages
	
	print "Flagging"
	quantile_steps = [i for i in _floatRange(0,1,.05)]   
	quantiles = lossTaxa["Avg"].quantile(quantile_steps)
	cutoff = quantiles.iloc[-1]
	lossTaxa["HGT_flag"] = lossTaxa["Avg"] >= cutoff
	
	print "Writing to file"
	lossTaxa.to_csv(OUTFILE)
