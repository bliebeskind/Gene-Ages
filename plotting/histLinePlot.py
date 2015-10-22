import gc
import numpy as np
import pandas as pd
from scipy import stats

## Functions for creating a histogram/line plot

def _myRange(start,stop,step):
    num = start
    while num < stop:
        yield num
        num += step

def _getHistBins(histScore,binNumber):
	descr = histScore.describe()
	max,min = descr["max"], descr["min"]
	step = (max - min)/binNumber
	halfStep = step/2.0
	return {i+halfStep:(i,i+step) for i in _myRange(min,max,step)}

def binLineScore(df,lineScoreCol,histScoreCol,binNumber=50):
	'''
	Given a DataFrame with two columns of scores, bin one score into bins based
	on a histogram of the second score. This way, the first score can be plotted
	as a line over the histogram. Builds a dictionary mapping each bin#, which is 
	the midway point for each bin from the histogram, to a list of scores from the 
	line column. Must be fed to get_stats to be plotted.
	
	df has scores as columns: "lineScoreCol" and "histScoreCol", and genes as index
	Set binNumber to the desired number of histogram bins.
	'''
	assert histScoreCol in df.columns, "Column: %s not found" % histScoreCol
	assert lineScoreCol in df.columns, "Column: %s not found" % lineScoreCol
	D = {}
	newDF = df.loc[:,[lineScoreCol,histScoreCol]]
	df = None
	gc.collect() # save some memory
	count = 0
	bins = _getHistBins(newDF[histScoreCol],binNumber)
	for gene,row in newDF.iterrows():
		lineScore,histScore = row[lineScoreCol], row[histScoreCol]
		if np.isnan(lineScore) or np.isnan(histScore): # skip NaNs
			continue
		for bin in bins:
			if bins[bin][0] <= histScore < bins[bin][1]:
				lineScoreBin = bin
				break
		else:
			raise Exception("no bin found for score %f, gene %s" % (lineScore,i))
		if lineScoreBin in D:
			D[lineScoreBin].append(lineScore)
		else:
			D[lineScoreBin] = [lineScore]
		count += 1
		if count %100 == 0:
			print count
	return D
	
def getLineScoreStats(df,lineScoreCol,histScoreCol,binNumber=50):
	'''Return a Dataframe of line score stats for each bin. Relevant
	one is probably the mean.'''
	D = {}
	binnedScores = binLineScore(df,lineScoreCol,histScoreCol,binNumber)
	for bin in binnedScores:
		L = binnedScores[bin]
		if len(L) <=1:
			mean,var,dev = L[0],0,0
			continue
		mean = stats.tmean(L)
		var = stats.tvar(L)
		stanD = stats.tstd(L)
		D[bin] = {"mean":mean,"var":var,"stanDev.": stanD}
	return pd.DataFrame(D).T