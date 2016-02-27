#! /usr/bin/env python

import pandas as pd
import sys

### This program is for making the final data files. It just concatenates the consensus binned output
### with the node statistics (node error and bimodality) and the HGT flagging information

CONSENSUS="Consensus/consensus_<SPECIES>.csv"
OUTFILE="main_<SPECIES>.csv"
STATS="NodeAges/nodeStats_<SPECIES>.csv"
HGTFLAGS="Errors/Losses/HGTFlag_<SPECIES>.csv"

### Don't change below here ###

sys.stderr.write("Reading files\n")
con = pd.read_csv(CONSENSUS,index_col=0)
stats = pd.read_csv(STATS,index_col=0)
hgt = pd.read_csv(HGTFLAGS,index_col=0)

con.sort(inplace=True)
stats.sort(inplace=True)
hgt.sort(inplace=True)
 
sys.stderr.write("Joining consensus and nodeStats files\n")
tmp = con.join(stats,how='inner')
assert len(tmp) == len(tmp), "nodeStats join was unsuccessful"

sys.stderr.write("Joining consensus and HGT flag file\n")
out = tmp.join(hgt["HGT_flag"],how='inner')
assert len(out) == len(tmp), "HGT flag join was unsuccessful"

sys.stderr.write("Writing output\n")
out.to_csv(OUTFILE)
