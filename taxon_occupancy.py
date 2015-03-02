#! /usr/bin/env python

import sys
from infer_ancNOGs import make_groups
	
def occupancy(member_file,infile_list):
	groupD,num_groups = make_groups(infile_list)
	with open(member_file) as f:
		kog = None
		groups = []
		group_counts = {g:0 for g in groupD}
		line = f.readline() # skip header
		for line in f:
			line = line.strip().split("\t")
			curr_kog = line[0].strip()
			if curr_kog == kog:
				pass
				#if len(groups) == num_groups: # all groups found
					#continue # skip to next iteration
			else: # new kog
				yield kog, groups, group_counts
				kog = curr_kog
				groups = []
				group_counts = {g:0 for g in groupD}
			taxID = line[1].split(".")[0]
			for i in groupD:
				if taxID in groupD[i]:
					group = i
					assert group in group_counts, \
					"group %s not found in group_counts" % group
					group_counts[group] +=1
					break
			else:
				sys.stderr.write("%s not found in any group\n" % taxID)
				continue
			if group in groups:
				continue
			else:
				groups.append(group)
				
if __name__ == '__main__':
	groups_list = ["Unikonts.txt","Bikonts.txt","Bacteria.txt","Archaea.txt"]
	outfileD = {"LECA":open("Leca_taxonDist.txt",'a'),
		"EukArch":open("Euk+Arch_taxonDist.txt",'a'),
		"EukBac":open("Euk+Bac_taxonDist.txt",'a'),
		"LCA":open("Lca_taxonDist.txt",'a')} # ignores Unikonts+Bacteria etc...
	for i in occupancy("NOG.members.txt",groups_list):
		groups = i[1]
		if len(groups) <2:
			continue
		elif len(groups) == 2:
			if sorted(groups) == ['Bikonts','Unikonts']:
				taxonCounts = i[2]
				try:
					logos = str(taxonCounts['Unikonts'] - taxonCounts['Bikonts'])
				except ZeroDivisionError:
					print str(taxonCounts)
					break
				outfileD["LECA"].write('\t'.join([i[0],logos])+'\n')
		elif len(groups) == 3:
			if sorted(groups) == ['Archaea','Bikonts','Unikonts']:
				taxonCounts = i[2]
				try:
					logos = \
					str((taxonCounts['Unikonts'] + \
					taxonCounts['Bikonts']) - taxonCounts['Archaea'])
				except ZeroDivisionError:
					print str(taxonCounts)
					break
				outfileD["EukArch"].write('\t'.join([i[0],logos])+'\n')
			elif sorted(groups) == ['Bacteria','Bikonts','Unikonts']:
				taxonCounts = i[2]
				try:
					logos = \
					str((taxonCounts['Unikonts'] + \
					taxonCounts['Bikonts']) - taxonCounts['Bacteria'])
				except ZeroDivisionError:
					print str(taxonCounts)
					break
				outfileD["EukBac"].write('\t'.join([i[0],logos])+'\n')
			else:
				continue # could be, e.g. ['Bacteria', 'Unikonts', 'Archaea']
		else:
			continue
#			if sorted(groups) == ['Archaea','Bacteria','Bikonts','Unikonts']:
#				outfileD["LCA"].write(i[0]+'\n')
#			else:
#				raise Exception("unexpected list %s" % str(groups))