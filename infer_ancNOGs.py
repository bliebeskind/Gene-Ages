### For classifying NOGs by their taxonomic distribution. The executable portion 
### below is hard-coded to separate NOGs into four groups: those that are in both 
### Unikonts and Bikonts(LECA-specific NOGs); Unikonts, Bikonts, and Archaea; Unikonts,
### Bikonts, and Eu-Bacteria; and all four. It writes the NOGs to four files
### corresponding to these distributions.

import sys

def make_groups(infile_list):
	'''Assumes a tab delimited file with taxon ids in the second column'''
	num_groups = len(infile_list)
	groupD = {}
	for infile in infile_list:
		with open(infile) as f:
			group = infile.split(".")[0]
			taxIDs = [line.split("\t")[1].strip() for line in f]
			groupD[group] = taxIDs
	return groupD, num_groups
	
def infer_gen(member_file,infile_list):
	'''For iterating over a eggNOG members file. Returns a generator yielding
	NOGs and a list of groups that are in them'''
	groupD,num_groups = make_groups(infile_list)
	with open(member_file) as f:
		kog = None
		groups = []
		line = f.readline() # skip header
		for line in f:
			line = line.strip().split("\t")
			curr_kog = line[0].strip()
			if curr_kog == kog:
				if len(groups) == num_groups: # all groups found
					continue # skip to next iteration
			else: # new kog
				yield kog, groups
				kog = curr_kog
				groups = []
			taxID = line[1].split(".")[0]
			for i in groupD:
				if taxID in groupD[i]:
					group = i
					break
			else:
				sys.stderr.write("%s not found in any group\n" % taxID)
			if group in groups:
				continue
			else:
				groups.append(group)
		
if __name__ == '__main__':
	groups_list = ["Unikonts.txt","Bikonts.txt","Bacteria.txt","Archaea.txt"]
	outfileD = {"LECA":open("LecaNOGs.txt",'a'),
		"EukArch":open("Euk+ArchNOGs.txt",'a'),
		"EukBac":open("Euk+BacNOGs.txt",'a'),
		"LCA":open("LcaNOGs.txt",'a')} # ignores Unikonts+Bacteria etc...
	for i in infer_gen("NOG.members.txt",groups_list):
		groups = i[1]
		if len(groups) <2:
			continue
		elif len(groups) == 2:
			if sorted(groups) == ['Bikonts','Unikonts']:
				outfileD["LECA"].write(i[0]+'\n')
		elif len(groups) == 3:
			if sorted(groups) == ['Archaea','Bikonts','Unikonts']:
				outfileD["EukArch"].write(i[0]+'\n')
			elif sorted(groups) == ['Bacteria','Bikonts','Unikonts']:
				outfileD["EukBac"].write(i[0]+'\n')
			else:
				continue # could be, e.g. ['Bacteria', 'Unikonts', 'Archaea']
		else:
			if sorted(groups) == ['Archaea','Bacteria','Bikonts','Unikonts']:
				outfileD["LCA"].write(i[0]+'\n')
			else:
				raise Exception("unexpected list %s" % str(groups))