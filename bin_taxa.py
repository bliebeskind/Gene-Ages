### Simple script to parse taxonomy files.

import sys

def get_taxonomic_level(infile,euk_level=False):
	euk_level = eval(euk_level)
	with open(infile) as f:
		for line in f:
			line = line.strip().split("\t")
			if len(line) < 2:
				continue
			else:
				lineageD = eval(line[2])
				kingdom = lineageD["Lineage"].split(';')[1].strip()
				if euk_level:
					if kingdom == "Eukaryota":
						level = lineageD["Lineage"].split(';')[2]
						print '\t'.join([line[0], line[1], level])
						continue
				else:
					print '\t'.join([line[0], line[1], kingdom])
					
if __name__ == '__main__':
	infile,euk = sys.argv[1],sys.argv[2]
	get_taxonomic_level(infile,euk)