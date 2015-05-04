import sys

def skip_header(infile):
	with open(infile) as f:
		skip = True
		for line in f:
			if skip:
				if line.startswith("AADNV"):
					skip = False
					yield line
				else:
					continue
			else:
				yield line

def parse_spec(src):
	cols = ["code","kingdom","taxonID","official name","common name","synonym"]
	line = []
	for i in src:
		if ":" in i:
			i = i.strip().split()
			i[2] = i[2].strip(":")
			if line != []:
				yield dict(map(None,cols,line))
				line = []
			line += i[:3]
			line += [' '.join(i[3:])[2:]]
		else:
			line.append(i[i.find("=")+1:].strip())
	yield dict(map(None,cols,line))
	
def print_csv(infile):
	gen = skip_header(infile)
	cols = ["code","kingdom","taxonID","official name","common name","synonym"]
	is_first = True
	for line in parse_spec(gen):
		if is_first:
			print ",".join(cols)
			is_first = False
		print ",".join([str(line[i]) for i in cols])
		
if __name__ == '__main__':
	print_csv(sys.argv[1])
			
		
		
