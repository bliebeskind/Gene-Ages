import sys

infile = sys.argv[1]

with open(infile) as f:
    count = 0
    for line in f:
        line = line.strip()
	if line == '':
	    continue
        line = line.split("\t")
        D = eval(line[2])
        lineage = "\t".join([i.strip() for i in D["Lineage"].split(";")])
        print "\t".join([line[0].strip(),line[1].strip(),lineage])
	count += 1
        if count % 100 == 0:
            sys.stderr.write(str(count)+"\n")
