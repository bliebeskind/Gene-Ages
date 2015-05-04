# coding: utf-8
def clean_mim(infile):
    gene = []
    with open(infile) as f:
        for line in f:
            if line.startswith("ENSG") and len(gene) >0:
                yield ''.join(gene)+"\n"
                gene = []
            gene.append(line.strip())
            
